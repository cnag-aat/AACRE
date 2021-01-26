"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module contains functions relating to Pilon polishing of a Unicycler assembly.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import os
import subprocess
import shutil
from collections import defaultdict
from .misc import load_fasta, reverse_complement, int_to_str, underline, get_percentile_sorted, dim
from .assembly_graph import AssemblyGraph
from .assembly_graph_segment import Segment
from .string_graph import StringGraph, StringGraphSegment
from . import settings
from . import log


class CannotPolish(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


def polish_with_pilon_multiple_rounds(graph, insert_size_graph, args, polish_dir,
                                      do_pilon_reassembly):
    """
    This function does multiple rounds of Pilon polishing, until either it stops making changes
    or the limit is reached. It takes two graphs - one for polishing and one for insert size
    determination. They can be the same graph, but when polishing long read unitigs, it's quicker
    to get the insert size from a short read graph.
    """
    if not os.path.exists(polish_dir):
        os.makedirs(polish_dir)

    # TO DO: I'd like to come up with a way to polish (maybe with Pilon, maybe with something else)
    # that has a sensitivity adjustment. That way early rounds of polishing can take care of the
    # obvious errors and later rounds can make more tough calls. I hope this means that long read
    # problems (e.g. homopolymer runs) get fixed first and tougher problems (e.g. SNPs in repeat
    # regions) get fixed later, when the sequence is cleaner and easier to align to.

    # To avoid issues with paths that contain spaces, we will move into the temporary directory
    # to run these commands.
    starting_dir = os.getcwd()
    os.chdir(polish_dir)

    fix_type = 'bases'
    insert_size_1st, insert_size_99th = get_insert_size_range(insert_size_graph, args, polish_dir)
    for i in range(settings.MAX_PILON_POLISH_COUNT):
        if i > 0:
            graph.rotate_circular_sequences()
        change_count = polish_with_pilon(graph, args, polish_dir, insert_size_1st,
                                         insert_size_99th, i+1, fix_type)
        if change_count == 0:
            if fix_type == 'bases' and do_pilon_reassembly:
                fix_type = 'all'
            else:
                break

    os.chdir(starting_dir)
    if args.keep < 3 and os.path.exists(polish_dir):
        shutil.rmtree(polish_dir, ignore_errors=True)

    return insert_size_1st, insert_size_99th


def get_insert_size_range(graph, args, polish_dir):
    """
    This function just does a quick alignment of the paired-end reads to figure out the 1st and
    99th percentiles for the insert size.
    """
    using_paired_reads = bool(args.short1) and bool(args.short2)
    if not using_paired_reads:
        return 0, 1000

    log.log('Aligning reads to find appropriate insert size range...')

    segments_to_polish = [x for x in graph.segments.values()
                          if x.get_length() >= args.min_polish_size]
    if not segments_to_polish:
        raise CannotPolish('segments are too short')

    fasta_filename = '0_insert_size_check.fasta'
    sam_filename = '0_alignments.sam'

    with open(fasta_filename, 'w') as polish_fasta:
        for segment in segments_to_polish:
            polish_fasta.write('>' + get_segment_name(segment) + '\n')
            polish_fasta.write(segment.forward_sequence)
            polish_fasta.write('\n')

    # Prepare the FASTA for Bowtie2 alignment.
    bowtie2_build_command = [args.bowtie2_build_path, fasta_filename, fasta_filename]
    log.log(dim('  ' + ' '.join(bowtie2_build_command)), 2)
    try:
        subprocess.check_output(bowtie2_build_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('bowtie2-build encountered an error:\n' + e.output.decode())
    if not any(x.endswith('.bt2') for x in os.listdir(polish_dir)):
        raise CannotPolish('bowtie2-build failed to build an index')

    # Perform the alignment with Bowtie2.
    bowtie2_command = [args.bowtie2_path, '-1', args.short1, '-2', args.short2,
                       '-x', fasta_filename, '--fast',
                       '--threads', str(args.threads), '-I', '0', '-X', '5000', '-S', sam_filename]
    log.log(dim('  ' + ' '.join(bowtie2_command)), 2)
    try:
        subprocess.check_output(bowtie2_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('Bowtie2 encountered an error:\n' + e.output.decode())

    insert_sizes = []
    with open(sam_filename, 'rt') as raw_sam:
        for sam_line in raw_sam:
            try:
                sam_parts = sam_line.split('\t')
                sam_flags = int(sam_parts[1])
                if sam_flags & 2:  # if read mapped in proper pair
                    insert_size = int(sam_parts[8])
                    if insert_size > 0.0:
                        insert_sizes.append(insert_size)
            except (ValueError, IndexError):
                pass
    if not insert_sizes:
        raise CannotPolish('no read pairs aligned')
    insert_sizes = sorted(insert_sizes)
    insert_size_1st = get_percentile_sorted(insert_sizes, 1.0)
    insert_size_99th = get_percentile_sorted(insert_sizes, 99.0)

    log.log('Insert size 1st percentile:  ' + str(insert_size_1st))
    log.log('Insert size 99th percentile: ' + str(insert_size_99th))
    log.log('')

    if args.keep < 3:
        for f in [fasta_filename, sam_filename, fasta_filename + '.1.bt2',
                  fasta_filename + '.2.bt2', fasta_filename + '.3.bt2', fasta_filename + '.4.bt2',
                  fasta_filename + '.rev.1.bt2', fasta_filename + '.rev.2.bt2']:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass

    return insert_size_1st, insert_size_99th


def run_bowtie_samtools_commands(args, bowtie2_command, sam_filename, bam_filename):

    # Run bowtie2 alignment command
    log.log(dim('  ' + ' '.join(bowtie2_command)), 2)
    try:
        subprocess.check_output(bowtie2_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('Bowtie2 encountered an error:\n' + e.output.decode())

    # Sort the alignments.
    samtools_sort_command = [args.samtools_path, 'sort', '-@', str(args.threads),
                             '-o', bam_filename, '-O', 'bam', '-T', 'temp', sam_filename]
    log.log(dim('  ' + ' '.join(samtools_sort_command)), 2)
    try:
        subprocess.check_output(samtools_sort_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('Samtools encountered an error:\n' + e.output.decode())

    # Always delete the SAM file, regardless of keep level (it's big).
    try:
        os.remove(sam_filename)
    except FileNotFoundError:
        pass

    # Index the alignments.
    samtools_index_command = [args.samtools_path, 'index', bam_filename]
    log.log(dim('  ' + ' '.join(samtools_index_command)), 2)
    try:
        subprocess.check_output(samtools_index_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('Samtools encountered an error:\n' + e.output.decode())


def polish_with_pilon(graph, args, polish_dir, insert_size_1st, insert_size_99th, round_num,
                      fix_type):
    """
    Runs Pilon on the graph to hopefully fix up small mistakes.
    """
    log.log(underline('Pilon polish round ' + str(round_num)))

    using_paired_reads = bool(args.short1) and bool(args.short2)
    using_unpaired_reads = bool(args.unpaired)
    assert using_paired_reads or using_unpaired_reads

    input_filename = str(round_num) + '_polish_input.fasta'
    output_prefix = str(round_num) + '_pilon'
    pilon_fasta_filename = str(round_num) + '_pilon.fasta'
    pilon_changes_filename = str(round_num) + '_pilon.changes'
    pilon_output_filename = str(round_num) + '_pilon.out'

    segments_to_polish = [x for x in graph.segments.values()
                          if x.get_length() >= args.min_polish_size]
    if not segments_to_polish:
        raise CannotPolish('no segments are long enough to polish')

    with open(input_filename, 'w') as polish_fasta:
        for segment in segments_to_polish:
            polish_fasta.write('>' + get_segment_name(segment) + '\n')
            polish_fasta.write(segment.forward_sequence)
            polish_fasta.write('\n')

    # Prepare the FASTA for Bowtie2 alignment.
    bowtie2_build_command = [args.bowtie2_build_path, input_filename, input_filename]
    log.log(dim('  ' + ' '.join(bowtie2_build_command)), 2)
    try:
        subprocess.check_output(bowtie2_build_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('bowtie2-build encountered an error:\n' + e.output.decode())
    if not any(x.endswith('.bt2') for x in os.listdir(polish_dir)):
        raise CannotPolish('bowtie2-build failed to build an index')

    # Perform the alignment with Bowtie2.
    bowtie2_command = [args.bowtie2_path, '--local', '--very-sensitive-local',
                       '--threads', str(args.threads), '-I', str(insert_size_1st),
                       '-X', str(insert_size_99th), '-x', input_filename]
    if using_paired_reads:
        paired_sam_filename = str(round_num) + '_paired_alignments.sam'
        paired_bam_filename = str(round_num) + '_paired_alignments.bam'
        this_bowtie2_command = bowtie2_command + ['-S', paired_sam_filename,
                                                  '-1', args.short1, '-2', args.short2]
        run_bowtie_samtools_commands(args, this_bowtie2_command, paired_sam_filename,
                                     paired_bam_filename)
    else:
        paired_bam_filename = ''

    if using_unpaired_reads:
        unpaired_sam_filename = str(round_num) + '_unpaired_alignments.sam'
        unpaired_bam_filename = str(round_num) + '_unpaired_alignments.bam'
        this_bowtie2_command = bowtie2_command + ['-S', unpaired_sam_filename, '-U', args.unpaired]
        run_bowtie_samtools_commands(args, this_bowtie2_command, unpaired_sam_filename,
                                     unpaired_bam_filename)
    else:
        unpaired_bam_filename = ''

    # Polish with Pilon.
    if args.pilon_path.endswith('.jar'):
        pilon_command = [args.java_path, '-jar', args.pilon_path]
    else:
        pilon_command = [args.pilon_path]
    pilon_command += ['--genome', input_filename, '--changes',
                      '--output', output_prefix, '--outdir', polish_dir, '--fix', fix_type]
    if using_paired_reads:
        pilon_command += ['--frags', paired_bam_filename]
    if using_unpaired_reads:
        pilon_command += ['--unpaired', unpaired_bam_filename]

    log.log(dim('  ' + ' '.join(pilon_command)), 2)
    try:
        pilon_stdout = subprocess.check_output(pilon_command, stderr=subprocess.STDOUT)
        with open(pilon_output_filename, 'wb') as pilon_out:
            pilon_out.write(pilon_stdout)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('Pilon encountered an error:\n' + e.output.decode())
    if not os.path.isfile(pilon_fasta_filename):
        raise CannotPolish('Pilon did not produce FASTA file')
    if not os.path.isfile(pilon_changes_filename):
        raise CannotPolish('Pilon did not produce changes file')

    # Display Pilon changes.
    change_count = defaultdict(int)
    change_lines = defaultdict(list)
    total_count = 0
    pilon_changes = open(pilon_changes_filename, 'rt')
    for line in pilon_changes:
        try:
            seg_name = line.split(':')[0]
            change_count[seg_name] += 1
            total_count += 1
            change_lines[seg_name].append(line.strip())
        except ValueError:
            pass
    if total_count == 0:
        log.log('No Pilon changes')
    else:
        log.log('Total number of changes: ' + int_to_str(total_count))
        log.log('', 2)
        seg_names = sorted(graph.segments)
        polish_input_seg_names = set(get_segment_name_or_number(x) for x in segments_to_polish)
        for seg_name in seg_names:
            if seg_name in polish_input_seg_names:
                count = change_count[str(seg_name)]
                if count < 1:
                    continue
                log.log('Segment ' + str(seg_name) + ' (' +
                        int_to_str(graph.segments[seg_name].get_length()) + ' bp): ' +
                        int_to_str(count) + ' change' +
                        ('s' if count > 1 else ''), 2)
                try:
                    changes = change_lines[str(seg_name)]
                    changes = sorted(changes, key=lambda x:
                                     int(x.replace(' ', ':').replace('-', ':').split(':')[1]))
                    for change in changes:
                        log.log('  ' + dim(change), 3)
                except (ValueError, IndexError):
                    pass

        # Replace segment sequences with Pilon-polished versions.
        pilon_results = load_fasta(pilon_fasta_filename)
        for header, sequence in pilon_results:
            if header.endswith('_pilon'):
                header = header[:-6]
            if isinstance(graph, AssemblyGraph):
                segment = graph.segments[int(header)]
            elif isinstance(graph, StringGraph):
                segment = graph.segments[header]
            else:
                assert False
            segment.forward_sequence = sequence
            segment.reverse_sequence = reverse_complement(sequence)

    log.log('')

    if args.keep < 3:
        list_of_files = [input_filename, pilon_fasta_filename, pilon_changes_filename,
                         input_filename + '.1.bt2', input_filename + '.2.bt2',
                         input_filename + '.3.bt2', input_filename + '.4.bt2',
                         input_filename + '.rev.1.bt2', input_filename + '.rev.2.bt2',
                         pilon_output_filename, paired_bam_filename, paired_bam_filename + '.bai',
                         unpaired_bam_filename, unpaired_bam_filename + '.bai']
        for f in list_of_files:
            try:
                os.remove(f)
            except (FileNotFoundError, OSError):
                pass

    return total_count


def get_segment_name(segment):
    if isinstance(segment, Segment):
        return str(segment.number)
    elif isinstance(segment, StringGraphSegment):
        return segment.full_name
    else:
        assert False


def get_segment_name_or_number(segment):
    if isinstance(segment, Segment):
        return segment.number
    elif isinstance(segment, StringGraphSegment):
        return segment.full_name
    else:
        assert False
