"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module contains functionality related to miniasm, which Unicycler uses to build an assembly
using both Illumina contigs and long reads.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import os
import shutil
import subprocess
import sys
import itertools
import collections
from .misc import green, red, line_iterator, print_table, int_to_str, float_to_str, \
    reverse_complement, gfa_path, racon_version
from .minimap_alignment import align_long_reads_to_assembly_graph, range_overlap_size, \
    load_minimap_alignments
from .string_graph import StringGraph, StringGraphSegment, \
    merge_string_graph_segments_into_unitig_graph
from .read_ref import load_references, load_long_reads
from .unicycler_align import semi_global_align_long_reads
from . import log
from . import settings

try:
    from .cpp_wrappers import minimap_align_reads, miniasm_assembly, start_seq_alignment, \
        end_seq_alignment
except AttributeError as att_err:
    sys.exit('Error when importing C++ library: ' + str(att_err) + '\n'
             'Have you successfully built the library file using make?')


class MiniasmFailure(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


def make_miniasm_string_graph(graph, read_dict, long_read_filename, scoring_scheme, read_nicknames,
                              counter, args, anchor_segments, existing_long_read_assembly):
    log.log_section_header('Assembling contigs and long reads with miniasm')
    if graph is not None:
        log.log_explanation('Unicycler uses miniasm to construct a string graph assembly using '
                            'both the short read contigs and the long reads. It will then use the '
                            'resulting string graph to produce bridges between contigs. This '
                            'method requires decent coverage of long reads and therefore may not '
                            'be fruitful if long reads are sparse. However, it does not rely on '
                            'the short read assembly graph having good connectivity and is able to '
                            'bridge an assembly graph even when it contains many dead ends.',
                            extra_empty_lines_after=0)
        log.log_explanation('Unicycler uses two types of "reads" as assembly input: anchor contigs '
                            'from the short-read assembly and actual long reads which '
                            'overlap two or more of these contigs. It then assembles them with '
                            'miniasm.')
        if existing_long_read_assembly:
            log.log_explanation('Since you directly provided a long read assembly (using the '
                                '--existing_long_read_assembly option), Unicycler will use that '
                                'instead of performing its own miniasm/Racon assembly. However, '
                                'Unicycler still runs miniasm as this can provide useful contig '
                                'trimming information.')

    miniasm_dir = os.path.join(args.out, 'miniasm_assembly')
    if not os.path.exists(miniasm_dir):
        os.makedirs(miniasm_dir)

    short_reads_available = graph is not None
    seg_nums_to_bridge = set(x.number for x in anchor_segments)

    assembly_reads_filename = os.path.join(miniasm_dir, '01_assembly_reads.fastq')
    mappings_filename = os.path.join(miniasm_dir, '02_mappings.paf')
    string_graph_filename = os.path.join(miniasm_dir, '10_final_string_graph.gfa')
    branching_paths_removed_filename = os.path.join(miniasm_dir, '11_branching_paths_removed.gfa')
    unitig_graph_filename = os.path.join(miniasm_dir, '12_unitig_graph.gfa')
    racon_polished_filename = os.path.join(miniasm_dir, '13_racon_polished.gfa')
    pilon_polished_filename = os.path.join(miniasm_dir, '14_pilon_polished.gfa')
    contigs_placed_filename = os.path.join(miniasm_dir, '15_contigs_placed.gfa')
    miniasm_read_list = os.path.join(miniasm_dir, 'all_reads.txt')

    # If the long read assembly already exists, then we can skip ahead quite a lot.
    if os.path.isfile(pilon_polished_filename) and os.path.isfile(miniasm_read_list):
        log.log('Long read assembly already exists: ' + pilon_polished_filename)
        log.log('')
        unitig_graph = StringGraph(pilon_polished_filename)

    # If not, we need to do all of the long read assembly steps now.
    else:
        assembly_read_names = get_miniasm_assembly_reads(graph, read_dict, long_read_filename,
                                                         miniasm_dir, args.threads)

        # TO DO: identify chimeric reads and throw them out. This was part of miniasm, but it was
        # removed due to 'not working as intended', so I pulled it out of my miniasm as well.

        save_assembly_reads_to_file(assembly_reads_filename, assembly_read_names, read_dict, graph,
                                    seg_nums_to_bridge)

        # Do an all-vs-all alignment of the assembly FASTQ, for miniasm input. Contig-contig
        # alignments are excluded (because single-copy contigs, by definition, should not
        # significantly overlap each other).
        log.log('Finding overlaps with minimap... ', end='')
        minimap_alignments_str = minimap_align_reads(assembly_reads_filename,
                                                     assembly_reads_filename,
                                                     args.threads, 0, 'read vs read')
        overlap_count = 0
        with open(mappings_filename, 'wt') as mappings:
            for minimap_alignment_str in line_iterator(minimap_alignments_str):
                if minimap_alignment_str.count('CONTIG_') < 2:
                    mappings.write(minimap_alignment_str)
                    mappings.write('\n')
                    overlap_count += 1

        if overlap_count == 0:
            log.log(red('failed'))
        else:
            log.log(green('success'))
            log.log('  ' + int_to_str(overlap_count) + ' overlaps\n')

        # TO DO: refine these overlaps? Perhaps using Unicycler-align? I suspect that the quality
        # of a miniasm assembly is highly dependent on the input overlaps.
        #
        # I could even try to do more sophisticated stuff, like using the alignments to identify
        # repeat regions, then finding these repeat regions in reads and tossing out alignments
        # which are contained only in repeat regions.

        # TO DO: Unicycler's mode (conservative, normal or bold) should appropriately affect the
        # miniasm settings.

        # Now actually do the miniasm assembly, which will create a GFA file of the string graph.
        log.log('Assembling reads with miniasm... ', end='')
        min_depth = 3
        miniasm_assembly(assembly_reads_filename, mappings_filename, miniasm_dir, min_depth)
        if not os.path.isfile(string_graph_filename):
            log.log(red('failed'))
            raise MiniasmFailure('miniasm failed to generate a string graph')
        string_graph = StringGraph(string_graph_filename)

        if len(string_graph.segments) == 0:
            log.log(red('empty result'))
            unitig_graph = None

        else:
            log.log(green('success'))
            log.log('  ' + int_to_str(len(string_graph.segments)) + ' segments, ' +
                    int_to_str(len(string_graph.links) // 2) + ' links\n')

            if not short_reads_available and args.keep > 0:
                shutil.copyfile(string_graph_filename,
                                gfa_path(args.out, next(counter), 'string_graph'))

            string_graph.remove_branching_paths()
            string_graph.save_to_gfa(branching_paths_removed_filename, include_depth=False)

            log.log('Merging segments into unitigs:')
            unitig_graph = merge_string_graph_segments_into_unitig_graph(string_graph,
                                                                         read_nicknames)
            circular_count = unitig_graph.get_circular_segment_count()
            linear_count = unitig_graph.get_linear_segment_count()
            unitig_graph_size = unitig_graph.get_total_segment_length()
            if circular_count > 0:
                log.log('  ' + int_to_str(circular_count) + ' circular unitig' +
                        ('' if circular_count == 1 else 's'))
            if linear_count > 0:
                log.log('  ' + int_to_str(linear_count) + ' linear unitig' +
                        ('' if linear_count == 1 else 's'))
            log.log('  total size = ' + int_to_str(unitig_graph_size) + ' bp')

            unitig_graph.save_to_gfa(unitig_graph_filename, include_depth=False)
            if not short_reads_available and args.keep > 0:
                unitig_graph.save_to_gfa(gfa_path(args.out, next(counter), 'unitig_graph'),
                                         include_depth=False)

            # If the miniasm assembly looks too small, then we don't bother polishing it or using
            # it for bridging.
            if short_reads_available:
                estimated_genome_size = graph.get_estimated_sequence_len()
                target_size = estimated_genome_size * \
                    settings.REQUIRED_MINIASM_ASSEMBLY_SIZE_FOR_BRIDGING
                if unitig_graph_size < target_size:
                    log.log('')
                    log.log(red('miniasm assembly too small for bridging'))
                    unitig_graph = None

            if unitig_graph is not None:
                # If the user supplied a long read assembly, then we use that instead of doing a
                # Racon polish.
                if existing_long_read_assembly:
                    log.log('')
                    log.log('Using provided long read assembly instead of running Racon: ' +
                            existing_long_read_assembly)
                    log.log('')
                    unitig_graph = StringGraph(existing_long_read_assembly)
                else:
                    polish_unitigs_with_racon(unitig_graph, miniasm_dir, read_dict, graph,
                                              args.racon_path, args.threads, scoring_scheme,
                                              seg_nums_to_bridge)
                    unitig_graph.save_to_gfa(racon_polished_filename)
                    if not short_reads_available and args.keep > 0:
                        unitig_graph.save_to_gfa(gfa_path(args.out, next(counter),
                                                          'racon_polished'))
                if short_reads_available and args.keep > 0:
                    unitig_graph.save_to_gfa(gfa_path(args.out, next(counter),
                                                      'long_read_assembly'))

    if unitig_graph is not None and short_reads_available:
        log.log('')
        trim_dead_ends_based_on_miniasm_trimming(graph, miniasm_read_list)
        unitig_graph = place_contigs(miniasm_dir, graph, unitig_graph, args.threads,
                                     scoring_scheme, seg_nums_to_bridge)
        unitig_graph.save_to_gfa(contigs_placed_filename, include_depth=False)

    if args.keep < 3:
        shutil.rmtree(miniasm_dir, ignore_errors=True)
    return unitig_graph


def get_miniasm_assembly_reads(graph, read_dict, long_read_filename, miniasm_dir, threads):
    if graph is not None:  # hybrid assembly
        minimap_alignments = align_long_reads_to_assembly_graph(graph, long_read_filename,
                                                                miniasm_dir, threads)
        miniasm_assembly_reads = []
        for read_name, alignments in minimap_alignments.items():
            if any(a.overlaps_reference() for a in alignments):
                miniasm_assembly_reads.append(read_name)
        return sorted(miniasm_assembly_reads)
    else:  # long-read-only assembly
        return sorted([x for x in read_dict.keys()])


def save_assembly_reads_to_file(read_filename, read_names, read_dict, graph, seg_nums_to_bridge,
                                contig_copy_count=1):
    qual = chr(settings.CONTIG_READ_QSCORE + 33)
    log.log('Saving to ' + read_filename + ':')

    with open(read_filename, 'wt') as fastq:
        if graph is not None:  # hybrid assembly
            # First save the Illumina contigs as 'reads'. They are given a constant high qscore to
            # reflect our confidence in them.
            seg_count = 0
            for seg in sorted(graph.segments.values(), key=lambda x: x.number):
                if segment_suitable_for_miniasm_assembly(graph, seg, seg_nums_to_bridge):
                    for i in range(contig_copy_count):
                        fastq_name = '@CONTIG_' + str(seg.number)
                        if contig_copy_count > 1:
                            fastq_name += '_' + str(i+1)
                        fastq.write(fastq_name + '\n')
                        if i % 2 == 0:  # evens
                            fastq.write(seg.forward_sequence)
                        else:  # odds
                            fastq.write(seg.reverse_sequence)
                        fastq.write('\n+\n')
                        fastq.write(qual * seg.get_length())
                        fastq.write('\n')
                    seg_count += 1
            if contig_copy_count > 0:
                message = '  ' + int_to_str(seg_count) + ' short-read contigs'
                if contig_copy_count > 1:
                    message += ' (duplicated ' + str(contig_copy_count) + ' times)'
                log.log(message)

        # Now save the actual long reads.
        for read_name in read_names:
            read = read_dict[read_name]
            seq = read.sequence
            if len(seq) < 100:
                continue
            quals = read.qualities
            fastq.write('@' + read_name + '\n')
            fastq.write(seq)
            fastq.write('\n+\n')
            fastq.write(quals)
            fastq.write('\n')
        log.log('  ' + int_to_str(len(read_names)) + ' long reads')
        log.log('')


def segment_suitable_for_miniasm_assembly(graph, segment, seg_nums_to_bridge):
    """
    Returns True if the segment is:
      1) one of the to-be-bridged segments
      2) not already circular and complete
    """
    if segment.number not in seg_nums_to_bridge:
        return False
    return not graph.is_component_complete([segment.number])


def polish_unitigs_with_racon(unitig_graph, miniasm_dir, read_dict, graph, racon_path, threads,
                              scoring_scheme, seg_nums_to_bridge):
    log.log_section_header('Polishing miniasm assembly with Racon')
    log.log_explanation('Unicycler now uses Racon to polish the miniasm assembly. It does '
                        'multiple rounds of polishing to get the best consensus. Circular unitigs '
                        'are rotated between rounds such that all parts (including the ends) are '
                        'polished well.')

    polish_dir = os.path.join(miniasm_dir, 'racon_polish')
    if not os.path.isdir(polish_dir):
        os.makedirs(polish_dir)

    # TO DO: I might want to filter the reads used for polishing. In particular, throw out
    # chimeric reads (if I come up with a good way of spotting them) and reads with a window that
    # drops below a quality threshold.

    polish_read_names = sorted([x for x in read_dict.keys()])
    polish_reads = os.path.join(polish_dir, 'polishing_reads.fastq')
    save_assembly_reads_to_file(polish_reads, polish_read_names, read_dict, graph,
                                seg_nums_to_bridge, settings.RACON_CONTIG_DUPLICATION_COUNT)

    col_widths = [6, 12, 14]
    racon_table_header = ['Polish round', 'Assembly size', 'Mapping quality']
    print_table([racon_table_header], fixed_col_widths=col_widths, left_align_header=False,
                alignments='LRR', indent=0)

    best_fasta = None
    best_unitig_sequences = {}
    best_mapping_quality = 0
    times_quality_failed_to_beat_best = 0

    counter = itertools.count(start=1)
    current_fasta = os.path.join(polish_dir, ('%03d' % next(counter)) + '_unpolished_unitigs.fasta')
    unitig_graph.save_to_fasta(current_fasta)

    if graph is None:  # Long-read-only assembly
        racon_loop_count = settings.RACON_POLISH_LOOP_COUNT_LONG_ONLY
    else:  # Hybrid assembly
        racon_loop_count = settings.RACON_POLISH_LOOP_COUNT_HYBRID

    # The Racon command will be different for older versions of Racon.
    old_racon_version = (racon_version(racon_path) == '-')

    for polish_round_count in range(racon_loop_count):

        mappings_filename = os.path.join(polish_dir, ('%03d' % next(counter)) + '_alignments.paf')
        racon_log = os.path.join(polish_dir, ('%03d' % next(counter)) + '_racon.log')
        polished_fasta = os.path.join(polish_dir, ('%03d' % next(counter)) + '_polished.fasta')
        fixed_fasta = os.path.join(polish_dir, ('%03d' % next(counter)) + '_fixed.fasta')
        rotated_fasta = os.path.join(polish_dir, ('%03d' % next(counter)) + '_rotated.fasta')

        mapping_quality, unitig_depths = \
            make_racon_polish_alignments(current_fasta, mappings_filename, polish_reads, threads)

        racon_table_row = ['begin' if polish_round_count == 0 else str(polish_round_count),
                           int_to_str(unitig_graph.get_total_segment_length()),
                           float_to_str(mapping_quality, 2)]
        print_table([racon_table_row], fixed_col_widths=col_widths, left_align_header=False,
                    alignments='LRR', indent=0, header_format='normal', bottom_align_header=False)

        # Do we have a new best?
        if mapping_quality > best_mapping_quality:
            best_mapping_quality = mapping_quality
            best_fasta = current_fasta
            best_unitig_sequences = {name: seg.forward_sequence
                                     for name, seg in unitig_graph.segments.items()}
            times_quality_failed_to_beat_best = 0
        else:
            times_quality_failed_to_beat_best += 1

        # If we've failed to improve on our best quality for a few rounds, then we're done!
        if times_quality_failed_to_beat_best > 2:
            break

        # Run Racon. It crashes sometimes, so repeat until its return code is 0.
        return_code = 1
        for t in range(100):  # Only try a fixed number of times, to prevent an infinite loop.

            # The old version of Racon takes the output file (polished fasta) as an argument.
            if old_racon_version:
                command = [racon_path, '--verbose', '9', '-t', str(threads), '--bq', '-1',
                           polish_reads, mappings_filename, current_fasta, polished_fasta]

            # The new version of Racon outputs the polished fasta to stdout.
            else:
                command = [racon_path, '-t', str(threads), polish_reads, mappings_filename,
                           current_fasta]

            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = process.communicate()
            return_code = process.returncode

            if old_racon_version:
                with open(racon_log, 'wb') as log_file:
                    log_file.write(out)
                    log_file.write(err)
            else:
                with open(racon_log, 'wb') as log_file:
                    log_file.write(err)
                with open(polished_fasta, 'wb') as out_file:
                    out_file.write(out)

            if return_code == 0 and os.path.isfile(polished_fasta):
                break
            if os.path.isfile(polished_fasta):
                os.remove(polished_fasta)
            if os.path.isfile(racon_log):
                os.remove(racon_log)

        # If even after all those tries Racon still didn't succeed, then we give up!
        if return_code != 0 or not os.path.isfile(polished_fasta):
            break

        unitig_graph.replace_with_polished_sequences(polished_fasta, scoring_scheme,
                                                     old_racon_version)
        unitig_graph.save_to_fasta(fixed_fasta)
        unitig_graph.rotate_circular_sequences()
        unitig_graph.save_to_fasta(rotated_fasta)
        current_fasta = rotated_fasta

    log.log('')
    if best_fasta:
        log.log('Best polish: ' + best_fasta)
        for unitig_name, unitig_seq in best_unitig_sequences.items():
            segment = unitig_graph.segments[unitig_name]
            segment.forward_sequence = unitig_seq
            segment.reverse_sequence = reverse_complement(unitig_seq)
            if unitig_name in unitig_depths:
                segment.depth = unitig_depths[unitig_name]
        unitig_graph.normalise_read_depths()
    else:
        log.log(red('Polishing failed'))


def place_contigs(miniasm_dir, assembly_graph, unitig_graph, threads, scoring_scheme,
                  seg_nums_to_bridge):
    log.log('', verbosity=1)
    log.log_explanation('Unicycler now places the single copy contigs back into the unitig '
                        'graph. This serves two purposes: a) it replaces long read assembly '
                        'sequences (which may be error prone) with Illumina assembly sequence '
                        '(which is probably quite accurate), improving the assembly quality, and '
                        'b) it defines inter-contig sequences for use in building bridges.',
                        verbosity=1)

    # Use semi-global alignment of contig ends to find the contigs in the unitig graph. We first
    # look for contigs using relatively large chunks of their ends and progress to smaller chunks
    # only if they aren't found.
    contig_positions = []
    contig_numbers = []
    for seg in sorted(assembly_graph.segments.values(), key=lambda x: x.number):
        if segment_suitable_for_miniasm_assembly(assembly_graph, seg, seg_nums_to_bridge):
            contig_numbers.append(seg.number)
    for contig_search_end_size in settings.CONTIG_SEARCH_END_SIZES:
        log.log('Searching for contigs using ' + int_to_str(contig_search_end_size) + ' bp of '
                'contig ends.')
        log.log('')

        contig_position_results, not_found_contig_numbers = \
            find_contig_starts_and_ends(miniasm_dir, assembly_graph, unitig_graph, threads,
                                        scoring_scheme, contig_search_end_size, contig_numbers)
        contig_positions += contig_position_results
        if not_found_contig_numbers:
            contig_numbers = not_found_contig_numbers
        else:
            break

    # Now we build a new string graph that consists of contigs and their bridging sequences!
    new_graph = StringGraph(None)
    bridge_num = itertools.count(start=1)
    for seg in sorted(unitig_graph.segments.values(), key=lambda x: x.get_length(), reverse=True):

        # We work through each contig separately.
        unitig_name = seg.full_name
        unitig_seq = seg.forward_sequence
        extended_unitig_seq = unitig_seq + unitig_seq
        unitig_length = len(seg.forward_sequence)
        circular_unitig = unitig_graph.segment_is_circular(unitig_name)

        # Do any contigs overlap each other a lot? It's very weird if they do, so throw them out!
        good_contig_positions = []
        unitig_contig_positions = [x for x in contig_positions if x[3] == seg.full_name]
        for i, contig_position in enumerate(unitig_contig_positions):
            other_positions = unitig_contig_positions[:i] + unitig_contig_positions[i+1:]
            if range_overlap_size(contig_position[:2], [x[:2] for x in other_positions]) <= \
                    settings.FOUND_CONTIG_MAX_OVERLAP_SIZE:
                good_contig_positions.append(contig_position)
        unitig_contig_positions = sorted(good_contig_positions)

        # In the unfortunate scenario of no contig placements, the unitig just goes into the new
        # graph as it is.
        segment_names = []
        if not unitig_contig_positions:
            seg_name = 'BRIDGE_' + str(next(bridge_num))
            new_graph.segments[seg_name] = StringGraphSegment(seg_name, unitig_seq)
            segment_names.append(seg_name + '+')

        # In the much more ideal scenario where there is at least one contig placement, we make
        # the contig and bridge segments in order of their position on the unitig.
        for i, contig_position in enumerate(unitig_contig_positions):
            start_pos, end_pos, rev_strand, _, contig_number = contig_position

            # If this is the first contig and the unitig is linear, make a 'bridge' for the start.
            if i == 0 and not circular_unitig:
                assert start_pos >= 0
                bridge_seq = unitig_seq[:start_pos]
                if bridge_seq:
                    seg_name = 'BRIDGE_' + str(next(bridge_num))
                    new_graph.segments[seg_name] = StringGraphSegment(seg_name, bridge_seq)
                    segment_names.append(seg_name + '+')

            # Make the contig segment.
            seg_name = 'CONTIG_' + str(contig_number)
            contig_seq = assembly_graph.segments[contig_number].forward_sequence
            new_graph.segments[seg_name] = StringGraphSegment(seg_name, contig_seq)
            if rev_strand:
                segment_names.append(seg_name + '-')
            else:
                segment_names.append(seg_name + '+')

            # Make a bridge to the next segment.
            not_last = i < len(unitig_contig_positions) - 1
            last_and_circular = (i == len(unitig_contig_positions) - 1 and circular_unitig)
            if not_last or last_and_circular:
                assert end_pos >= 0
                bridge_start = end_pos
                if not_last:
                    bridge_end = unitig_contig_positions[i+1][0]
                else:   # last_and_circular
                    bridge_end = unitig_contig_positions[0][0] + unitig_length

                if bridge_end >= bridge_start:  # Non-overlapping bridges (length of 0 or more)
                    bridge_seq = extended_unitig_seq[bridge_start:bridge_end]
                    seg_name = 'BRIDGE_' + str(next(bridge_num))
                else:  # Overlapping bridges
                    bridge_seq = extended_unitig_seq[bridge_end:bridge_start]
                    seg_name = 'OVERLAPPING_BRIDGE_' + str(next(bridge_num))
                new_graph.segments[seg_name] = StringGraphSegment(seg_name, bridge_seq)
                segment_names.append(seg_name + '+')

            # If this is the last contig and it's linear, then make a bridge to the unitig's end.
            if i == len(unitig_contig_positions) - 1 and not circular_unitig:
                bridge_seq = unitig_seq[end_pos:unitig_length]
                if bridge_seq:
                    seg_name = 'BRIDGE_' + str(next(bridge_num))
                    new_graph.segments[seg_name] = StringGraphSegment(seg_name, bridge_seq)
                    segment_names.append(seg_name + '+')

        # Create links between all the new segments.
        if circular_unitig:
            segment_names.append(segment_names[0])
        for i, seg_2 in enumerate(segment_names[1:]):
            seg_1 = segment_names[i]

            seg_1_overlapping_bridge = seg_1.startswith('OVERLAPPING_BRIDGE')
            seg_2_overlapping_bridge = seg_2.startswith('OVERLAPPING_BRIDGE')
            assert not (seg_1_overlapping_bridge and seg_2_overlapping_bridge)

            seg_1_seq = new_graph.seq_from_signed_seg_name(seg_1)
            seg_2_seq = new_graph.seq_from_signed_seg_name(seg_2)

            if seg_1_overlapping_bridge:
                seg_1_overlap = len(seg_1_seq)
                seg_2_overlap = start_seq_alignment(seg_1_seq, seg_2_seq, scoring_scheme)
                new_graph.add_link(seg_1, seg_2, seg_1_overlap, seg_2_overlap)
            elif seg_2_overlapping_bridge:
                seg_1_overlap = len(seg_1_seq) - end_seq_alignment(seg_2_seq, seg_1_seq,
                                                                   scoring_scheme)
                seg_2_overlap = len(seg_2_seq)
                new_graph.add_link(seg_1, seg_2, seg_1_overlap, seg_2_overlap)
            else:  # neither is an overlapping bridge
                new_graph.add_link(seg_1, seg_2, 0, 0)

    return new_graph


def find_contig_starts_and_ends(miniasm_dir, assembly_graph, unitig_graph, threads, scoring_scheme,
                                contig_search_end_size, contig_numbers):
    """
    This function searches for contig start and end points in the unitig graph.
    """
    if not contig_numbers:
        return [], []

    contig_search_dir = os.path.join(miniasm_dir, 'contig_search')
    if not os.path.isdir(contig_search_dir):
        os.makedirs(contig_search_dir)

    # Save the contigs to a FASTA file.
    contigs_fasta = os.path.join(contig_search_dir, 'contigs.fasta')
    longest_contig_seq_len = 0
    smallest_contig_seq_len = float('inf')
    with open(contigs_fasta, 'wt') as fasta:
        for contig_number in contig_numbers:
            seg = assembly_graph.segments[contig_number]
            seq = seg.forward_sequence

            longest_contig_seq_len = max(longest_contig_seq_len, len(seq))
            smallest_contig_seq_len = min(smallest_contig_seq_len, len(seq))

            contig_name = 'CONTIG_' + str(seg.number)

            # If the contig is long, we just look for the ends.
            if len(seq) >= contig_search_end_size * 2:
                start_seq = seq[:contig_search_end_size]
                end_seq = seq[-contig_search_end_size:]
                fasta.write('>' + contig_name + '_START\n')
                fasta.write(start_seq)
                fasta.write('\n')
                fasta.write('>' + contig_name + '_END\n')
                fasta.write(end_seq)
                fasta.write('\n')

            # If the contig is short, we look for the whole thing.
            else:
                fasta.write('>' + contig_name + '_WHOLE\n')
                fasta.write(seq)
                fasta.write('\n')

    # Save the unitigs to a FASTA file. If the segment is circular, then we save it extra sequence
    # on the end, so we can find contigs which loop.
    unitigs_fasta = os.path.join(contig_search_dir, 'unitigs.fasta')
    with open(unitigs_fasta, 'wt') as fasta:
        for seg in sorted(unitig_graph.segments.values(), key=lambda x: x.get_length(),
                          reverse=True):
            seg_seq = seg.forward_sequence
            if unitig_graph.segment_is_circular(seg.full_name):
                if len(seg_seq) <= longest_contig_seq_len:
                    seg_seq += seg_seq
                else:
                    seg_seq += seg_seq[:longest_contig_seq_len]
            fasta.write('>' + seg.full_name + '\n')
            fasta.write(seg_seq)
            fasta.write('\n')

    # Now we align the contigs using the semi-global read aligning method. In this case, the
    # 'reads' are the contigs and the 'references' are the unitigs.
    references = load_references(unitigs_fasta, section_header=None, show_progress=False)
    read_dict, read_names, read_filename = load_long_reads(contigs_fasta, silent=True)
    min_alignment_len = min(contig_search_end_size * 0.9, smallest_contig_seq_len * 0.9)
    semi_global_align_long_reads(references, unitigs_fasta, read_dict, read_names, read_filename,
                                 threads, scoring_scheme, [None], False, min_alignment_len, None,
                                 None, 10, 0, None, verbosity=0)
    start_positions = {}
    end_positions = {}
    for contig_name in read_names:
        contig_number = int(contig_name.split('_')[1])
        contig = read_dict[contig_name]
        if contig.alignments:
            a = sorted(contig.alignments, key=lambda x: x.scaled_score)[-1]
            if a.percent_identity < settings.CONTIG_SEARCH_MIN_IDENTITY:
                continue

            unitig_name = a.ref.name
            if contig_name.endswith('_START') or contig_name.endswith('_WHOLE'):
                if a.rev_comp:  # reverse strand
                    start_positions[contig_number] = (unitig_name, a.ref_end_pos, a.rev_comp)
                else:  # forward strand
                    start_positions[contig_number] = (unitig_name, a.ref_start_pos, a.rev_comp)
            if contig_name.endswith('_END') or contig_name.endswith('_WHOLE'):
                if a.rev_comp:  # reverse strand
                    end_positions[contig_number] = (unitig_name, a.ref_start_pos, a.rev_comp)
                else:  # forward strand
                    end_positions[contig_number] = (unitig_name, a.ref_end_pos, a.rev_comp)

    # Gather up the start/end positions for all found contigs.
    contig_positions = []
    for contig_number in contig_numbers:
        if contig_number in start_positions and contig_number in end_positions:
            start_unitig_name, start_pos, start_rev_comp = start_positions[contig_number]
            end_unitig_name, end_pos, end_rev_comp = end_positions[contig_number]

            # Only accept start/ends on the same strand of the same unitig.
            if start_unitig_name != end_unitig_name:
                continue
            if start_rev_comp != end_rev_comp:
                continue

            unitig_name = start_unitig_name
            unitig_length = unitig_graph.segments[unitig_name].get_length()
            circular_unitig = unitig_graph.segment_is_circular(unitig_name)
            rev_comp = start_rev_comp

            # Make sure each position is in the actual unitig coordinates, not in the overlap past
            # the end.
            if start_pos >= unitig_length:
                start_pos -= unitig_length
            if end_pos >= unitig_length:
                end_pos -= unitig_length

            # Fix up start/ends that overlap the end of a circular unitig.
            if start_pos > end_pos and not rev_comp and circular_unitig:
                start_pos -= unitig_length
            if end_pos > start_pos and rev_comp and circular_unitig:
                end_pos -= unitig_length

            if rev_comp:
                start_pos, end_pos = end_pos, start_pos

            contig = assembly_graph.segments[contig_number]
            length_of_contig_in_unitig = end_pos - start_pos
            length_ratio = length_of_contig_in_unitig / contig.get_length()
            if length_ratio < settings.FOUND_CONTIG_MIN_RATIO or \
                    length_ratio > settings.FOUND_CONTIG_MAX_RATIO:
                continue
            contig_positions.append((start_pos, end_pos, rev_comp, unitig_name, contig_number))

    # Display the search results in a table.
    not_found_contig_numbers = []
    contig_search_table = [['Contig', 'Result', 'Start pos', 'End pos', 'Strand']]
    for contig_number in contig_numbers:
        result = [x for x in contig_positions if x[4] == contig_number]
        assert len(result) == 0 or len(result) == 1
        if result:
            result = result[0]
            contig_search_table.append([str(contig_number), 'found in unitig ' + result[3],
                                        str(result[0]), str(result[1]), '-' if result[2] else '+'])
        else:
            contig_search_table.append([str(contig_number), 'not found', '', '', ''])
            not_found_contig_numbers.append(contig_number)
    print_table(contig_search_table, alignments='RLRRR', indent=0, sub_colour={'not found': 'red'})
    log.log('')

    return contig_positions, not_found_contig_numbers


def make_racon_polish_alignments(current_fasta, mappings_filename, polish_reads, threads):
    mapping_quality = 0
    unitig_depths = collections.defaultdict(float)

    minimap_alignments_str = minimap_align_reads(current_fasta, polish_reads, threads, 3,
                                                 preset_name='find contigs')
    alignments_by_read = load_minimap_alignments(minimap_alignments_str,
                                                 filter_overlaps=True, allowed_overlap=10,
                                                 filter_by_minimisers=True)
    with open(mappings_filename, 'wt') as mappings:
        for read_name in sorted(alignments_by_read.keys()):
            for a in alignments_by_read[read_name]:
                mappings.write(a.paf_line)
                mappings.write('\n')
                mapping_quality += a.matching_bases / a.num_bases
                unitig_depths[a.ref_name] += a.fraction_ref_aligned()

    return mapping_quality, unitig_depths


def trim_dead_ends_based_on_miniasm_trimming(assembly_graph, miniasm_read_list):
    log.log_explanation('Contigs in the short-read assembly graph which end in dead ends may '
                        'contain bogus sequence near the dead end. Unicycler therefore uses the '
                        'read clipping values from the miniasm assembly to trim these dead ends '
                        'to only the parts which aligned well to long reads.')

    dead_end_trim_table = [['Segment', 'Previous length', 'Trimmed from start', 'Trimmed from end',
                            'Final length']]

    with open(miniasm_read_list, 'rt') as miniasm_read_names:
        for line in miniasm_read_names:
            if line.startswith('CONTIG_'):
                line = line.strip()
                contig_number = int(line.split('CONTIG_')[1].split(':')[0])
                start_dead_end = assembly_graph.starts_with_dead_end(contig_number)
                end_dead_end = assembly_graph.ends_with_dead_end(contig_number)
                if start_dead_end or end_dead_end:
                    contig = assembly_graph.segments[contig_number]
                    starting_length = contig.get_length()
                    contig_range = line.split(':')[1]
                    contig_start_pos, contig_end_pos = [int(x) for x in contig_range.split('-')]

                    contig_start_trim = contig_start_pos - 1  # 1-based range to Python 0-base range
                    contig_end_trim = contig.get_length() - contig_end_pos

                    if not start_dead_end:
                        contig_start_trim = 0
                    if not end_dead_end:
                        contig_end_trim = 0
                    if contig_start_trim > settings.MAX_MINIASM_DEAD_END_TRIM_SIZE:
                        contig_start_trim = 0
                    if contig_end_trim > settings.MAX_MINIASM_DEAD_END_TRIM_SIZE:
                        contig_end_trim = 0

                    if contig_start_trim and start_dead_end:
                        contig.trim_from_start(contig_start_trim)
                    if contig_end_trim and end_dead_end:
                        contig.trim_from_end(contig_end_trim)

                    ending_length = contig.get_length()
                    table_row = [int_to_str(contig_number),
                                 int_to_str(starting_length),
                                 int_to_str(contig_start_trim) if contig_start_trim else '-',
                                 int_to_str(contig_end_trim) if contig_end_trim else '-',
                                 int_to_str(ending_length)]
                    if contig_start_trim or contig_end_trim:
                        dead_end_trim_table.append(table_row)

    if len(dead_end_trim_table) > 1:
        print_table(dead_end_trim_table, fixed_col_widths=[7, 10, 7, 7, 10], alignments='LRRRR',
                    indent=0, left_align_header=False)
    else:
        log.log('No dead ends required trimming.')
