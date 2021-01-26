#!/usr/bin/env python3
"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module makes Unicycler's semi-global long read alignment available as a stand-alone program.
It is executed when a user runs `unicycler_align` (after installation) or
`unicycler_align-runner.py`.

Semi-global alignment does not penalise end gaps, but the alignment will continue until one of the
two sequences ends. This includes cases where the two sequences overlap and cases where one
sequence is contained within the other:

  TAGAA        GTGCCGGAACA         GGCCACAC     AGTAAGAT
  |||||          |||||||           |||||           |||||
ACTAGAACG        GCCGGAA       GGCTGGCCA           AAGATCTTG

This tool is intended for cases where the reads and reference are expected to match perfectly (or
at least as perfectly as error-prone long reads can match). An example of an appropriate case would
be if the reference sequences are assembled contigs of a bacterial strain and the long reads are
from the same strain.

Required inputs:
  1) FASTA file of one or more reference sequences
  2) FASTQ file of long reads

Output: SAM file of alignments

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import sys
import os
import argparse
import time
import random
import shutil
import math
from multiprocessing.dummy import Pool as ThreadPool
import threading
from .misc import int_to_str, float_to_str, check_file_exists, quit_with_error, \
    weighted_average_list, get_sequence_file_type, MyHelpFormatter, dim, magenta, colour,\
    get_default_thread_count
from .read_ref import load_references, load_long_reads
from .alignment import Alignment, AlignmentScoringScheme
from . import settings
from .minimap_alignment import load_minimap_alignments
from . import log

try:
    from .cpp_wrappers import semi_global_alignment, new_ref_seqs, add_ref_seq, \
        delete_ref_seqs, get_random_sequence_alignment_mean_and_std_dev, minimap_align_reads
except AttributeError as e:
    sys.exit('Error when importing C++ library: ' + str(e) + '\n'
             'Have you successfully built the library file using make?')


# Used to ensure that multiple threads writing to the same SAM file don't write at the same time.
SAM_WRITE_LOCK = threading.Lock()

# VERBOSITY controls how much the script prints to the screen.
# 0 = nothing is printed
# 1 = a relatively simple output is printed
# 2 = a more thorough output is printed, including details on each Seqan alignment
# 3 = even more output is printed, including stuff from the C++ code
# 4 = tons of stuff is printed, including all k-mer positions in each Seqan alignment
VERBOSITY = 0


def main():
    """
    If this script is run on its own, execution starts here.
    """
    # Fix the random seed so the program produces the same output every time it's run.
    random.seed(0)

    full_command = ' '.join(sys.argv)
    args = get_arguments()
    check_file_exists(args.ref)
    check_file_exists(args.reads)

    references = load_references(args.ref)
    read_dict, read_names, read_filename = load_long_reads(args.reads)
    scoring_scheme = AlignmentScoringScheme(args.scores)

    semi_global_align_long_reads(references, args.ref, read_dict, read_names, read_filename,
                                 args.threads, scoring_scheme, [args.low_score], args.keep_bad,
                                 args.min_len, args.sam, full_command, args.allowed_overlap,
                                 args.sensitivity, args.contamination, VERBOSITY)
    sys.exit(0)


def get_arguments():
    """
    Specifies the command line arguments required by the script.
    """
    terminal_width = shutil.get_terminal_size().columns
    os.environ['COLUMNS'] = str(terminal_width)

    parser = argparse.ArgumentParser(description='Unicycler align - a sensitive semi-global long '
                                                 'read aligner',
                                     formatter_class=MyHelpFormatter)

    parser.add_argument('--ref', type=str, required=True,
                        help='FASTA file containing one or more reference sequences')
    parser.add_argument('--reads', type=str, required=True,
                        help='FASTQ or FASTA file of long reads')
    parser.add_argument('--sam', type=str, required=True,
                        help='SAM file of resulting alignments')

    add_aligning_arguments(parser, True)

    parser.add_argument('--keep_bad', action='store_true',
                        help='Include alignments in the results even if they are below the low '
                             'score threshold (default: low-scoring alignments are discarded)')
    parser.add_argument('--sensitivity', type=int, default=0,
                        help='A number from 0 (least sensitive) to 3 (most sensitive)')
    parser.add_argument('--threads', type=int, required=False, default=get_default_thread_count(),
                        help='Number of threads used (default: number of CPUs, up to ' +
                             str(settings.MAX_AUTO_THREAD_COUNT) + ')')
    parser.add_argument('--verbosity', type=int, required=False, default=1,
                        help='Level of stdout information (0 to 4)')
    parser.add_argument('--min_len', type=int, required=False, default=100,
                        help='Minimum alignment length (bp) - exclude alignments shorter than this '
                             'length')
    parser.add_argument('--allowed_overlap', type=int, required=False, default=100,
                        help='Allow this much overlap between alignments in a single read')

    args = parser.parse_args()

    global VERBOSITY
    VERBOSITY = args.verbosity
    log.logger = log.Log(log_filename=None, stdout_verbosity_level=VERBOSITY)

    fix_up_arguments(args)

    return args


def add_aligning_arguments(parser, show_help):
    """
    Adds some aligning-specific arguments to the parser. These are in a separate function because
    other tools (e.g. unicycler.py, unicycler_check.py) can use this function too.
    """
    parser.add_argument('--contamination', required=False,
                        help='FASTA file of known contamination in long reads'
                             if show_help else argparse.SUPPRESS)
    parser.add_argument('--scores', type=str, required=False, default='3,-6,-5,-2',
                        help='Comma-delimited string of alignment scores: match, mismatch, '
                             'gap open, gap extend'
                             if show_help else argparse.SUPPRESS)
    parser.add_argument('--low_score', type=float, required=False,
                        help='Score threshold - alignments below this are considered poor '
                             '(default: set threshold automatically)'
                             if show_help else argparse.SUPPRESS)


def fix_up_arguments(args):
    # If the user just said 'lambda' for the contamination, then we use the lambda phage FASTA
    # included with Unicycler.
    if args.contamination == 'lambda':
        lambda_phage_fasta = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                          'gene_data', 'lambda_phage.fasta')
        if not os.path.isfile(lambda_phage_fasta):
            quit_with_error('Could not find lambda_phage.fasta - please reinstall')
        args.contamination = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                          'gene_data', 'lambda_phage.fasta')

    # If a contamination file was given, make sure it's FASTA.
    if args.contamination:
        args.contamination = os.path.abspath(args.contamination)
        contamination_type = get_sequence_file_type(args.contamination)
        if contamination_type != 'FASTA':
            quit_with_error('Long read contamination file must be FASTA format')


def semi_global_align_long_reads(references, ref_fasta, read_dict, read_names, reads_fastq,
                                 threads, scoring_scheme, low_score_threshold_list, keep_bad,
                                 min_align_length, sam_filename, full_command, allowed_overlap,
                                 sensitivity_level, contamination_fasta, verbosity=None,
                                 stdout_header='Aligning reads', display_low_score=True,
                                 single_copy_segment_names=None):
    """
    This function does the primary work of this module: aligning long reads to references in an
    end-gap-free, semi-global manner. It returns a dictionary of Read objects which contain their
    alignments.
    The low score threshold is taken as a list so the function can alter it and the caller can
    get the altered value.
    """
    if sensitivity_level is None:
        sensitivity_level = 0
    if verbosity:
        global VERBOSITY
        VERBOSITY = verbosity

    if single_copy_segment_names is None:
        single_copy_segment_names = set()

    # If the user supplied a low score threshold, we use that. Otherwise, we'll use the median
    # score minus three times the MAD.
    if display_low_score and verbosity > 0:
        log.log_section_header('Determining low score threshold')
        log.log_explanation('Before conducting semi-global alignment of the long reads to the '
                            'assembly graph, Unicycler must determine a minimum alignment '
                            'score threshold such that nonsense alignments are excluded. '
                            'To choose a threshold automatically, it examines alignments between '
                            'random sequences and selects a score a few standard deviations above '
                            'the mean.')
    low_score_threshold = low_score_threshold_list[0]
    if low_score_threshold is not None:
        if display_low_score and verbosity > 0:
            log.log('Using user-supplied threshold: ' + float_to_str(low_score_threshold, 2))
    else:
        if display_low_score and verbosity > 0:
            log.log('Automatically choosing a threshold using random alignment scores.\n')
        std_devs_over_mean = settings.AUTO_SCORE_STDEV_ABOVE_RANDOM_ALIGNMENT_MEAN
        low_score_threshold, rand_mean, rand_std_dev = get_auto_score_threshold(scoring_scheme,
                                                                                std_devs_over_mean)
        low_score_threshold_list[0] = low_score_threshold
        if display_low_score and verbosity > 0:
            log.log('Random alignment mean score: ' + float_to_str(rand_mean, 2))
            log.log('         standard deviation: ' + float_to_str(rand_std_dev, 2, rand_mean))
            log.log('        Low score threshold: ' + float_to_str(rand_mean, 2) + ' + (' +
                    str(std_devs_over_mean) + ' x ' + float_to_str(rand_std_dev, 2) + ') = ' +
                    float_to_str(low_score_threshold, 2))

    using_contamination = contamination_fasta is not None
    if using_contamination:
        references += load_references(contamination_fasta, contamination=True)

    reference_dict = {x.name: x for x in references}

    if verbosity > 0:
        log.log_section_header('Aligning reads with minimap', verbosity=2)
    minimap_alignments_str = minimap_align_reads(ref_fasta, reads_fastq, threads, 0, 'default')
    minimap_alignments = load_minimap_alignments(minimap_alignments_str)
    if verbosity > 0:
        log.log('', 3)
        log.log('Done! ' + str(len(minimap_alignments)) + ' out of ' +
                str(len(read_dict)) + ' reads aligned', 2)

    # Create the SAM file.
    if sam_filename:
        with open(sam_filename, 'w') as sam_file:
            # Header line.
            sam_file.write('@HD' + '\t')
            sam_file.write('VN:1.5' + '\t')
            sam_file.write('SO:unknown' + '\n')

            # Reference lines.
            for ref in references:
                sam_file.write('@SQ' + '\t')
                sam_file.write('SN:' + ref.name + '\t')
                sam_file.write('LN:' + str(ref.get_length()) + '\n')

            # Program line.
            sam_file.write('@PG' + '\t')
            sam_file.write('ID:' + 'unicycler_align')
            if full_command:
                sam_file.write('\tCL:' + full_command + '\t')
            sam_file.write('SC:' + str(scoring_scheme) + '\n')

    reads_to_align = [read_dict[x] for x in read_names]

    num_alignments = len(reads_to_align)
    if verbosity > 0:
        log.log_section_header(stdout_header)
    if VERBOSITY == 1:
        log.log_progress_line(0, num_alignments)
    completed_count = 0

    # Create a C++ ReferenceSeqs object and add each reference sequence.
    ref_seqs_ptr = new_ref_seqs()
    for ref in references:
        add_ref_seq(ref_seqs_ptr, ref.name, ref.sequence)

    # If single-threaded, just do the work in a simple loop.
    if threads == 1:
        for read in reads_to_align:
            output = seqan_alignment(read, reference_dict, scoring_scheme, ref_seqs_ptr,
                                     low_score_threshold, keep_bad, min_align_length,
                                     sam_filename, allowed_overlap, minimap_alignments[read.name],
                                     sensitivity_level, single_copy_segment_names)
            completed_count += 1
            if VERBOSITY == 1:
                log.log_progress_line(completed_count, num_alignments)
            if VERBOSITY > 1:
                fraction = str(completed_count) + '/' + str(num_alignments) + ': '
                log.log(fraction + output + '\n', 2, end='')

    # If multi-threaded, use a thread pool.
    else:
        pool = ThreadPool(threads)
        arg_list = []
        for read in reads_to_align:
            arg_list.append((read, reference_dict, scoring_scheme, ref_seqs_ptr,
                             low_score_threshold, keep_bad, min_align_length,
                             sam_filename, allowed_overlap, minimap_alignments[read.name],
                             sensitivity_level, single_copy_segment_names))

        # If the verbosity is 1, then the order doesn't matter, so use imap_unordered to deliver
        # the results evenly. If the verbosity is higher, deliver the results in order with imap.
        if VERBOSITY > 1:
            imap_function = pool.imap
        else:
            imap_function = pool.imap_unordered

        for output in imap_function(seqan_alignment_one_arg, arg_list):
            completed_count += 1
            if VERBOSITY == 1:
                log.log_progress_line(completed_count, num_alignments)
            if VERBOSITY > 1:
                fraction = str(completed_count) + '/' + str(num_alignments) + ': '
                log.log(fraction + output + '\n', 2, end='')

    # We're done with the C++ ReferenceSeqs object, so delete it now.
    delete_ref_seqs(ref_seqs_ptr)

    if VERBOSITY == 1:
        log.log_progress_line(completed_count, completed_count, end_newline=True)

    if verbosity > 0:
        print_alignment_summary_table(read_dict, VERBOSITY, using_contamination)
    return read_dict


def get_percent_contamination(read_dict):
    """
    Returns the number and percentage of reads which mostly align to contamination, both by base
    count and read count.
    """
    contamination_count, some_alignment_count = 0, 0
    contamination_bases, some_alignment_bases = 0, 0
    for read in read_dict.values():
        if read.get_fraction_aligned() > 0.0:
            some_alignment_count += 1
            some_alignment_bases += read.get_length()
            if read.mostly_aligns_to_contamination():
                contamination_count += 1
                contamination_bases += read.get_length()

    if some_alignment_count == 0:
        percentage_by_count = 0.0
    else:
        percentage_by_count = 100.0 * contamination_count / some_alignment_count

    if some_alignment_bases == 0:
        percentage_by_bases = 0.0
    else:
        percentage_by_bases = 100.0 * contamination_bases / some_alignment_bases

    return contamination_count, percentage_by_count, contamination_bases, percentage_by_bases


def print_alignment_summary_table(read_dict, verbosity, using_contamination):
    """
    Outputs a summary of the reads' alignments, grouping them by fully aligned, partially aligned
    and unaligned.
    """
    fully_aligned, partially_aligned, unaligned = group_reads_by_fraction_aligned(read_dict)
    ref_bases_aligned = 0
    for read in read_dict.values():
        ref_bases_aligned += read.get_reference_bases_aligned()
    log.log_section_header('Read alignment summary', single_newline=(verbosity > 1))
    max_v = max(len(read_dict), ref_bases_aligned)

    if using_contamination:
        contaminant_reads, contaminant_read_per, contaminant_bases, contaminant_base_per = \
            get_percent_contamination(read_dict)
    else:
        contaminant_reads, contaminant_read_per, contaminant_bases, contaminant_base_per = \
            None, None, None, None
    log.log('Total read count:        ' + int_to_str(len(read_dict), max_v))
    log.log('Fully aligned reads:     ' + int_to_str(len(fully_aligned), max_v))
    log.log('Partially aligned reads: ' + int_to_str(len(partially_aligned), max_v))
    if partially_aligned:
        log.log(dim(', '.join([x.name for x in partially_aligned])), 3)
        log.log('', 3)
    log.log('Unaligned reads:         ' + int_to_str(len(unaligned), max_v))
    if unaligned:
        log.log(dim(', '.join([x.name for x in unaligned])), 3)
        log.log('', 3)

    if using_contamination:
        log.log('Contaminant reads:       ' + int_to_str(contaminant_reads, max_v))
        log.log('Contaminant reads:       ' + float_to_str(contaminant_read_per, 1, max_v) + '%')

    log.log('Total bases aligned:     ' + int_to_str(ref_bases_aligned, max_v) + ' bp')
    if using_contamination:
        log.log('Contaminant bases:       ' + int_to_str(contaminant_bases, max_v) + ' bp')
        log.log('Contaminant bases:       ' + float_to_str(contaminant_base_per, 1, max_v) + '%')

    identities = []
    lengths = []
    for read in fully_aligned + partially_aligned:
        identities += [x.percent_identity for x in read.alignments]
        lengths += [x.get_aligned_ref_length() for x in read.alignments]
    mean_identity = weighted_average_list(identities, lengths)
    log.log('Mean alignment identity: ' + float_to_str(mean_identity, 1, max_v) + '%')


def load_sam_alignments(sam_filename, read_dict, reference_dict, scoring_scheme):
    """
    This function returns a list of Alignment objects from the given SAM file.
    """
    log.log_section_header('Loading alignments')

    sam_lines = []
    sam_file = open(sam_filename, 'rt')
    for line in sam_file:
        line = line.strip()
        if line and not line.startswith('@') and line.split('\t', 3)[2] != '*':
            sam_lines.append(line)
    num_alignments = sum(1 for line in open(sam_filename) if not line.startswith('@'))
    if not num_alignments:
        return []
    log.log_progress_line(0, num_alignments)

    if not sam_lines:
        log.log('No alignments to load')
        return []

    sam_alignments = []
    last_progress = 0.0
    step = settings.LOADING_ALIGNMENTS_PROGRESS_STEP
    for line in sam_lines:
        sam_alignments.append(Alignment(sam_line=line, read_dict=read_dict,
                                        reference_dict=reference_dict,
                                        scoring_scheme=scoring_scheme))
        progress = 100.0 * len(sam_alignments) / num_alignments
        progress_rounded_down = math.floor(progress / step) * step
        if progress == 100.0 or progress_rounded_down > last_progress:
            log.log_progress_line(len(sam_alignments), num_alignments)
            last_progress = progress_rounded_down

    # At this point, we should have loaded num_alignments alignments. But check to make sure and
    # fix up the progress line if any didn't load.
    if len(sam_alignments) < num_alignments:
        log.log_progress_line(len(sam_alignments), len(sam_alignments))
    log.log('')

    return sam_alignments


def seqan_alignment_one_arg(all_args):
    """
    This is just a one-argument version of seqan_alignment to make it easier to use that function
    in a thread pool.
    """
    read, reference_dict, scoring_scheme, ref_seqs_ptr, low_score_threshold, keep_bad, \
        min_align_length, sam_filename, allowed_overlap, minimap_alignments, \
        sensitivity_level, single_copy_segment_names = all_args
    return seqan_alignment(read, reference_dict, scoring_scheme, ref_seqs_ptr,
                           low_score_threshold, keep_bad, min_align_length,
                           sam_filename, allowed_overlap, minimap_alignments, sensitivity_level,
                           single_copy_segment_names)


def seqan_alignment(read, reference_dict, scoring_scheme, ref_seqs_ptr, low_score_threshold,
                    keep_bad, min_align_length, sam_filename, allowed_overlap,
                    minimap_alignments, sensitivity_level, single_copy_segment_names):
    """
    Aligns a single read against all reference sequences using Seqan.
    """
    start_time = time.time()
    output = ''

    # Don't bother trying to align reads too short to have a good alignment.
    if read.get_length() < min_align_length:
        if VERBOSITY > 1:
            output += '  too short to align\n'
    else:
        minimap_alignments_str = ';'.join([x.get_concise_string() for x in minimap_alignments])
        alignment_strings = []

        # Try at each sensitivity up to the current level. E.g. if sensitivity level is 2,
        # we try the Seqan alignment at levels 0, 1 and 2. This produces many redundant
        # alignments but we'll filter them out later.
        for sensitivity in range(0, sensitivity_level+1):
            results = semi_global_alignment(read.name, read.sequence, VERBOSITY,
                                            minimap_alignments_str, ref_seqs_ptr,
                                            scoring_scheme.match, scoring_scheme.mismatch,
                                            scoring_scheme.gap_open, scoring_scheme.gap_extend,
                                            low_score_threshold, keep_bad, sensitivity).split(';')

            alignment_strings += results[:-1]
            output += results[-1]

            for alignment_string in alignment_strings:
                alignment = Alignment(seqan_output=alignment_string, read=read,
                                      reference_dict=reference_dict, scoring_scheme=scoring_scheme)
                read.alignments.append(alignment)

        if VERBOSITY > 2:
            if not alignment_strings:
                output += '  None\n'
            else:
                output += 'All Seqan alignments (time to align = ' + \
                          float_to_str(time.time() - start_time, 3) + ' s):\n'
                output += read.get_alignment_table()

        read.remove_conflicting_alignments(allowed_overlap)
        if not keep_bad:
            read.remove_low_score_alignments(low_score_threshold)
        read.remove_short_alignments(min_align_length)

        if VERBOSITY > 2:
            output += 'Final alignments:\n'
        if VERBOSITY > 1:
            if read.alignments:
                output += read.get_alignment_table()
            else:
                output += '  None\n'

        # Write alignments to SAM.
        if sam_filename and read.alignments:
            SAM_WRITE_LOCK.acquire()
            sam_file = open(sam_filename, 'a')
            for alignment in read.alignments:
                if not alignment.ref.name.startswith('CONTAMINATION_'):
                    sam_file.write(alignment.get_sam_line())
            sam_file.close()
            SAM_WRITE_LOCK.release()

    # Colour the output title based on the alignment quality.
    if read.mostly_aligns_to_contamination() or not read.alignments:
        title_colour = 'red'
    elif read.aligns_to_multiple_single_copy_segments(single_copy_segment_names) and \
            read.get_fraction_aligned() > 0.8:
        title_colour = 'green'
    else:
        title_colour = 'normal'
    output_title = colour(str(read), title_colour) + '\n'

    formatted_output = ''.join(magenta(x[7:]) if x.startswith('R_code:') else dim(x)
                               for x in output.splitlines(True))

    return output_title + formatted_output


def group_reads_by_fraction_aligned(read_dict):
    """
    Groups reads into three lists:
      1) Fully aligned
      2) Partially aligned
      3) Unaligned
    """
    fully_aligned_reads = []
    partially_aligned_reads = []
    unaligned_reads = []
    for read in read_dict.values():
        fraction_aligned = read.get_fraction_aligned()
        if fraction_aligned == 1.0:
            fully_aligned_reads.append(read)
        elif fraction_aligned == 0.0:
            unaligned_reads.append(read)
        else:
            partially_aligned_reads.append(read)
    return fully_aligned_reads, partially_aligned_reads, unaligned_reads


def get_auto_score_threshold(scoring_scheme, std_devs_over_mean):
    """
    This function determines a good low score threshold for the alignments. To do this it examines
    the distribution of scores acquired by aligning random sequences.
    """
    # If the scoring scheme is a typical one, don't actually do the random alignments now - just
    # use precomputed values (made with a lot of iterations so they should be pretty good).
    scoring_scheme_str = str(scoring_scheme)
    if scoring_scheme_str == '1,0,0,0':
        mean, std_dev = 50.225667, 2.467919
    elif scoring_scheme_str == '0,-1,-1,-1':
        mean, std_dev = 49.024927, 2.724548
    elif scoring_scheme_str == '1,-1,-1,-1':
        mean, std_dev = 51.741783, 2.183467
    elif scoring_scheme_str == '5,-4,-8,-6':   # GraphMap
        mean, std_dev = 42.707636, 2.435548
    elif scoring_scheme_str == '5,-6,-10,0':   # BLASR
        mean, std_dev = 58.65047, 0.853201
    elif scoring_scheme_str == '2,-5,-2,-1':   # BWA-MEM
        mean, std_dev = 72.712148, 0.95266
    elif scoring_scheme_str == '1,-3,-5,-2':   # CUSHAW2 / blastn-short
        mean, std_dev = 46.257408, 2.162765
    elif scoring_scheme_str == '5,-11,-2,-4':  # proovread
        mean, std_dev = 73.221967, 1.363692
    elif scoring_scheme_str == '3,-6,-5,-2':   # Unicycler-align
        mean, std_dev = 61.656918, 1.314624
    elif scoring_scheme_str == '2,-3,-5,-2':   # blastn / dc-megablast
        mean, std_dev = 47.453862, 1.985947
    elif scoring_scheme_str == '1,-2,0,0':     # megablast
        mean, std_dev = 81.720641, 0.77204
    elif scoring_scheme_str == '0,-6,-5,-3':   # Bowtie2 end-to-end
        mean, std_dev = 62.647055, 1.738603
    elif scoring_scheme_str == '2,-6,-5,-3':   # Bowtie2 local
        mean, std_dev = 59.713806, 1.641191
    elif scoring_scheme_str == '1,-4,-6,-1':   # BWA
        mean, std_dev = 60.328393, 1.176776

    # If scheme doesn't match any of the above, then we have to actually do the random alignments.
    else:
        mean, std_dev = get_random_sequence_alignment_mean_and_std_dev(100, 25000, scoring_scheme)

    threshold = mean + (std_devs_over_mean * std_dev)

    # Keep the threshold bounded to sane levels.
    threshold = min(threshold, 95.0)
    threshold = max(threshold, 50.0)

    return threshold, mean, std_dev
