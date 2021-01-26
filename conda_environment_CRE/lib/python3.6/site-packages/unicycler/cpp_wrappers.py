"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

Unicycler makes use of several C++ functions which are in cpp_functions.so. This module uses ctypes
to wrap them in similarly named Python functions.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import os
from ctypes import CDLL, cast, c_char_p, c_int, c_uint, c_ulong, c_double, c_void_p, c_bool, \
    c_float, POINTER
from .misc import quit_with_error


SO_FILE = 'cpp_functions.so'
SO_FILE_FULL = os.path.join(os.path.dirname(os.path.realpath(__file__)), SO_FILE)
if not os.path.isfile(SO_FILE_FULL):
    quit_with_error('could not find ' + SO_FILE + '\n' +
                    "Please reinstall Unicycler or run make from Unicycler's source directory")
C_LIB = CDLL(SO_FILE_FULL)



# This is the big semi-global C++ Seqan alignment function at the heart of the aligner.
C_LIB.semiGlobalAlignment.argtypes = [c_char_p,  # Read name
                                      c_char_p,  # Read sequence
                                      c_int,     # Verbosity
                                      c_char_p,  # Minimap alignment info
                                      c_void_p,  # KmerPositions pointer
                                      c_int,     # Match score
                                      c_int,     # Mismatch score
                                      c_int,     # Gap open score
                                      c_int,     # Gap extension score
                                      c_double,  # Low score threshold
                                      c_bool,    # Return bad alignments
                                      c_int]     # Sensitivity level
C_LIB.semiGlobalAlignment.restype = c_void_p     # String describing alignments

def semi_global_alignment(read_name, read_sequence, verbosity, minimap_alignments_str,
                          kmer_positions_ptr, match_score, mismatch_score, gap_open_score,
                          gap_extend_score, low_score_threshold, keep_bad, sensitivity_level):
    ptr = C_LIB.semiGlobalAlignment(read_name.encode('utf-8'), read_sequence.encode('utf-8'),
                                    verbosity, minimap_alignments_str.encode('utf-8'),
                                    kmer_positions_ptr, match_score, mismatch_score,
                                    gap_open_score, gap_extend_score, low_score_threshold,
                                    keep_bad, sensitivity_level)
    return c_string_to_python_string(ptr)



# This function does an exhaustive semi-global alignment (nothing fancy, only suitable for short
# sequences).
C_LIB.semiGlobalAlignmentExhaustive.argtypes = [c_char_p,  # Sequence 1
                                                c_char_p,  # Sequence 2
                                                c_int,  # Match score
                                                c_int,  # Mismatch score
                                                c_int,  # Gap open score
                                                c_int]  # Gap extension score
C_LIB.semiGlobalAlignmentExhaustive.restype = c_void_p  # String describing alignment

def semi_global_alignment_exhaustive(sequence_1, sequence_2, scoring_scheme):
    ptr = C_LIB.semiGlobalAlignmentExhaustive(sequence_1.encode('utf-8'),
                                              sequence_2.encode('utf-8'),
                                              scoring_scheme.match, scoring_scheme.mismatch,
                                              scoring_scheme.gap_open, scoring_scheme.gap_extend)
    return c_string_to_python_string(ptr)



# This is the global alignment function mainly used to compare read consensus sequences to assembly
# graph paths.
C_LIB.fullyGlobalAlignment.argtypes = [c_char_p,  # Sequence 1
                                       c_char_p,  # Sequence 2
                                       c_int,  # Match score
                                       c_int,  # Mismatch score
                                       c_int,  # Gap open score
                                       c_int,  # Gap extension score
                                       c_bool,  # Use banding
                                       c_int]  # Band size
C_LIB.fullyGlobalAlignment.restype = c_void_p  # String describing alignment

def fully_global_alignment(sequence_1, sequence_2, scoring_scheme, use_banding, band_size):
    ptr = C_LIB.fullyGlobalAlignment(sequence_1.encode('utf-8'), sequence_2.encode('utf-8'),
                                     scoring_scheme.match, scoring_scheme.mismatch,
                                     scoring_scheme.gap_open, scoring_scheme.gap_extend,
                                     use_banding, band_size)
    return c_string_to_python_string(ptr)



# This is the mostly-global alignment function mainly used to compare potential path sequences to
# a read consensus. It is 'mostly-global' because there are free end gaps in the first sequence,
# so the path isn't penalised for not being complete.
C_LIB.pathAlignment.argtypes = [c_char_p,  # Sequence 1
                                c_char_p,  # Sequence 2
                                c_int,  # Match score
                                c_int,  # Mismatch score
                                c_int,  # Gap open score
                                c_int,  # Gap extension score
                                c_bool,  # Use banding
                                c_int]  # Band size
C_LIB.pathAlignment.restype = c_void_p  # String describing alignment

def path_alignment(partial_seq, full_seq, scoring_scheme, use_banding, band_size):
    ptr = C_LIB.pathAlignment(partial_seq.encode('utf-8'), full_seq.encode('utf-8'),
                              scoring_scheme.match, scoring_scheme.mismatch,
                              scoring_scheme.gap_open, scoring_scheme.gap_extend,
                              use_banding, band_size)
    return c_string_to_python_string(ptr)



# This function cleans up the heap memory for the C strings returned by the other C functions. It
# must be called after them.
C_LIB.freeCString.argtypes = [c_void_p]
C_LIB.freeCString.restype = None

def c_string_to_python_string(c_string):
    """
    This function casts a C string to a Python string and then calls a function to delete the C
    string from the heap.
    """
    python_string = cast(c_string, c_char_p).value.decode()
    C_LIB.freeCString(c_string)
    return python_string



# These functions make/delete a C++ object that will hold reference sequences for quick access.
C_LIB.newRefSeqs.argtypes = []
C_LIB.newRefSeqs.restype = c_void_p

def new_ref_seqs():
    return C_LIB.newRefSeqs()

C_LIB.addRefSeq.argtypes = [c_void_p,  # SeqMap pointer
                            c_char_p,  # Name
                            c_char_p]  # Sequence
C_LIB.addRefSeq.restype = None

def add_ref_seq(ref_seqs_ptr, name, sequence):
    C_LIB.addRefSeq(ref_seqs_ptr, name.encode('utf-8'), sequence.encode('utf-8'))

C_LIB.deleteRefSeqs.argtypes = [c_void_p]
C_LIB.deleteRefSeqs.restype = None

def delete_ref_seqs(ref_seqs_ptr):
    C_LIB.deleteRefSeqs(ref_seqs_ptr)



# This function gets the mean and standard deviation of alignments between random sequences.
C_LIB.getRandomSequenceAlignmentScores.argtypes = [c_int,  # Random sequence length
                                                   c_int,  # Count
                                                   c_int,  # Match score
                                                   c_int,  # Mismatch score
                                                   c_int,  # Gap open score
                                                   c_int]  # Gap extension score
C_LIB.getRandomSequenceAlignmentScores.restype = c_void_p

def get_random_sequence_alignment_mean_and_std_dev(seq_length, count, scoring_scheme):
    ptr = C_LIB.getRandomSequenceAlignmentScores(seq_length, count,
                                                 scoring_scheme.match, scoring_scheme.mismatch,
                                                 scoring_scheme.gap_open, scoring_scheme.gap_extend)
    return_str = c_string_to_python_string(ptr)
    return_parts = return_str.split(',')
    return float(return_parts[0]), float(return_parts[1])



# This function gets the mean and standard deviation of alignments between random sequences.
C_LIB.getRandomSequenceAlignmentErrorRates.argtypes = [c_int,  # Random sequence length
                                                       c_int,  # Count
                                                       c_int,  # Match score
                                                       c_int,  # Mismatch score
                                                       c_int,  # Gap open score
                                                       c_int]  # Gap extension score
C_LIB.getRandomSequenceAlignmentErrorRates.restype = c_void_p

def get_random_sequence_alignment_error_rates(seq_length, count, scoring_scheme):

    ptr = C_LIB.getRandomSequenceAlignmentErrorRates(seq_length, count, scoring_scheme.match,
                                                     scoring_scheme.mismatch,
                                                     scoring_scheme.gap_open,
                                                     scoring_scheme.gap_extend)
    return c_string_to_python_string(ptr)



# This function gets the mean and standard deviation of alignments between random sequences.
C_LIB.simulateDepths.argtypes = [POINTER(c_int),  # Alignment lengths
                                 c_int,  # Alignment count
                                 c_int,  # Reference length
                                 c_int,  # Iterations
                                 c_int]  # Threads
C_LIB.simulateDepths.restype = c_void_p

def simulate_depths(read_lengths, ref_length, iterations, threads):
    # noinspection PyCallingNonCallable
    read_lengths_array = (c_int * len(read_lengths))(*read_lengths)
    ptr = C_LIB.simulateDepths(read_lengths_array, len(read_lengths), ref_length, iterations,
                               threads)
    return c_string_to_python_string(ptr)



# This function gets the mean and standard deviation of alignments between random sequences.
C_LIB.multipleSequenceAlignment.argtypes = [POINTER(c_char_p),  # Sequences
                                            POINTER(c_char_p),  # Qualities
                                            c_ulong,  # Count
                                            c_uint,  # Bandwidth
                                            c_int,  # Match score
                                            c_int,  # Mismatch score
                                            c_int,  # Gap open score
                                            c_int]  # Gap extension score
C_LIB.multipleSequenceAlignment.restype = c_void_p

def consensus_alignment(sequences, qualities, scoring_scheme, bandwidth=1000):
    count = len(sequences)
    if not count:  # At least one sequence is required.
        return "", []

    if len(qualities) < len(sequences):
        qualities += [""] * (len(sequences) - len(qualities))

    sequences = [x.encode('utf-8') for x in sequences]
    qualities = [x.encode('utf-8') for x in qualities]

    # noinspection PyCallingNonCallable
    sequences = (c_char_p * len(sequences))(*sequences)
    # noinspection PyCallingNonCallable
    qualities = (c_char_p * len(qualities))(*qualities)

    ptr = C_LIB.multipleSequenceAlignment(sequences, qualities, count, bandwidth,
                                          scoring_scheme.match, scoring_scheme.mismatch,
                                          scoring_scheme.gap_open, scoring_scheme.gap_extend)
    result = c_string_to_python_string(ptr)
    result_parts = result.split(';')
    consensus = result_parts[0]
    scores = [float(x) for x in result_parts[1].split(',')]
    return consensus, scores



# These functions conduct a minimap alignment between reads and reference.
C_LIB.minimapAlignReads.argtypes = [c_char_p,  # Reference FASTA filename
                                    c_char_p,  # Reads FASTQ filename
                                    c_int,     # Threads
                                    c_int,     # Sensitivity level
                                    c_int]     # Settings preset
C_LIB.minimapAlignReads.restype = c_void_p     # String describing alignments

def minimap_align_reads(reference_fasta, reads_fastq, threads, sensitivity_level,
                        preset_name='default'):
    preset = 0  # default
    if preset_name == 'read vs read':
        preset = 1
    elif preset_name == 'find contigs':
        preset = 2
    if preset_name == 'scrub reads with reads':
        preset = 1
    if preset_name == 'scrub assembly with reads':
        preset = 2
    ptr = C_LIB.minimapAlignReads(reference_fasta.encode('utf-8'), reads_fastq.encode('utf-8'),
                                  threads, sensitivity_level, preset)
    return c_string_to_python_string(ptr)

C_LIB.minimapAlignReadsWithSettings.argtypes = [c_char_p,  # Reference FASTA filename
                                                c_char_p,  # Reads FASTQ filename
                                                c_int,     # Threads
                                                c_bool,    # Whether an all vs all alignment (-S)
                                                c_int,     # K-mer size (-k)
                                                c_int,     # Minimiser size (-w)
                                                c_float,   # Merge fraction (-m)
                                                c_int,     # Minimum match length (-L)
                                                c_int,     # Maximum minimiser gap (-g)
                                                c_int,     # Bandwidth radius (-r)
                                                c_int]     # Minimum minimiser count (-c)
C_LIB.minimapAlignReadsWithSettings.restype = c_void_p     # String describing alignments


def minimap_align_reads_with_settings(reference_fasta, reads_fastq, threads, all_vs_all=False,
                                      kmer_size=15, minimiser_size=10, merge_fraction=0.5,
                                      min_match_len=40, max_gap=10000, bandwidth=500, min_count=4):
    ptr = C_LIB.minimapAlignReadsWithSettings(reference_fasta.encode('utf-8'),
                                              reads_fastq.encode('utf-8'), threads, all_vs_all,
                                              kmer_size, minimiser_size, merge_fraction,
                                              min_match_len, max_gap, bandwidth, min_count)
    return c_string_to_python_string(ptr)



# This function conducts a miniasm assembly
C_LIB.miniasmAssembly.argtypes = [c_char_p,  # Reads FASTQ filename
                                  c_char_p,  # Overlaps PAF filename
                                  c_char_p,  # Output GFA filename
                                  c_int]     # Min depth
C_LIB.miniasmAssembly.restype = None         # No return value (function creates GFA files)

def miniasm_assembly(reads_fastq, overlaps_paf, output_gfa, min_depth):
    C_LIB.miniasmAssembly(reads_fastq.encode('utf-8'), overlaps_paf.encode('utf-8'),
                          output_gfa.encode('utf-8'), min_depth)



# This is the overlap alignment function to see how much sequence 1 and sequence 2 overlap.
C_LIB.overlapAlignment.argtypes = [c_char_p,  # Sequence 1
                                   c_char_p,  # Sequence 2
                                   c_int,     # Match score
                                   c_int,     # Mismatch score
                                   c_int,     # Gap open score
                                   c_int,     # Gap extension score
                                   c_int]     # Guess overlap
C_LIB.overlapAlignment.restype = c_void_p     # String describing alignment

def overlap_alignment(sequence_1, sequence_2, scoring_scheme, guess_overlap):
    ptr = C_LIB.overlapAlignment(sequence_1.encode('utf-8'), sequence_2.encode('utf-8'),
                                 scoring_scheme.match, scoring_scheme.mismatch,
                                 scoring_scheme.gap_open, scoring_scheme.gap_extend, guess_overlap)
    result = c_string_to_python_string(ptr)
    return (int(x) for x in result.split(','))


# When s1 is expected to be at the start of s2, this function will align them to give the s2
# position where s1 ends.
C_LIB.startAlignment.argtypes = [c_char_p,  # Sequence 1
                                 c_char_p,  # Sequence 2
                                 c_int,     # Match score
                                 c_int,     # Mismatch score
                                 c_int,     # Gap open score
                                 c_int]     # Gap extension score
C_LIB.startAlignment.restype = c_int        # Seq 2 position at seq 1 end

def start_seq_alignment(sequence_1, sequence_2, scoring_scheme):
    return C_LIB.startAlignment(sequence_1.encode('utf-8'), sequence_2.encode('utf-8'),
                                scoring_scheme.match, scoring_scheme.mismatch,
                                scoring_scheme.gap_open, scoring_scheme.gap_extend)


# When s1 is expected to be at the end of s2, this function will align them to give the s2 position
# where s1 starts.
C_LIB.endAlignment.argtypes = [c_char_p,  # Sequence 1
                               c_char_p,  # Sequence 2
                               c_int,     # Match score
                               c_int,     # Mismatch score
                               c_int,     # Gap open score
                               c_int]     # Gap extension score
C_LIB.endAlignment.restype = c_int        # Seq 2 position at seq 1 start

def end_seq_alignment(sequence_1, sequence_2, scoring_scheme):
    return C_LIB.endAlignment(sequence_1.encode('utf-8'), sequence_2.encode('utf-8'),
                              scoring_scheme.match, scoring_scheme.mismatch,
                              scoring_scheme.gap_open, scoring_scheme.gap_extend)


# This function does the Unicycler-scrub splitting in C++ (for speed).
C_LIB.splitSequences.argtypes = [c_char_p,  # Alignments string
                                 c_int,     # Sequence length
                                 c_double,  # Starting score
                                 c_int,     # Positive score feather size
                                 c_int,     # Negative score feather size
                                 c_double,  # Positive score scaling factor
                                 c_int]     # Split adjustment
C_LIB.splitSequences.restype = c_void_p     # String describing pos/neg ranges of sequence

def split_sequences_cpp(alignments, seq_len, parameters):
    alignments_string = ';'.join(a.get_string_for_cpp_scrub() for a in alignments)
    ptr = C_LIB.splitSequences(alignments_string.encode('utf-8'), seq_len,
                               parameters.starting_score, parameters.pos_score_feather_size,
                               parameters.neg_score_feather_size,
                               parameters.pos_score_scaling_factor, parameters.split_adjustment)
    ranges_str = c_string_to_python_string(ptr)
    pos_ranges_str, neg_ranges_str = ranges_str.split(';')
    pos_ranges, neg_ranges = [], []
    if pos_ranges_str:
        for range in pos_ranges_str.split(','):
            range_parts = range.split('-')
            pos_ranges.append((int(range_parts[0]), int(range_parts[1])))
    if neg_ranges_str:
        for range in neg_ranges_str.split(','):
            range_parts = range.split('-')
            neg_ranges.append((int(range_parts[0]), int(range_parts[1])))
    return pos_ranges, neg_ranges
