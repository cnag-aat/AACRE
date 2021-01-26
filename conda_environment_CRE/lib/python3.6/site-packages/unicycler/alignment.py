"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module contains classes for alignments and related functions.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import re
from .misc import get_nice_header, reverse_complement, float_to_str


class AlignmentScoringScheme(object):
    """
    This class holds an alignment scoring scheme.
    """

    def __init__(self, scheme_string):
        scheme_parts = scheme_string.split(',')

        # Default scoring scheme
        self.match = 3
        self.mismatch = -6
        self.gap_open = -5
        self.gap_extend = -2

        if len(scheme_parts) == 4:
            self.match = int(scheme_parts[0])
            self.mismatch = int(scheme_parts[1])
            self.gap_open = int(scheme_parts[2])
            self.gap_extend = int(scheme_parts[3])

    def __repr__(self):
        return str(self.match) + ',' + str(self.mismatch) + ',' + str(self.gap_open) + ',' + \
               str(self.gap_extend)

    def get_full_string(self):
        """
        Returns a string verbosely describing the scheme.
        """
        return 'match = ' + str(self.match) + ', mismatch = ' + str(self.mismatch) + \
               ', gap open = ' + str(self.gap_open) + ', gap extend = ' + str(self.gap_extend)


class Alignment(object):
    """
    This class describes an alignment between a long read and a contig.
    It can be constructed either from a SAM line or from the C++ Seqan output.
    """

    def __init__(self,
                 sam_line=None, read_dict=None,
                 seqan_output=None, read=None,
                 reference_dict=None, scoring_scheme=None):

        # Make sure we have the appropriate inputs for one of the two ways to construct an
        # alignment.
        assert (sam_line and read_dict) or (seqan_output and read)

        # Some inputs are required for both types of construction.
        assert scoring_scheme and reference_dict

        # Read details
        self.read = None
        self.read_start_pos = None
        self.read_end_pos = None
        self.read_end_gap = None

        # Reference details
        self.ref = None
        self.ref_start_pos = None
        self.ref_end_pos = None

        # Alignment details
        self.rev_comp = None
        self.cigar_parts = None
        self.match_count = None
        self.mismatch_count = None
        self.insertion_count = None
        self.deletion_count = None
        self.alignment_length = None
        self.edit_distance = None
        self.percent_identity = None
        self.raw_score = None
        self.scaled_score = None
        self.milliseconds = None

        # How some of the values are gotten depends on whether this alignment came from SAM
        # or a Seqan alignment.
        if seqan_output:
            self.setup_using_seqan_output(seqan_output, read, reference_dict)
        elif sam_line:
            self.setup_using_sam(sam_line, read_dict, reference_dict)

        self.tally_up_score_and_errors(scoring_scheme)

    def setup_using_seqan_output(self, seqan_output, read, reference_dict):
        """
        This function sets up the Alignment using the Seqan results. This kind of alignment has
        complete details about the alignment.
        """
        seqan_parts = seqan_output.split(',', 9)
        assert len(seqan_parts) >= 10

        self.rev_comp = (seqan_parts[1] == '-')
        self.cigar_parts = re.findall(r'\d+\w', seqan_parts[9])
        self.milliseconds = int(seqan_parts[8])

        self.read = read
        self.read_start_pos = int(seqan_parts[2])
        self.read_end_pos = int(seqan_parts[3])
        self.read_end_gap = self.read.get_length() - self.read_end_pos

        self.ref = reference_dict[get_nice_header(seqan_parts[0])]
        self.ref_start_pos = int(seqan_parts[4])
        self.ref_end_pos = int(seqan_parts[5])

    def setup_using_sam(self, sam_line, read_dict, reference_dict):
        """
        This function sets up the Alignment using a SAM line.
        """
        sam_parts = sam_line.split('\t', 6)
        self.rev_comp = bool(int(sam_parts[1]) & 0x10)
        self.cigar_parts = re.findall(r'\d+\w', sam_parts[5])

        self.read = read_dict[sam_parts[0]]
        self.read_start_pos = self.get_start_soft_clips()
        self.read_end_pos = self.read.get_length() - self.get_end_soft_clips()
        self.read_end_gap = self.get_end_soft_clips()

        self.ref = reference_dict[get_nice_header(sam_parts[2])]
        self.ref_start_pos = int(sam_parts[3]) - 1
        self.ref_end_pos = self.ref_start_pos
        for cigar_part in self.cigar_parts:
            self.ref_end_pos += get_ref_shift_from_cigar_part(cigar_part)

        # If all is good with the CIGAR, then we should never end up with a ref_end_pos out of the
        # reference range. But we check just to be safe.
        if self.ref_end_pos > len(self.ref.sequence):
            self.ref_end_pos = len(self.ref.sequence)

    def tally_up_score_and_errors(self, scoring_scheme):
        """
        This function steps through the CIGAR string for the alignment to get the score, identity
        and count/locations of errors.
        """
        # Clear any existing tallies.
        self.match_count = 0
        self.mismatch_count = 0
        self.insertion_count = 0
        self.deletion_count = 0
        self.percent_identity = 0.0
        self.raw_score = 0

        # Remove the soft clipping parts of the CIGAR string for tallying.
        cigar_parts = self.cigar_parts[:]
        if cigar_parts[0][-1] == 'S':
            cigar_parts.pop(0)
        if cigar_parts and cigar_parts[-1][-1] == 'S':
            cigar_parts.pop()
        if not cigar_parts:
            return

        read_len = self.read.get_length()
        if self.rev_comp:
            read_seq = reverse_complement(self.read.sequence)
        else:
            read_seq = self.read.sequence
        read_i = self.read_start_pos

        ref_len = self.ref.get_length()
        ref_seq = self.ref.sequence
        ref_i = self.ref_start_pos
        align_i = 0

        for cigar_part in cigar_parts:
            cigar_count = int(cigar_part[:-1])
            cigar_type = cigar_part[-1]
            ins_del_cigar_score = scoring_scheme.gap_open + \
                ((cigar_count - 1) * scoring_scheme.gap_extend)
            if cigar_type == 'I':
                cigar_score = ins_del_cigar_score
                self.insertion_count += cigar_count
                read_i += cigar_count
            elif cigar_type == 'D':
                cigar_score = ins_del_cigar_score
                self.deletion_count += cigar_count
                ref_i += cigar_count
            else:  # match/mismatch
                cigar_score = 0
                for _ in range(cigar_count):
                    # If all is good with the CIGAR, then we should never end up with a sequence
                    # index out of the sequence range. But a CIGAR error can cause this, so check
                    # here.
                    if read_i >= read_len or ref_i >= ref_len:
                        break
                    read_base = read_seq[read_i]
                    ref_base = ref_seq[ref_i]
                    if read_base == ref_base:
                        self.match_count += 1
                        cigar_score += scoring_scheme.match
                    else:
                        self.mismatch_count += 1
                        cigar_score += scoring_scheme.mismatch
                    read_i += 1
                    ref_i += 1

            self.raw_score += cigar_score
            align_i += cigar_count

        self.percent_identity = 100.0 * self.match_count / align_i
        self.edit_distance = self.mismatch_count + self.insertion_count + self.deletion_count
        self.alignment_length = align_i
        perfect_score = scoring_scheme.match * self.alignment_length
        worst_score = scoring_scheme.mismatch * self.alignment_length
        self.scaled_score = 100.0 * (self.raw_score - worst_score) / (perfect_score - worst_score)

    def __repr__(self):
        read_start, read_end = self.read_start_end_positive_strand()
        return_str = self.read.name + ' (' + str(read_start) + '-' + str(read_end) + ', '
        if self.rev_comp:
            return_str += 'strand: -), '
        else:
            return_str += 'strand: +), '
        return_str += self.ref.name + ' (' + str(self.ref_start_pos) + '-' + \
            str(self.ref_end_pos) + ')'
        if self.scaled_score is not None:
            return_str += ', raw score = ' + str(self.raw_score)
            return_str += ', scaled score = ' + float_to_str(self.scaled_score, 2)
        if self.percent_identity is not None:
            return_str += ', ' + float_to_str(self.percent_identity, 2) + '% ID'
        return return_str

    def get_aligned_ref_length(self):
        """
        Returns the length of the reference used in this alignment. Could be the whole reference
        length or just a part of it.
        """
        return self.ref_end_pos - self.ref_start_pos

    def get_aligned_read_length(self):
        """
        Returns the length of the read used in this alignment. Could be the whole read length or
        just a part of it.
        """
        return self.read_end_pos - self.read_start_pos

    def get_ref_to_read_ratio(self):
        """
        Returns the length ratio between the aligned parts of the reference and read.
        """
        return self.get_aligned_ref_length() / self.get_aligned_read_length()

    def get_read_to_ref_ratio(self):
        """
        Returns the length ratio between the aligned parts of the read and reference.
        """
        return 1.0 / self.get_ref_to_read_ratio()

    def read_start_end_positive_strand(self):
        """
        This function returns the read start/end coordinates for the positive strand of the read.
        For alignments on the positive strand, this is just the normal start/end. But for
        alignments on the negative strand, the coordinates are flipped to the other side.
        """
        return self.read_start_positive_strand(), self.read_end_positive_strand()

    def read_start_positive_strand(self):
        """
        This function returns the read start coordinates for the positive strand of the read.
        """
        if self.rev_comp:
            return self.read.get_length() - self.read_end_pos
        else:
            return self.read_start_pos

    def read_end_positive_strand(self):
        """
        This function returns the read end coordinates for the positive strand of the read.
        """
        if self.rev_comp:
            return self.read.get_length() - self.read_start_pos
        else:
            return self.read_end_pos

    def get_start_soft_clips(self):
        """
        Returns the number of soft-clipped bases at the start of the alignment.
        """
        if self.cigar_parts[0][-1] == 'S':
            return int(self.cigar_parts[0][:-1])
        else:
            return 0

    def get_end_soft_clips(self):
        """
        Returns the number of soft-clipped bases at the start of the alignment.
        """
        if self.cigar_parts[-1][-1] == 'S':
            return int(self.cigar_parts[-1][:-1])
        else:
            return 0

    def get_sam_line(self):
        """
        Returns a SAM alignment line.
        """
        sam_parts = [self.read.name]  # Query template name
        if self.rev_comp:
            sam_parts.append('16')  # Bitwise flag
        else:
            sam_parts.append('0')  # Bitwise flag
        sam_parts.append(self.ref.name)  # Reference sequence name
        sam_parts.append(str(self.ref_start_pos + 1))  # 1-based leftmost mapping position
        sam_parts.append('255')  # Mapping quality (255 means unavailable)
        sam_parts.append(''.join(self.cigar_parts))  # CIGAR string
        sam_parts.append('*')  # Ref. name of the mate/next read (* means unavailable)
        sam_parts.append('0')  # Position of the mate/next read (0 means unavailable)
        sam_parts.append('0')  # Observed template length (0 means unavailable)

        if self.rev_comp:
            sam_parts.append(reverse_complement(self.read.sequence))  # Segment sequence
            sam_parts.append(self.read.qualities[::-1])  # ASCII of Phred-scaled base quality+33
        else:
            sam_parts.append(self.read.sequence)  # Segment sequence
            sam_parts.append(self.read.qualities)  # ASCII of Phred-scaled base quality+33

        sam_parts.append('AS:i:' + str(self.raw_score))  # Alignment score generated by aligner

        edit_distance = self.mismatch_count + self.insertion_count + self.deletion_count
        sam_parts.append('NM:i:' + str(edit_distance))  # Edit distance to the reference, including
        # ambiguous bases but excluding clipping
        return '\t'.join(sam_parts) + '\n'

    def is_very_similar(self, other):
        """
        Returns true if this alignment and the other alignment seem to be redundant.
        Specifically, the have to be from the same read, the same reference and overlap by 90% or
        more.
        """
        if self.read.name != other.read.name:
            return False
        if self.ref.name != other.ref.name:
            return False
        if self.rev_comp != other.rev_comp:
            return False

        this_start, this_end = self.read_start_end_positive_strand()
        other_start, other_end = other.read_start_end_positive_strand()
        if other_start > this_end or this_start > other_end:
            return False

        # If the code got here then the alignments are overlapping.
        overlap_size = min(this_end, other_end) - max(this_start, other_start)
        smaller_alignment_length = min(this_end - this_start, other_end - other_start)
        if smaller_alignment_length == 0:
            return False
        return overlap_size / smaller_alignment_length >= 0.9

    def get_signed_ref_num(self):
        """
        If the reference is in SPAdes contig format, then this function returns the number of the
        contig with the sign from the alignment.
        """
        if self.rev_comp:
            return -self.ref.number
        else:
            return self.ref.number


def get_ref_shift_from_cigar_part(cigar_part):
    """
    This function returns how much a given cigar moves on a reference.
    Examples:
      * '5M' returns 5
      * '5S' returns 0
      * '5D' returns 5
      * '5I' returns 0
    """
    if cigar_part[-1] == 'M':
        return int(cigar_part[:-1])
    if cigar_part[-1] == 'I':
        return 0
    if cigar_part[-1] == 'D':
        return int(cigar_part[:-1])
    if cigar_part[-1] == 'S':
        return 0
