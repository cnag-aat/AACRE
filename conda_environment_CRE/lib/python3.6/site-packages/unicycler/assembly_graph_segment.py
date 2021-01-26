"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module describes the segments (a.k.a. contigs, a.k.a. nodes) which hold the sequences in an
assembly graph.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import textwrap
from .misc import reverse_complement, add_line_breaks_to_sequence
from .bridge_long_read import LongReadBridge
from .bridge_spades_contig import SpadesContigBridge
from .bridge_loop_unroll import LoopUnrollingBridge
from .bridge_long_read_simple import SimpleLongReadBridge
from .bridge_miniasm import MiniasmBridge


# noinspection PyAugmentAssignment
class Segment(object):
    """
    This hold a graph segment with a number, depth, direction and sequence.
    """
    def __init__(self, number, depth, sequence, positive, bridge=None, graph_path=None,
                 original_depth=True):
        self.number = number
        self.depth = depth
        self.original_depth = original_depth
        self.forward_sequence = ''
        self.reverse_sequence = ''
        self.bridge = bridge
        self.graph_path = graph_path
        if positive:
            self.forward_sequence = sequence
        else:
            self.reverse_sequence = sequence
        self.used_in_bridges = []

    def __repr__(self):
        if len(self.forward_sequence) > 6:
            seq_string = self.forward_sequence[:3] + '...' + self.forward_sequence[-3:]
        else:
            seq_string = self.forward_sequence
        return str(self.number) + ' (' + seq_string + ')'

    def add_sequence(self, sequence, positive):
        if positive:
            self.forward_sequence = sequence
        else:
            self.reverse_sequence = sequence

    def build_other_sequence_if_necessary(self):
        if not self.forward_sequence:
            self.forward_sequence = reverse_complement(self.reverse_sequence)
        if not self.reverse_sequence:
            self.reverse_sequence = reverse_complement(self.forward_sequence)

    def get_length(self):
        return len(self.forward_sequence)

    def get_length_no_overlap(self, overlap):
        return len(self.forward_sequence) - overlap

    def is_homopolymer(self):
        """
        Returns True if the segment's sequence is made up of only one base.
        """
        if len(self.forward_sequence) == 0:
            return False
        first_base = self.forward_sequence[0].lower()
        for base in self.forward_sequence[1:]:
            if base.lower() != first_base:
                return False
        return True

    def gfa_segment_line(self):
        """
        Returns an entire S line for GFA output, including the newline.
        """
        s_line = 'S\t'
        s_line += str(self.number) + '\t'
        s_line += self.forward_sequence + '\t'
        s_line += 'LN:i:' + str(self.get_length()) + '\t'
        s_line += 'dp:f:' + str(self.depth) + '\n'
        return s_line

    def get_fasta_name_and_description_line(self, circular_seg_nums=None):
        """
        Returns the segment's fasta line, including the '>' and the newline.
        """
        name_and_description = ''.join(['>', str(self.number), ' length=', str(self.get_length()),
                                        ' depth=', '%.2f' % self.depth, 'x'])
        if circular_seg_nums and self.number in circular_seg_nums:
            name_and_description += ' circular=true'
        return name_and_description + '\n'

    def save_to_fasta(self, fasta_filename):
        """
        Saves the segment's sequence to FASTA file.
        """
        fasta = open(fasta_filename, 'w')
        fasta.write(self.get_fasta_name_and_description_line())
        fasta.write(add_line_breaks_to_sequence(self.forward_sequence))
        fasta.close()

    def get_seg_type_label(self):
        """
        Given a particular segment, this function returns a label string based its type.
        """
        if self.bridge is None:
            return ''
        if isinstance(self.bridge, SpadesContigBridge):
            label = 'SPAdes contig bridge'
        elif isinstance(self.bridge, LoopUnrollingBridge):
            label = 'Loop unrolling bridge'
        elif isinstance(self.bridge, LongReadBridge):
            label = 'Long read bridge'
        elif isinstance(self.bridge, SimpleLongReadBridge):
            label = 'Simple long read bridge'
        elif isinstance(self.bridge, MiniasmBridge):
            label = 'Miniasm bridge'
        else:
            raise TypeError("unknown bridge type")
        if self.graph_path:
            graph_path_str = ', '.join([str(x) for x in self.graph_path])
            graph_path_str = '\\n'.join(textwrap.wrap(graph_path_str, 40))
            label += ':\\n' + graph_path_str
        return label

    def trim_from_end(self, amount):
        """
        Removes the specified number of bases from the end of the segment sequence.
        """
        assert self.get_length() >= amount
        if amount == 0:
            return
        self.forward_sequence = self.forward_sequence[:-amount]
        self.reverse_sequence = self.reverse_sequence[amount:]

    def trim_from_start(self, amount):
        """
        Removes the specified number of bases from the end of the segment sequence.
        """
        assert self.get_length() >= amount
        if amount == 0:
            return
        self.forward_sequence = self.forward_sequence[amount:]
        self.reverse_sequence = self.reverse_sequence[:-amount]

    def append_to_forward_sequence(self, additional_seq):
        """
        Adds the given sequence to the end of the forward sequence (and updates the reverse
        sequence accordingly).
        """
        self.forward_sequence = self.forward_sequence + additional_seq
        self.reverse_sequence = reverse_complement(self.forward_sequence)

    def append_to_reverse_sequence(self, additional_seq):
        """
        Adds the given sequence to the end of the reverse sequence (and updates the forward
        sequence accordingly).
        """
        self.reverse_sequence = self.reverse_sequence + additional_seq
        self.forward_sequence = reverse_complement(self.reverse_sequence)

    def prepend_to_forward_sequence(self, additional_seq):
        """
        Adds the given sequence to the end of the forward sequence (and updates the reverse
        sequence accordingly).
        """
        self.forward_sequence = additional_seq + self.forward_sequence
        self.reverse_sequence = reverse_complement(self.forward_sequence)

    def prepend_to_reverse_sequence(self, additional_seq):
        """
        Adds the given sequence to the end of the reverse sequence (and updates the forward
        sequence accordingly).
        """
        self.reverse_sequence = additional_seq + self.reverse_sequence
        self.forward_sequence = reverse_complement(self.reverse_sequence)

    def remove_sequence(self):
        """
        Gets rid of the segment sequence entirely, turning it into a zero-length segment.
        """
        self.forward_sequence = ''
        self.reverse_sequence = ''

    def rotate_sequence(self, start_pos, flip):
        """
        Rotates the sequence so it begins at start_pos. If flip is True, it also switches the
        forward and reverse strands. This function assumes that the segment is a circular
        completed replicon with no overlap.
        """
        unrotated_seq = self.forward_sequence
        rotated_seq = unrotated_seq[start_pos:] + unrotated_seq[:start_pos]
        rev_comp_rotated_seq = reverse_complement(rotated_seq)

        if flip:
            self.forward_sequence = rev_comp_rotated_seq
            self.reverse_sequence = rotated_seq
        else:
            self.forward_sequence = rotated_seq
            self.reverse_sequence = rev_comp_rotated_seq
