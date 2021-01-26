"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

Long read bridges are the big important type of bridge for a hybrid Unicycler assembly. They are
made using long reads which align to multiple segments in the graph.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import time
import math
from .bridge_common import get_bridge_str, get_mean_depth, get_depth_agreement_factor, \
    get_bridge_table_parameters, print_bridge_table_header, print_bridge_table_row
from .misc import float_to_str
from . import settings
from . import log
from .path_finding import get_best_paths_for_seq


class MiniasmBridge(object):
    """
    This class describes a bridge created from long read alignments.
    """
    def __init__(self, graph, start, end, bridge_sequence, start_overlap, end_overlap,
                 scoring_scheme, output, do_path_search=True):

        # The numbers of the two single copy segments which are being bridged.
        self.start_segment = start
        self.end_segment = end
        output += [str(start), str(end), '']

        # In some miniasm bridges where the contigs are very close, there will be overlap between
        # the contigs and the bridge. I.e. when the bridge is applied, contig sequence will be
        # trimmed off.
        self.start_overlap = start_overlap
        self.end_overlap = end_overlap

        # The bridge depth, a weighted mean of the start and end depths.
        self.depth = get_mean_depth(graph.segments[abs(self.start_segment)],
                                    graph.segments[abs(self.end_segment)], graph)

        self.segments_reduced_depth = []

        # When a MiniasmBridge is being made by splitting an already-made bridge, then we don't
        # do the path search.
        if not do_path_search:
            self.bridge_sequence = bridge_sequence
            self.all_paths = []
            self.graph_path = []
            self.quality = 1.0

        else:
            # Look for a graph path corresponding to the bridge sequence.
            target_path_length = len(bridge_sequence)
            output += [str(target_path_length), '', str(target_path_length)]
            path_start_time = time.time()
            self.all_paths, progressive_path_search = \
                get_best_paths_for_seq(graph, self.start_segment, self.end_segment,
                                       target_path_length, bridge_sequence, scoring_scheme, 90.0)
            path_time = time.time() - path_start_time

            output.append(str(len(self.all_paths)))
            output.append('progressive' if progressive_path_search else 'exhaustive')
            output.append(float_to_str(path_time, 1))

            if self.all_paths:
                self.graph_path = self.all_paths[0][0]
                raw_score = self.all_paths[0][1]
                scaled_score = self.all_paths[0][3]
                len_discrepancy = self.all_paths[0][2]
                if self.graph_path:
                    output.append(', '.join(str(x) for x in self.graph_path))
                else:
                    output.append('direct connection')
                best_path_len = graph.get_bridge_path_length(self.graph_path)
                output.append(str(best_path_len))
                output.append(float_to_str(raw_score, 1))
                output.append(float_to_str(scaled_score, 2))
                output.append(str(len_discrepancy))

            else:
                self.graph_path = []
                output += ['', '', '', '', '']
                scaled_score = 0.0

            # If a very good match was found, use the graph path and give the bridge a very high
            # quality score.
            if scaled_score > settings.MINIASM_BRIDGE_SCALED_SCORE_TO_USE_GRAPH_PATH:
                self.bridge_sequence = graph.get_path_sequence(self.graph_path)
                self.quality = settings.MINIASM_BRIDGE_QUAL_WITH_GRAPH_PATH

            # Otherwise, just use the miniasm bridge sequence for the bridge.
            else:
                self.bridge_sequence = bridge_sequence
                if graph.ends_with_dead_end(self.start_segment) or \
                        graph.starts_with_dead_end(self.end_segment):
                    self.quality = settings.MINIASM_BRIDGE_QUAL_WITH_DEAD_END
                else:
                    self.quality = settings.MINIASM_BRIDGE_QUAL_WITHOUT_PATH_OR_DEAD_END

            # Depth agreement affects bridge quality.
            start_seg = graph.segments[abs(self.start_segment)]
            end_seg = graph.segments[abs(self.end_segment)]
            self.quality *= get_depth_agreement_factor(start_seg.depth, end_seg.depth)

            # Bridge length affects quality too: short bridges are better.
            bridge_len = max(0, len(self.bridge_sequence))
            half_qual_len = settings.MINIASM_BRIDGE_HALF_QUAL_LENGTH
            self.quality *= half_qual_len / (bridge_len + half_qual_len)

            self.quality = 100.0 * math.sqrt(self.quality)
            output.append(self.quality)

    def __repr__(self):
        return 'miniasm bridge: ' + get_bridge_str(self) + \
               ' (quality = ' + float_to_str(self.quality, 2) + ')'

    @staticmethod
    def get_type_score():
        """
        Returns a score indicating the relative importance of the bridge types:
        LongReadBridge = 2, SpadesContigBridge = 1, LoopUnrollingBridge = 0
        """
        return 2

    @staticmethod
    def get_type_name():
        """
        Returns the name of the bridge type.
        """
        return 'miniasm'


def create_miniasm_bridges(graph, string_graph, anchor_segments, scoring_scheme, verbosity,
                           min_bridge_qual):
    """
    Makes bridges between single copy segments using the miniasm string graph.
    """
    log.log_section_header('Creating miniasm/Racon bridges')
    log.log_explanation('Now that the miniasm/Racon string graph is complete, Unicycler will '
                        'use it to build bridges between anchor segments.', verbosity=1)
    bridges = []
    anchor_seg_nums = set(x.number for x in anchor_segments)

    string_graph_bridge_segments = sorted([x for x in string_graph.segments
                                           if x.startswith('BRIDGE_') or
                                           x.startswith('OVERLAPPING_BRIDGE_')])

    filtered_string_graph_bridge_segments = []
    for bridge_seg_name in string_graph_bridge_segments:
        pos_seg_name = bridge_seg_name + '+'
        preceding_segments = string_graph.get_preceding_segments(pos_seg_name)
        following_segments = string_graph.get_following_segments(pos_seg_name)
        if len(preceding_segments) != 1:
            continue
        if len(following_segments) != 1:
            continue
        preceding_seg_name = preceding_segments[0]
        following_seg_name = following_segments[0]
        if not preceding_seg_name.startswith('CONTIG_') or \
                not following_seg_name.startswith('CONTIG_'):
            continue
        filtered_string_graph_bridge_segments.append(bridge_seg_name)

    # We want to display this table one row at a time, so we have to fix all of the column widths
    # at the start.
    bridge_count = len(filtered_string_graph_bridge_segments)
    alignments, col_widths = get_bridge_table_parameters(graph, bridge_count, verbosity,
                                                         'MiniasmBridge')
    print_bridge_table_header(alignments, col_widths, verbosity, 'MiniasmBridge')
    completed_count = 0

    for bridge_seg_name in filtered_string_graph_bridge_segments:
        bridge_seg = string_graph.segments[bridge_seg_name]
        pos_seg_name = bridge_seg_name + '+'
        preceding_seg_name = string_graph.get_preceding_segments(pos_seg_name)[0]
        following_seg_name = string_graph.get_following_segments(pos_seg_name)[0]

        first_link = string_graph.links[(preceding_seg_name, pos_seg_name)]
        second_link = string_graph.links[(pos_seg_name, following_seg_name)]
        preceding_seg_name = preceding_seg_name[7:]
        following_seg_name = following_seg_name[7:]
        preceding_segment_number = int(preceding_seg_name[:-1]) * \
            (1 if preceding_seg_name[-1] == '+' else -1)
        following_segment_number = int(following_seg_name[:-1]) * \
            (1 if following_seg_name[-1] == '+' else -1)
        assert abs(preceding_segment_number) in anchor_seg_nums
        assert abs(following_segment_number) in anchor_seg_nums

        start_overlap = first_link.seg_1_overlap
        end_overlap = second_link.seg_2_overlap
        output = []
        bridge = MiniasmBridge(graph, preceding_segment_number, following_segment_number,
                               bridge_seg.forward_sequence, start_overlap, end_overlap,
                               scoring_scheme, output)
        bridges.append(bridge)

        completed_count += 1
        print_bridge_table_row(alignments, col_widths, output, completed_count,
                               bridge_count, min_bridge_qual, verbosity, 'MiniasmBridge')

    # Now that the bridges are finalised, we split bridges that contain anchor segments in their
    # path such that all bridges start and end on an anchor segment but contain no anchor segments
    # in their path.
    split_bridges = []
    for bridge in bridges:
        if not bridge.graph_path:
            split_bridges.append(bridge)
        else:  # bridge has a path
            if not any(abs(x) in anchor_seg_nums for x in bridge.graph_path):  # already good
                split_bridges.append(bridge)

            # If the bridge path contains one or more anchor segments, it must be split!
            else:

                # A bridge with a path shouldn't have overlaps.
                assert bridge.start_overlap == 0
                assert bridge.end_overlap == 0

                full_path = [bridge.start_segment] + bridge.graph_path + [bridge.end_segment]
                anchor_indices = []
                for i, seg_num in enumerate(full_path):
                    if abs(seg_num) in anchor_seg_nums:
                        anchor_indices.append(i)
                anchor_indices = sorted(anchor_indices)
                for i in range(len(anchor_indices) - 1):
                    start_i = anchor_indices[i]
                    end_i = anchor_indices[i+1]
                    start_seg_num = full_path[start_i]
                    end_seg_num = full_path[end_i]
                    new_path = full_path[start_i+1:end_i]
                    bridge_sequence = graph.get_path_sequence(new_path)
                    split_bridge = MiniasmBridge(graph, start_seg_num, end_seg_num, bridge_sequence,
                                                 0, 0, scoring_scheme, [])
                    split_bridge.graph_path = new_path
                    split_bridge.all_paths = [new_path]
                    split_bridge.quality = bridge.quality
                    split_bridges.append(split_bridge)

    return bridges
