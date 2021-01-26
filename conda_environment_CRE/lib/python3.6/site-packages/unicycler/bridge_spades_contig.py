"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

SPAdes contig bridges are bridges created from the contigs.paths file made by SPAdes.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import math
from .bridge_common import get_bridge_str, get_mean_depth, get_depth_agreement_factor
from .misc import float_to_str, get_num_agreement, get_right_arrow, print_table
from . import log


class SpadesContigBridge(object):
    """
    This class describes a bridge created from the contigs.paths file made by SPAdes.

    Quality is affected by:
      * How well the start and end segments' depths agree.
      * The depth consistency within the path (only applies to bridges where the path segments
        exclusively lead to the start/end segments).
    """
    def __init__(self, graph, spades_contig_path):

        # The numbers of the two single copy segments which are being bridged.
        self.start_segment = None
        self.end_segment = None

        # The path through the unbridged graph.
        self.graph_path = []

        # The bridge sequence, gotten from the graph path.
        self.bridge_sequence = ''

        # The bridge depth, a weighted mean of the start and end depths.
        self.depth = 0.0

        # A score used to determine the order of bridge application. SPAdes contig bridges don't
        # start at 1.0 and so cannot get as high as long read bridges potentially can.
        self.quality = 0.4

        # When a bridge is applied, the segments in the bridge may have their depth reduced
        # accordingly. This member stores which segments have had their depth reduced and by how
        # much due to this bridge's application. It is stored so if this bridge is later deleted,
        # we can restore the depth to the segments.
        self.segments_reduced_depth = []

        # The first and last values in spades_contig_path are the start and end segments. The
        # values in between are the path.
        self.graph_path = spades_contig_path
        self.start_segment = self.graph_path.pop(0)
        self.end_segment = self.graph_path.pop()

        # If there are segments in between the start and end (there usually will be), then they
        # provide the bridge sequence.
        self.bridge_sequence = graph.get_path_sequence(self.graph_path)

        # The start segment and end segment should agree in depth. If they don't, that's bad.
        start_seg = graph.segments[abs(self.start_segment)]
        end_seg = graph.segments[abs(self.end_segment)]
        self.quality *= get_depth_agreement_factor(start_seg.depth, end_seg.depth)

        # If the segments in the path exclusively lead to the start and end segments (i.e they
        # cannot lead to any another segment), then we can also scale the quality based on the
        # depth consistency of the path. E.g. if a bridge path contains a segment 3 times and that
        # segment's depth also suggests about 3 times, that's good. If they don't agree, that's
        # bad.
        self.depth = get_mean_depth(start_seg, end_seg, graph)
        if path_is_self_contained(self.graph_path, self.start_segment, self.end_segment, graph):
            graph_path_pos_nums = list(set([abs(x) for x in self.graph_path]))
            for path_segment in graph_path_pos_nums:
                actual_depth = graph.segments[path_segment].depth
                expected_depth = graph_path_pos_nums.count(path_segment) * self.depth
                agreement = get_num_agreement(actual_depth, expected_depth)
                self.quality *= agreement

        # The quality of a SPAdes contig bridge should decline sharply if the bridging sequence
        # is too large compared to the short read insert size. E.g. if the insert size is 500 bp
        # and the SPAdes contig bridge is 2 kb long, we should not believe it!
        if self.graph_path:
            bridge_length = len(self.bridge_sequence)
            if bridge_length <= graph.insert_size_mean:
                bridge_length_factor = 1.0
            else:
                bridge_length_factor = graph.insert_size_deviation / (bridge_length -
                                                                      graph.insert_size_mean +
                                                                      graph.insert_size_deviation)
            self.quality *= bridge_length_factor

        # We finalise the quality to a range of 0 to 100. We also use the sqrt function to pull
        # the scores up a bit (otherwise they tend to hang near the bottom of the range).
        self.quality = 100.0 * math.sqrt(self.quality)

    def __repr__(self):
        return 'SPAdes bridge: ' + get_bridge_str(self) + \
               ' (quality = ' + float_to_str(self.quality, 2) + ')'

    @staticmethod
    def get_type_score():
        """
        Returns a score indicating the relative importance of the bridge types:
        LongReadBridge = 2, SpadesContigBridge = 1, LoopUnrollingBridge = 0
        """
        return 1

    @staticmethod
    def get_type_name():
        """
        Returns the name of the bridge type.
        """
        return 'SPAdes'


def create_spades_contig_bridges(graph, anchor_segments):
    """
    Builds graph bridges using the SPAdes contig paths.
    """
    log.log_section_header('Creating SPAdes contig bridges')
    log.log_explanation('SPAdes uses paired-end information to perform repeat resolution (RR) and '
                        'produce contigs from the assembly graph. SPAdes saves the graph paths '
                        'corresponding to these contigs in the contigs.paths file. When one of '
                        'these paths contains two or more anchor contigs, Unicycler can '
                        'create a bridge from the path.', verbosity=1)

    bridge_path_set = set()
    single_copy_numbers = [x.number for x in anchor_segments]
    for segment in anchor_segments:
        for path in graph.paths.values():
            flipped_path = [-x for x in reversed(path)]
            contig_bridges = find_contig_bridges(segment.number, path, single_copy_numbers)
            contig_bridges += find_contig_bridges(segment.number, flipped_path, single_copy_numbers)
            for contig_bridge in contig_bridges:
                flipped_contig_bridge = [-x for x in reversed(contig_bridge)]
                contig_bridge_str = ','.join([str(x) for x in contig_bridge])
                flipped_contig_bridge_str = ','.join([str(x) for x in flipped_contig_bridge])
                if contig_bridge_str not in bridge_path_set and \
                        flipped_contig_bridge_str not in bridge_path_set:
                    if contig_bridge[0] < 0 and contig_bridge[-1] < 0:
                        bridge_path_set.add(flipped_contig_bridge_str)
                    else:
                        bridge_path_set.add(contig_bridge_str)

    bridge_path_list = sorted(list([[int(y) for y in x.split(',')] for x in bridge_path_set]))

    # If multiple bridge paths start with or end with the same segment, that implies a conflict
    # between SPADes' paths and our single copy determination. Throw these bridges out.
    bridge_paths_by_start = {}
    bridge_paths_by_end = {}
    for path in bridge_path_list:
        start = path[0]
        end = path[-1]
        if start not in bridge_paths_by_start:
            bridge_paths_by_start[start] = []
        if end not in bridge_paths_by_end:
            bridge_paths_by_end[end] = []
        if -end not in bridge_paths_by_start:
            bridge_paths_by_start[-end] = []
        if -start not in bridge_paths_by_end:
            bridge_paths_by_end[-start] = []
        bridge_paths_by_start[start].append(path)
        bridge_paths_by_end[end].append(path)
        bridge_paths_by_start[-end].append(path)
        bridge_paths_by_end[-start].append(path)
    conflicting_paths = []
    for grouped_paths in bridge_paths_by_start.values():
        if len(grouped_paths) > 1:
            conflicting_paths += grouped_paths
    for grouped_paths in bridge_paths_by_end.values():
        if len(grouped_paths) > 1:
            conflicting_paths += grouped_paths
    conflicting_paths_no_dupes = []
    for path in conflicting_paths:
        if path not in conflicting_paths_no_dupes:
            conflicting_paths_no_dupes.append(path)
    conflicting_paths = conflicting_paths_no_dupes
    final_bridge_paths = [x for x in bridge_path_list if x not in conflicting_paths]

    bridges = [SpadesContigBridge(spades_contig_path=x, graph=graph) for x in final_bridge_paths]

    if bridges:
        bridge_table = [['Start', 'Path', 'End', 'Bridge quality']]
        spaded_arrow = ' ' + get_right_arrow() + ' '
        max_path_str_len = 5
        for bridge in bridges:
            path_str = spaded_arrow.join(str(seg) for seg in bridge.graph_path)
            max_path_str_len = max(max_path_str_len, len(path_str))
            bridge_table.append([str(bridge.start_segment), path_str, str(bridge.end_segment),
                                 float_to_str(bridge.quality, 1)])
        print_table(bridge_table, alignments='RCLR', left_align_header=False, indent=0,
                    fixed_col_widths=[5, max_path_str_len, 5, 7])
    else:
        log.log('No SPAdes contig bridges')

    return bridges


def find_contig_bridges(segment_num, path, single_copy_numbers):
    """
    This function returns a list of lists: every part of the path which starts on the segment_num
    and ends on any of the single_copy_numbers.
    """
    bridge_paths = []
    indices = [i for i, x in enumerate(path) if abs(x) == segment_num]
    for index in indices:
        bridge_path = [path[index]]
        for i in range(index + 1, len(path)):
            bridge_path.append(path[i])
            if path[i] in single_copy_numbers or -path[i] in single_copy_numbers:
                break
        else:
            bridge_path = []
        if bridge_path:
            bridge_paths.append(bridge_path)
    return bridge_paths


def path_is_self_contained(path, start, end, graph):
    """
    Returns True if the path segments are only connected to each other and the start/end segments.
    If they are connected to anything else, it returns False.
    """
    all_numbers_in_path = set()
    all_numbers_in_path.add(abs(start))
    all_numbers_in_path.add(abs(end))
    for segment in path:
        all_numbers_in_path.add(abs(segment))
    for segment in path:
        connected_segments = graph.get_connected_segments(segment)
        for connected_segment in connected_segments:
            if connected_segment not in all_numbers_in_path:
                return False
    return True
