"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

Loop unrolling bridges strive to unroll simple loops in the graph. They require a SPAdes contig
that connects the loop sequence with the base sequence, in order to be reasonably sure that we are
dealing with a loop and not a separate circular sequence. They use the relative depths of the
segments to determine how many times to traverse the loop.

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
from .misc import float_to_str, weighted_average, print_table
from . import log


class LoopUnrollingBridge(object):
    """
    This class describes a bridge created from unrolling an assembly graph loop.

    Quality is affected by:
      * How well the start and end segments' depths agree.
      * How close the determined loop count is to a whole number.
      * The final loop count (higher counts get lower quality).
    """
    def __init__(self, graph, start, end, middle, repeat):
        """
        This constructor assumes the the start, end, middle and repeat segments form a simple loop
        in the graph supported by either a SPAdes contig or a long read alignment. It will use
        segment depths to determine the loop count and score the bridge's quality.
        """
        # The numbers of the two single copy segments which are being bridged.
        self.start_segment = start
        self.end_segment = end

        self.middle_segment = middle
        self.repeat_segment = repeat

        # The path through the unbridged graph.
        self.graph_path = []

        # The bridge sequence, gotten from the graph path.
        self.bridge_sequence = ''

        # The bridge depth, a weighted mean of the start and end depths.
        self.depth = 0.0

        # A score used to determine the order of bridge application. This value starts at the
        # maximum for a loop unrolling bridge and can only decrease as the constructor continues.
        self.quality = 0.2

        # When a bridge is applied, the segments in the bridge may have their depth reduced
        # accordingly. This member stores which segments have had their depth reduced and by how
        # much due to this bridge's application. It is stored so if this bridge is later deleted,
        # we can restore the depth to the segments.
        self.segments_reduced_depth = []

        # Get the actual segments from the numbers. Since we are assuming they do form a simple
        # loop, we don't care about directionality.
        start_seg = graph.segments[abs(start)]
        end_seg = graph.segments[abs(end)]
        middle_seg = graph.segments[abs(middle)]
        repeat_seg = graph.segments[abs(repeat)]

        # The start segment and end segment should agree in depth. If they don't, that's very bad,
        # so depth_disagreement is applied to quality twice (squared effect).
        self.quality *= get_depth_agreement_factor(start_seg.depth, end_seg.depth)

        # We'll use a mean loop count that's weighted by the middle and repeat segment lengths.
        self.depth = get_mean_depth(start_seg, end_seg, graph)
        self.loop_count_by_middle = middle_seg.depth / self.depth
        self.loop_count_by_repeat = max((repeat_seg.depth - self.depth) / self.depth, 0.0)
        mean_loop_count = weighted_average(self.loop_count_by_middle, self.loop_count_by_repeat,
                                           middle_seg.get_length_no_overlap(graph.overlap),
                                           repeat_seg.get_length_no_overlap(graph.overlap))

        # If the average loop count is near a whole number, that's better. If it's near 0.5, that's
        # very bad because we don't know whether to round up or down.
        if mean_loop_count < 1.0:
            self.loop_count = 1
            closeness_to_whole_num = mean_loop_count
        else:
            self.loop_count = int(round(mean_loop_count))
            fractional_part = mean_loop_count % 1
            distance_from_whole_num = min(fractional_part, 1.0 - fractional_part)
            closeness_to_whole_num = 1.0 - (2.0 * distance_from_whole_num)
        self.quality *= closeness_to_whole_num

        # Finally, we reduce the quality for higher loop counts, as those are harder to call.
        loop_count_penalty = 1 / (2 ** (self.loop_count - 1))
        self.quality *= loop_count_penalty

        self.graph_path = [repeat]
        for _ in range(self.loop_count):
            self.graph_path += [middle, repeat]
        self.bridge_sequence = graph.get_path_sequence(self.graph_path)

        # We finalise the quality to a range of 0 to 100. We also use the sqrt function to pull
        # the scores up a bit (otherwise they tend to hang near the bottom of the range).
        self.quality = 100.0 * math.sqrt(self.quality)

    def __repr__(self):
        return 'loop bridge: ' + get_bridge_str(self) + \
               ' (quality = ' + float_to_str(self.quality, 2) + ')'

    @staticmethod
    def get_type_score():
        """
        Returns a score indicating the relative importance of the bridge types:
        LongReadBridge = 2, SpadesContigBridge = 1, LoopUnrollingBridge = 0
        """
        return 0

    @staticmethod
    def get_type_name():
        """
        Returns the name of the bridge type.
        """
        return 'loop'


def create_loop_unrolling_bridges(graph, anchor_segments):
    """
    This function creates loop unrolling bridges using the information in SPAdes paths.
    """
    log.log_section_header('Creating loop unrolling bridges')
    log.log_explanation('When a SPAdes contig path connects an anchor contig with the '
                        'middle contig of a simple loop, Unicycler concludes that the sequences '
                        'are contiguous (i.e. the loop is not a separate piece of DNA). It then '
                        'uses the read depth of the middle and repeat contigs to guess the number '
                        'of times to traverse the loop and makes a bridge.', verbosity=1)

    bridges = []
    simple_loops = graph.find_all_simple_loops()
    simple_loops = [x for x in simple_loops if x[2] is not None]
    seg_nums_to_bridge = set(x.number for x in anchor_segments)

    # A simple loop can either be caused by a repeat in one sequence (probably more typical) or by
    # a separate circular sequence which has some common sequence (less typical, but still very
    # possible: plasmids). We only want to unroll the former group, so we look for cases where the
    # loop's start or end is in a SPAdes contig path along with the middle. That implies that they
    # are on the same piece of DNA and can be unrolled.
    for start, end, middle, repeat in simple_loops:

        # We only want where the start and end are to-be-bridged segments but the middle is not.
        if abs(start) not in seg_nums_to_bridge:
            continue
        if abs(end) not in seg_nums_to_bridge:
            continue
        if abs(repeat) in seg_nums_to_bridge:
            continue

        joined = False
        for path in graph.paths.values():
            flipped_path = [-x for x in reversed(path)]
            if (start in path and middle in path) or \
                    (end in path and middle in path) or \
                    (start in flipped_path and middle in flipped_path) or \
                    (end in flipped_path and middle in flipped_path):
                joined = True
                break

        # If we've found evidence the simply loop is a single piece of DNA, then we'll make a loop
        # unrolling bridge!
        if joined:
            bridges.append(LoopUnrollingBridge(graph, start, end, middle, repeat))

    if bridges:
        bridge_table = [['Start', 'Repeat', 'Middle', 'End', 'Loop count by repeat',
                         'Loop count by middle', 'Loop count', 'Bridge quality']]
        for bridge in bridges:
            bridge_table.append([str(bridge.start_segment), str(bridge.repeat_segment),
                                 str(bridge.middle_segment), str(bridge.end_segment),
                                 float_to_str(bridge.loop_count_by_repeat, 2),
                                 float_to_str(bridge.loop_count_by_middle, 2),
                                 str(bridge.loop_count), float_to_str(bridge.quality, 1)])
        print_table(bridge_table, alignments='RRRRRRRR', left_align_header=False, indent=0,
                    fixed_col_widths=[5, 6, 6, 5, 10, 10, 5, 7])
    else:
        log.log('No loop unrolling bridges made')

    return bridges
