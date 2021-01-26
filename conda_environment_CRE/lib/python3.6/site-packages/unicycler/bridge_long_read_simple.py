"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module takes care of simple long read bridging. This is where parts of the graph have a very
straightforward repeat and minimap alignments can conclusively support one over the alternatives.
Simple long read bridging is done before other, more complex bridging processes.

Specifically, there are two kinds of simple bridging:
 * two way junctions
 * simple loops

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
import sys
import math
from collections import defaultdict
import itertools
from multiprocessing.dummy import Pool as ThreadPool
from .minimap_alignment import align_long_reads_to_assembly_graph, build_start_end_overlap_sets
from .misc import print_table, get_right_arrow, float_to_str
from .bridge_common import get_bridge_str, get_mean_depth, get_depth_agreement_factor
from . import log
from . import settings

try:
    from .cpp_wrappers import fully_global_alignment
except AttributeError as att_err:
    sys.exit('Error when importing C++ library: ' + str(att_err) + '\n'
             'Have you successfully built the library file using make?')


class SimpleLongReadBridge(object):
    def __init__(self, graph, start, end, path, votes_for, votes_against):

        # The numbers of the two single copy segments which are being bridged.
        self.start_segment = start
        self.end_segment = end

        # The path through the unbridged graph.
        self.graph_path = path

        # The bridge depth, a weighted mean of the start and end depths.
        self.depth = get_mean_depth(graph.segments[abs(self.start_segment)],
                                    graph.segments[abs(self.end_segment)], graph)

        # A score used to determine the order of bridge application. Simple long read bridges are
        # reliable and start with a high quality.
        self.quality = 1.0

        # When a bridge is applied, the segments in the bridge may have their depth reduced
        # accordingly. This member stores which segments have had their depth reduced and by how
        # much due to this bridge's application. It is stored so if this bridge is later deleted,
        # we can restore the depth to the segments.
        self.segments_reduced_depth = []

        # If there are segments in between the start and end (there usually will be), then they
        # provide the bridge sequence.
        self.bridge_sequence = graph.get_path_sequence(self.graph_path)

        # The start segment and end segment should agree in depth. If they don't, that's bad.
        start_seg = graph.segments[abs(self.start_segment)]
        end_seg = graph.segments[abs(self.end_segment)]
        self.quality *= get_depth_agreement_factor(start_seg.depth, end_seg.depth)

        # More 'landslidey' votes lead to high qualities. Close-call votes lead to low qualities.
        try:
            vote_proportion = votes_for / (votes_for + votes_against)
        except ZeroDivisionError:
            vote_proportion = 0.0
        vote_quality = max(0.0, 2.0 * (vote_proportion - 0.5))
        self.quality *= vote_quality

        # Low numbers of votes get penalised.
        if votes_for == 1:
            self.quality *= 0.5
        elif votes_for == 2:
            self.quality *= 0.75

        # We finalise the quality to a range of 0 to 100. We also use the sqrt function to pull
        # the scores up a bit (otherwise they tend to hang near the bottom of the range).
        self.quality = 100.0 * math.sqrt(self.quality)

    def __repr__(self):
        return 'Simple long read bridge: ' + get_bridge_str(self) + \
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
        return 'simple long read'


def create_simple_long_read_bridges(graph, out_dir, keep, threads, read_dict, long_read_filename,
                                    scoring_scheme, anchor_segments):
    """
    Create and return simple long read bridges.
    """
    log.log_section_header('Creating simple long read bridges')
    log.log_explanation('Unicycler uses long read alignments (from minimap) to resolve simple '
                        'repeat structures in the graph. This takes care of some "low-hanging '
                        'fruit" of the graph simplification.')

    bridging_dir = os.path.join(out_dir, 'simple_bridging')
    if not os.path.exists(bridging_dir):
        os.makedirs(bridging_dir)
    minimap_alignments = align_long_reads_to_assembly_graph(graph, long_read_filename,
                                                            bridging_dir, threads)
    start_overlap_reads, end_overlap_reads = build_start_end_overlap_sets(minimap_alignments)
    bridges = simple_bridge_two_way_junctions(graph, start_overlap_reads, end_overlap_reads,
                                              minimap_alignments, anchor_segments)
    bridges += simple_bridge_loops(graph, start_overlap_reads, end_overlap_reads,
                                   minimap_alignments, read_dict, scoring_scheme, threads,
                                   anchor_segments)
    if keep < 3:
        shutil.rmtree(bridging_dir, ignore_errors=True)
    return bridges


def simple_bridge_two_way_junctions(graph, start_overlap_reads, end_overlap_reads,
                                    minimap_alignments, segments_to_bridge):
    bridges = []
    c_with_arrows = get_right_arrow() + 'C' + get_right_arrow()
    log.log_explanation('Two-way junctions are defined as cases where two graph contigs (A and B) '
                        'join together (C) and then split apart again (D and E). This usually '
                        'represents a simple 2-copy repeat, and there are two possible options for '
                        'its resolution: (A' + c_with_arrows + 'D and B' + c_with_arrows + 'E) or '
                        '(A' + c_with_arrows + 'E and B' + c_with_arrows + 'D). '
                        'Each read which spans such a junction gets to "vote" for option 1, '
                        'option 2 or neither. Unicycler creates a bridge at each junction for '
                        'the most voted for option.')

    two_way_junctions_table = [['Junction', 'Option 1', 'Option 2', 'Op. 1 votes', 'Op. 2 votes',
                                'Neither votes', 'Final op.', 'Bridge quality']]

    junctions = graph.find_simple_two_way_junctions(segments_to_bridge)

    if not junctions:
        log.log('No suitable two-way junctions present')
        log.log('')
        return []

    for junction in junctions:
        table_row = [junction]

        inputs = graph.reverse_links[junction]
        outputs = graph.forward_links[junction]

        # There are two possible ways to resolve a simple two-way junction:
        #   1) inputs[0], junction, outputs[0]
        #      inputs[1], junction, outputs[1]
        #   2) inputs[0], junction, outputs[1]
        #      inputs[1], junction, outputs[0]]

        spaced_arrow = ' ' + get_right_arrow() + ' '
        option_1_1_str = spaced_arrow.join(str(x) for x in [inputs[0], junction, outputs[0]])
        option_1_2_str = spaced_arrow.join(str(x) for x in [inputs[1], junction, outputs[1]])
        table_row.append(option_1_1_str + ', ' + option_1_2_str)
        option_2_1_str = spaced_arrow.join(str(x) for x in [inputs[0], junction, outputs[1]])
        option_2_2_str = spaced_arrow.join(str(x) for x in [inputs[1], junction, outputs[0]])
        table_row.append(option_2_1_str + ', ' + option_2_2_str)

        # Gather up all reads which could possibly span this junction.
        relevant_reads = list(end_overlap_reads[inputs[0]] | end_overlap_reads[inputs[1]] |
                              end_overlap_reads[-outputs[0]] | end_overlap_reads[-outputs[1]] |
                              start_overlap_reads[outputs[0]] | start_overlap_reads[outputs[1]] |
                              start_overlap_reads[-inputs[0]] | start_overlap_reads[-inputs[1]])

        # Each read now casts a vote in one of four ways:
        #   1) Option 1: the read supports the connection of option 1
        #   2) Option 2: the read supports the connection of option 2
        #   3) Neither option: the read supports some other strange connection
        #   4) No vote: the read doesn't support any connection
        option_1_votes = 0
        option_2_votes = 0
        neither_option_votes = 0

        # The lists in expected_next_seg each hold three segment numbers:
        #   1) The starting segment
        #   2) The segment that should follow if option 1 is correct.
        #   3) The segment that should follow if option 2 is correct.
        expected_next_seg = [[inputs[0], outputs[0], outputs[1]],
                             [inputs[1], outputs[1], outputs[0]],
                             [-outputs[0], -inputs[0], -inputs[1]],
                             [-outputs[1], -inputs[1], -inputs[0]]]

        for r in relevant_reads:

            # Turn the alignments into a list of segment numbers, excluding the junction segment.
            alignments = [int(x.ref_name) * (-1 if x.read_strand == '-' else 1)
                          for x in minimap_alignments[r] if x.ref_name != str(junction)]
            alignments = [k for k, g in itertools.groupby(alignments)]  # remove adjacent duplicates

            for p in expected_next_seg:
                start, option_1_end, option_2_end = p[0], p[1], p[2]
                try:
                    after_start = alignments[alignments.index(start) + 1]
                    if after_start == option_1_end:
                        option_1_votes += 1
                    elif after_start == option_2_end:
                        option_2_votes += 1
                    else:
                        neither_option_votes += 1
                except (ValueError, IndexError):
                    pass

        table_row += [str(option_1_votes), str(option_2_votes), str(neither_option_votes)]

        # If there aren't any votes at all, that's a shame! The long reads were probably too few
        # and/or too short to span the junction.
        if option_1_votes == 0 and option_2_votes == 0:
            table_row += ['none', 'no reads']

        # If the best option isn't sufficiently better than the second best option, we don't bridge
        # anything.
        elif option_1_votes == option_2_votes:
            table_row += ['none', 'tie vote']

        else:
            # If we got here, then we're good to bridge!
            start_1 = inputs[0]
            start_2 = inputs[1]
            if option_1_votes > option_2_votes:
                table_row.append('1')
                end_1 = outputs[0]
                end_2 = outputs[1]
                votes_for = option_1_votes
                votes_against = option_2_votes + neither_option_votes
            else:  # option 2:
                table_row.append('2')
                end_1 = outputs[1]
                end_2 = outputs[0]
                votes_for = option_2_votes
                votes_against = option_1_votes + neither_option_votes
            bridges.append(SimpleLongReadBridge(graph, start_1, end_1, [junction], votes_for,
                                                votes_against))
            bridges.append(SimpleLongReadBridge(graph, start_2, end_2, [junction], votes_for,
                                                votes_against))
            table_row.append(float_to_str(bridges[-1].quality, 1))

        two_way_junctions_table.append(table_row)

    max_op_1_len = max(max(len(y) for y in x[1].split(', ')) for x in two_way_junctions_table) + 1
    max_op_2_len = max(max(len(y) for y in x[2].split(', ')) for x in two_way_junctions_table) + 1
    print_table(two_way_junctions_table, alignments='RCCRRRRR', left_align_header=False, indent=0,
                fixed_col_widths=[8, max_op_1_len, max_op_2_len, 5, 5, 7, 5, 7],
                sub_colour={'no reads': 'red', 'tie vote': 'red'})
    log.log('')
    return bridges


def simple_bridge_loops(graph, start_overlap_reads, end_overlap_reads, minimap_alignments,
                        read_dict, scoring_scheme, threads, segments_to_bridge):
    bridges = []
    ra = get_right_arrow()
    zero_loops = 'A' + ra + 'C' + ra + 'B'
    one_loop = 'A' + ra + 'C' + ra + 'D' + ra + 'C' + ra + 'B'
    two_loops = 'A' + ra + 'C' + ra + 'D' + ra + 'C' + ra + 'D' + ra + 'C' + ra + 'B'
    log.log_explanation('Simple loops are parts of the graph where two contigs (A and B) '
                        'are connected via a repeat (C) which loops back to itself (via D). It '
                        'is possible to traverse the loop zero times (' + zero_loops + '), one '
                        'time (' + one_loop + '), two times (' + two_loops + '), etc. '
                        'Long reads which span the loop inform which is the correct number of '
                        'times through. In this step, such reads are found and each is aligned '
                        'against alternative loop counts. A reads casts its "vote" for the loop '
                        'count it agrees best with, and Unicycler creates a bridge using the '
                        'most voted for count.')

    loops = sorted(graph.find_all_simple_loops())
    seg_nums_to_bridge = set(x.number for x in segments_to_bridge)
    loops = [x for x in loops
             if abs(x[0]) in seg_nums_to_bridge and abs(x[1]) in seg_nums_to_bridge and
             abs(x[3]) not in seg_nums_to_bridge]
    loops = [x for x in loops if abs(x[0]) != abs(x[1])]
    if not loops:
        log.log('No suitable simple loops present')
        return []

    col_widths = [5, 6, 6, 5, 5, 18, 5, 7]
    loop_table_header = ['Start', 'Repeat', 'Middle', 'End', 'Read count', 'Read votes',
                         'Loop count', 'Bridge quality']
    print_table([loop_table_header], fixed_col_widths=col_widths, left_align_header=False,
                alignments='RRRRRLRR', indent=0)

    for start, end, middle, repeat in loops:
        if middle is None:
            loop_table_row = [start, repeat, '', end]
        else:
            loop_table_row = [start, repeat, middle, end]

        forward_strand_reads = end_overlap_reads[start] & start_overlap_reads[end]
        reverse_strand_reads = end_overlap_reads[-end] & start_overlap_reads[-start]

        all_reads = list(forward_strand_reads) + list(reverse_strand_reads)
        strands = ['F'] * len(forward_strand_reads) + ['R'] * len(reverse_strand_reads)
        loop_table_row.append(len(all_reads))

        # This dictionary will collect the votes. The key is the number of times through the loop
        # and the value is the vote count. Votes for -1 times through the loop occur for reads
        # which don't conform to the loop assumption.
        votes = defaultdict(int)

        # We'll try a range of repeat counts. The segment depth gives us a first guess as to the
        # repeat count, which guides how high we should test.
        mean_start_end_depth = (graph.segments[abs(start)].depth +
                                graph.segments[abs(end)].depth) / 2
        if middle is None:
            repeat_depth = graph.segments[abs(repeat)].depth
            best_repeat_guess = int(round(repeat_depth / mean_start_end_depth)) - 1
        else:
            middle_depth = graph.segments[abs(middle)].depth
            best_repeat_guess = int(round(middle_depth / mean_start_end_depth))
        best_repeat_guess = max(1, best_repeat_guess)
        max_tested_loop_count = (best_repeat_guess + 1) * 2

        # Use a simple loop if we only have one thread.
        if threads == 1:
            for read, strand in zip(all_reads, strands):
                vote = get_read_loop_vote(start, end, middle, repeat, strand, minimap_alignments,
                                          read, read_dict, graph, max_tested_loop_count,
                                          scoring_scheme)
                votes[vote] += 1

        # Use a thread pool if we have more than one thread.
        else:
            pool = ThreadPool(threads)
            arg_list = []
            for read, strand in zip(all_reads, strands):
                arg_list.append((start, end, middle, repeat, strand, minimap_alignments,
                                 read, read_dict, graph, max_tested_loop_count, scoring_scheme))
            for vote in pool.imap_unordered(get_read_loop_vote_one_arg, arg_list):
                votes[vote] += 1

        # Format the vote totals nicely for the table.
        vote_str = ''
        for loop_count in sorted(votes.keys()):
            if loop_count == -1:
                vote_str += 'bad: '
            elif loop_count == 1:
                vote_str += '1 loop: '
            else:
                vote_str += str(loop_count) + ' loops: '
            vote_count = votes[loop_count]
            vote_str += str(vote_count) + ' vote' + ('s' if vote_count != 1 else '') + '    '
        loop_table_row.append(vote_str.strip())

        # Determine the repeat count which wins!
        results = sorted(list(votes.items()), key=lambda x: x[1], reverse=True)
        if not results:
            loop_table_row += ['no reads', '']
        else:
            winning_loop_count = results[0][0]
            winning_votes = results[0][1]
            if len(results) == 1:
                second_best_votes = 0
                votes_against = 0
            else:
                second_best_votes = results[1][1]
                votes_against = sum(r[1] for r in results) - winning_votes
            if winning_loop_count == -1:
                loop_table_row += ['bad reads', '']
            elif winning_votes == second_best_votes:
                loop_table_row += ['tie vote', '']
            else:
                # If we got here, then we're good to bridge!
                loop_table_row.append(str(winning_loop_count))

                bridge_path = [repeat]
                for _ in range(winning_loop_count):
                    if middle is not None:
                        bridge_path.append(middle)
                    bridge_path.append(repeat)

                bridges.append(SimpleLongReadBridge(graph, start, end, bridge_path, winning_votes,
                                                    votes_against))
                loop_table_row.append(float_to_str(bridges[-1].quality, 1))

        print_table([loop_table_row], fixed_col_widths=col_widths, header_format='normal',
                    alignments='RRRRRLRR', left_align_header=False, bottom_align_header=False,
                    sub_colour={'bad reads': 'red', 'no reads': 'red', 'tie vote': 'red'}, indent=0)
    return bridges


def get_read_loop_vote_one_arg(all_args):
    start, end, middle, repeat, strand, minimap_alignments, read, read_dict, graph, \
        max_tested_loop_count, scoring_scheme = all_args
    return get_read_loop_vote(start, end, middle, repeat, strand, minimap_alignments, read,
                              read_dict, graph, max_tested_loop_count, scoring_scheme)


def get_read_loop_vote(start, end, middle, repeat, strand, minimap_alignments, read, read_dict,
                       graph, max_tested_loop_count, scoring_scheme):
    if strand == 'F':
        s, e, m, r = start, end, middle, repeat
    else:  # strand == 'R'
        if middle is None:
            s, e, m, r = -end, -start, None, -repeat
        else:
            s, e, m, r = -end, -start, -middle, -repeat
    alignments = minimap_alignments[read]

    last_index_of_start = -1
    for i, a in enumerate(alignments):
        if a.get_signed_ref_name() == str(s):
            last_index_of_start = i
    first_index_of_end = -1
    for i in range(last_index_of_start + 1, len(alignments)):
        a = alignments[i]
        if a.get_signed_ref_name() == str(e):
            first_index_of_end = i
            break

    # We should now have the indices of the alignments around the repeat.
    if last_index_of_start == -1 or first_index_of_end == -1:
        return -1  # vote for bad read

    # If there are any alignments in between the start and end segments, they are only
    # allowed to be the middle or repeat segments.
    bad_middle = False
    for i in range(last_index_of_start + 1, first_index_of_end):
        ref_name = alignments[i].get_signed_ref_name()
        if m is None:
            if ref_name != str(r):
                bad_middle = True
        else:
            if ref_name != str(m) and ref_name != str(r):
                bad_middle = True
    if bad_middle:
        return -1  # vote for bad read

    start_alignment = alignments[last_index_of_start]
    end_alignment = alignments[first_index_of_end]

    # Now that we have the alignments, we can extract the relevant part of the read...
    read_start_pos, read_end_pos = start_alignment.read_start, end_alignment.read_end
    read_seq = read_dict[read].sequence[read_start_pos:read_end_pos]

    # ... and the relevant parts of the start/end segments.
    if start_alignment.read_strand == '+':
        start_seg_start_pos = start_alignment.ref_start
    else:  # start_alignment.read_strand == '-'
        start_seg_start_pos = start_alignment.ref_length - start_alignment.ref_end
    if end_alignment.read_strand == '+':
        end_seg_end_pos = end_alignment.ref_end
    else:  # end_alignment.read_strand == '-'
        end_seg_end_pos = end_alignment.ref_length - end_alignment.ref_start
    start_seg_seq = graph.seq_from_signed_seg_num(s)[start_seg_start_pos:]
    end_seg_seq = graph.seq_from_signed_seg_num(e)[:end_seg_end_pos]

    if m is None:
        middle_seq = ''
    else:
        middle_seq = graph.seq_from_signed_seg_num(m)
    repeat_seq = graph.seq_from_signed_seg_num(r)

    best_score = None
    best_count = None

    loop_count = 0
    while True:
        test_seq = start_seg_seq + repeat_seq
        for _ in range(loop_count):
            test_seq += middle_seq + repeat_seq
        test_seq += end_seg_seq
        alignment_result = fully_global_alignment(read_seq, test_seq, scoring_scheme, True,
                                                  settings.SIMPLE_REPEAT_BRIDGING_BAND_SIZE)
        if alignment_result:
            seqan_parts = alignment_result.split(',', 9)
            test_seq_score = int(seqan_parts[6])
            if best_score is None or test_seq_score > best_score:
                best_score = test_seq_score
                best_count = loop_count

        # Break when we've hit the max loop count. But if the max isn't our best, then we keep
        # trying higher.
        if loop_count >= max_tested_loop_count and loop_count != best_count:
            break

        # Just in case to prevent an infinite loop.
        if loop_count > max_tested_loop_count * 10:
            break

        loop_count += 1

    # This read now casts its vote for the best repeat count!
    if best_count is not None:
        return best_count
