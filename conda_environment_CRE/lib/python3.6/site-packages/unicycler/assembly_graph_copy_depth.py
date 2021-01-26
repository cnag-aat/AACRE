"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module describes the logic used to assign copy depths to assembly graph segments.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

from .misc import print_table, get_right_arrow
from . import settings
from . import log


def determine_copy_depth(graph):
    """
    Assigns a copy depth to each segment in the graph.
    """
    # Reset any existing copy depths.
    graph.copy_depths = {}

    log.log_section_header('Determining graph multiplicity')
    log.log_explanation('Multiplicity is the number of times a '
                        'sequence occurs in the underlying sequence. Single-copy contigs '
                        '(those with a multiplicity of one, occurring only once in the underlying '
                        'sequence) are particularly useful.',
                        verbosity=1, extra_empty_lines_after=0)

    log.log_explanation('Multiplicity determination begins by identifying contigs which are '
                        'clearly single-copy because they are of low depth and do not have more '
                        'than one link per side.',
                        verbosity=2)

    single_copy_depth = graph.get_single_copy_depth()

    # Assign single copy status to segments within the tolerance of the single copy depth.
    # Also, if the graph has manually set multiplicity values (using the ML tag), we also use them
    # to assign single copy status.
    max_depth = single_copy_depth + settings.INITIAL_SINGLE_COPY_TOLERANCE
    initial_single_copy_segments = []
    for segment in sorted([x for x in graph.segments.values()],
                          key=lambda x: x.get_length(), reverse=True):
        num = segment.number
        depth = segment.depth
        if (depth <= max_depth and okay_for_initial_single_copy(graph, segment)) or \
                (num in graph.manual_multiplicity and graph.manual_multiplicity[num] == 1):
            graph.copy_depths[segment.number] = [segment.depth]
            initial_single_copy_segments.append(segment.number)

    if initial_single_copy_segments:
        log.log('\nInitial single copy segments:', 2)
        log.log_number_list(initial_single_copy_segments, 2)
    else:
        log.log('Initial single copy segments: none', 2)

    log.log('', verbosity=2)
    log.log_explanation('Unicycler now uses a greedy algorithm to propagate multiplicity '
                        'through the graph. For example, if two single-copy contigs merge '
                        'together, the resulting contig will get a multiplicity of two. '
                        'When no more propagation is possible, additional single-copy contigs '
                        'are added and the process is repeated. This allows for multiplicity to be '
                        'correctly assigned to the chromosome (at the median depth) but also for '
                        'plasmids (which may be higher or lower in depth).',
                        verbosity=2)

    # Propagate copy depth as much as possible using those initial assignments.
    copy_depth_table = [['Input', '', 'Output']]
    determine_copy_depth_part_2(graph, settings.COPY_PROPAGATION_TOLERANCE, copy_depth_table)

    # Assign single copy to the largest available segment, propagate and repeat.
    while True:
        assignments = assign_single_copy_depth(graph, settings.MIN_SINGLE_COPY_LENGTH,
                                               copy_depth_table)
        determine_copy_depth_part_2(graph, settings.COPY_PROPAGATION_TOLERANCE, copy_depth_table)
        if not assignments:
            break

    # Now propagate with no tolerance threshold to complete the remaining segments.
    if log.logger.stdout_verbosity_level >= 3:
        copy_depth_table.append(['REMOVING PROPAGATION TOLERANCE', '', ''])
    determine_copy_depth_part_2(graph, 1.0, copy_depth_table)

    print_table(copy_depth_table, alignments='RLL', max_col_width=999, hide_header=True,
                indent=0, col_separation=1, verbosity=2)


def determine_copy_depth_part_2(graph, tolerance, copy_depth_table):
    """
    Propagates copy depth repeatedly until assignments stop.
    """
    if log.logger.stdout_verbosity_level >= 3:
        copy_depth_table.append(['MERGING MULTIPLICITY', '', ''])
    while merge_copy_depths(graph, tolerance, copy_depth_table):
        pass
    if log.logger.stdout_verbosity_level >= 3:
        copy_depth_table.append(['SPLITTING MULTIPLICITY', '', ''])
    if redistribute_copy_depths(graph, tolerance, copy_depth_table):
        determine_copy_depth_part_2(graph, tolerance, copy_depth_table)


def assign_single_copy_depth(graph, min_single_copy_length, copy_depth_table):
    """
    This function assigns a single copy to the longest available segment.
    """
    if log.logger.stdout_verbosity_level >= 3:
        copy_depth_table.append(['FINDING NEW SINGLE-COPY', '', ''])
    segments = sorted(get_segments_without_copies(graph), key=lambda x: x.get_length(),
                      reverse=True)
    for segment in segments:
        if segment.get_length() < min_single_copy_length:
            continue
        num = segment.number
        if num in graph.manual_multiplicity and graph.manual_multiplicity[num] != 1:
            continue
        if exactly_one_link_per_end(graph, segment):
            name_depth_before = get_seg_name_depth_str(graph, segment.number)
            graph.copy_depths[segment.number] = [segment.depth]
            name_depth_after = get_seg_name_depth_str(graph, segment.number)
            add_to_copy_depth_table(name_depth_before, name_depth_after, copy_depth_table)
            return 1
    return 0


def merge_copy_depths(graph, error_margin, copy_depth_table):
    """
    This function looks for segments where they have input on one end where:
      1) All input segments have copy depth assigned.
      2) All input segments exclusively input to this segment.
    All such cases are evaluated, and the segment with the lowest error (if that error is below
    the allowed error margin) is assigned copy depths, scaling the inputs so their sum
    exactly matches the segment's depth.
    """
    segments = get_segments_without_copies(graph)
    if not segments:
        return 0

    best_segment_num = None
    best_source_nums = None
    best_new_depths = []
    lowest_error = float('inf')

    for segment in segments:
        num = segment.number
        exclusive_inputs = graph.get_exclusive_inputs(num)
        exclusive_outputs = graph.get_exclusive_outputs(num)
        in_depth_possible = exclusive_inputs and all_have_copy_depths(graph, exclusive_inputs)
        out_depth_possible = exclusive_outputs and all_have_copy_depths(graph, exclusive_outputs)
        if in_depth_possible:
            depths, error = scale_copy_depths_from_source_segments(graph, num, exclusive_inputs)
            conflict = (num in graph.manual_multiplicity and
                        graph.manual_multiplicity[num] != len(depths))
            if error < lowest_error and not conflict:
                lowest_error = error
                best_segment_num = num
                best_source_nums = exclusive_inputs
                best_new_depths = depths
        if out_depth_possible:
            depths, error = scale_copy_depths_from_source_segments(graph, num, exclusive_outputs)
            conflict = (num in graph.manual_multiplicity and
                        graph.manual_multiplicity[num] != len(depths))
            if error < lowest_error and not conflict:
                lowest_error = error
                best_segment_num = num
                best_source_nums = exclusive_outputs
                best_new_depths = depths
    if best_segment_num and lowest_error < error_margin:
        graph.copy_depths[best_segment_num] = best_new_depths
        add_to_copy_depth_table(' + '.join(get_seg_name_depth_str(graph, x)
                                           for x in best_source_nums),
                                get_seg_name_depth_str(graph, best_segment_num), copy_depth_table)
        return 1
    else:
        return 0


def get_seg_name_depth_str(graph, segment_num):
    if segment_num in graph.copy_depths and len(graph.copy_depths[segment_num]) > 0:
        copy_str = '(' + str(len(graph.copy_depths[segment_num])) + 'x)'
    else:
        copy_str = '(none)'
    return str(segment_num) + copy_str


def add_to_copy_depth_table(before_str, after_str, copy_depth_table):
    while len(before_str) > settings.COPY_DEPTH_PROPAGATION_TABLE_ROW_WIDTH:
        before_str_parts = before_str.split(' + ')
        partial_before_str = before_str_parts[0]
        before_str_parts.pop(0)
        while before_str_parts:
            if (len(partial_before_str) + 3 + len(before_str_parts[0]) <=
                    settings.COPY_DEPTH_PROPAGATION_TABLE_ROW_WIDTH):
                partial_before_str += ' + ' + before_str_parts[0]
                before_str_parts.pop(0)
            else:
                break
        copy_depth_table.append([partial_before_str, '', ''])
        before_str = ' + ' + ' + '.join(before_str_parts)

    copy_depth_table.append([before_str, get_right_arrow(), after_str])


def redistribute_copy_depths(graph, error_margin, copy_depth_table):
    """
    This function deals with the easier case of copy depth redistribution: where one segments
    with copy depth leads exclusively to multiple segments without copy depth.
    We will then try to redistribute the source segment's copy depths among the destination
    segments.  If it can be done within the allowed error margin, the destination segments will
    get their copy depths.
    """
    segments = get_segments_with_two_or_more_copies(graph)
    if not segments:
        return 0

    for segment in segments:
        num = segment.number
        connections = graph.get_exclusive_inputs(num)
        if not connections or all_have_copy_depths(graph, connections):
            connections = graph.get_exclusive_outputs(num)
        if not connections or all_have_copy_depths(graph, connections):
            continue

        # If we got here, then we can try to redistribute the segment's copy depths to its
        # connections which are lacking copy depth.
        copy_depths = graph.copy_depths[num]
        bins = [[]] * len(connections)
        targets = [None if x not in graph.copy_depths else len(graph.copy_depths[x])
                   for x in connections]

        # For cases where there are many copy depths being distributed to many segments, there
        # will be too many combinations, so we don't bother trying.
        arrangement_count = len(bins) ** len(copy_depths)
        if arrangement_count > settings.MAX_COPY_DEPTH_DISTRIBUTION_ARRANGEMENTS:
            continue
        arrangements = shuffle_into_bins(copy_depths, bins, targets)
        if not arrangements:
            continue

        lowest_error = float('inf')
        best_arrangement = None
        for i, arrangement in enumerate(arrangements):
            error = get_error_for_multiple_segments_and_depths(graph, connections, arrangement)
            if i == 0 or error < lowest_error:
                lowest_error = error
                best_arrangement = arrangement

        # Make sure this redistribution of copy depths does not conflict with any manually assigned
        # multiplicities.
        conflict = False
        if best_arrangement is not None:
            for connection_num, connection_depths in zip(connections, best_arrangement):
                if (connection_num in graph.manual_multiplicity and
                        graph.manual_multiplicity[connection_num] != len(connection_depths)):
                    conflict = True

        if lowest_error < error_margin and not conflict:
            if assign_copy_depths_where_needed(graph, connections, best_arrangement, error_margin):
                add_to_copy_depth_table(get_seg_name_depth_str(graph, num),
                                        ' + '.join(get_seg_name_depth_str(graph, x)
                                                   for x in connections),
                                        copy_depth_table)
                return 1
    return 0


def okay_for_initial_single_copy(graph, segment):
    """
    Returns True if the given segment's links don't preclude calling this a single copy segment
    for the initial round of copy depth assignment.
    """
    num = segment.number
    forward_count, reverse_count = 0, 0
    if num in graph.forward_links:
        forward_count = len(graph.forward_links[num])
    if num in graph.reverse_links:
        reverse_count = len(graph.reverse_links[num])

    # If a segment is short, then we want to be particularly strict about assigning initial
    # single copy status. The only way for a short segment to pass is by having exactly one
    # link per side and for the segments it is linked to not be single copy themselves.
    if segment.get_length() < settings.MIN_SINGLE_COPY_LENGTH:
        if forward_count != 1 or reverse_count != 1:
            return False
        downstream_seg = abs(graph.forward_links[num][0])
        if downstream_seg in graph.copy_depths and len(graph.copy_depths[downstream_seg]) == 1:
            return False
        upstream_seg = abs(graph.reverse_links[num][0])
        if upstream_seg in graph.copy_depths and len(graph.copy_depths[upstream_seg]) == 1:
            return False
        return True

    # If the code got here, then the segment is longer and we're a bit more lenient with
    # allowing initial single copy status.

    # Having 1 link per side is ideal, but 0 links are okay too.
    forward_okay = forward_count <= 1
    reverse_okay = reverse_count <= 1

    # If either direction has too many links, we'll still accept it if the copy depths are
    # very inconsistent. This is because very inconsistent depths (e.g. 1 + 1 = 1) tend to
    # indicate bogus connections.
    if not forward_okay:
        exclusive_outputs = graph.get_exclusive_outputs(num)
        if exclusive_outputs:
            downstream_depth_sum = sum(graph.segments[x].depth for x in exclusive_outputs)
            error = get_error(downstream_depth_sum, segment.depth)
            if error > settings.COPY_PROPAGATION_TOLERANCE:
                forward_okay = True
    if not reverse_okay:
        exclusive_inputs = graph.get_exclusive_inputs(num)
        if exclusive_inputs:
            upstream_depth_sum = sum(graph.segments[x].depth for x in exclusive_inputs)
            error = get_error(upstream_depth_sum, segment.depth)
            if error > settings.COPY_PROPAGATION_TOLERANCE:
                reverse_okay = True

    # Both sides have to be okay for the segment to pass.
    return forward_okay and reverse_okay


def exactly_one_link_per_end(graph, segment):
    """
    Returns True if the given segment has exactly one link on either end.
    """
    num = segment.number
    if num in graph.forward_links and len(graph.forward_links[num]) != 1:
        return False
    if num in graph.reverse_links and len(graph.reverse_links[num]) != 1:
        return False
    return True


def all_have_copy_depths(graph, segment_numbers):
    """
    Takes a list of segment numbers and returns whether every segment in the list has copy
    depths assigned.
    """
    for num in segment_numbers:
        if num not in graph.copy_depths:
            return False
    return True


def scale_copy_depths_from_source_segments(graph, segment_number, source_segment_numbers):
    """
    Using a list of segments which are the source of copy depth, this function scales them so
    that their sum matches the depth of the given segment.
    It returns:
      1) a list of depth numbers
      2) the error (i.e. the degree of scaling which had to occur)
    It assumes that all of the source segments definitely have copy depths.
    """
    source_depths = []
    for num in source_segment_numbers:
        source_depths += graph.copy_depths[num]
    target_depth = graph.segments[segment_number].depth
    return scale_copy_depths(target_depth, source_depths)


def scale_copy_depths(target_depth, source_depths):
    """
    This function takes the source depths and scales them so their sum matches the target
    depth.  It returns the scaled depths and the error.
    """
    source_depth_sum = sum(source_depths)
    scaling_factor = target_depth / source_depth_sum
    scaled_depths = sorted([scaling_factor * x for x in source_depths], reverse=True)
    error = get_error(source_depth_sum, target_depth)
    return scaled_depths, error


def get_segments_without_copies(graph):
    """
    Returns a list of the graph segments lacking copy depth information.
    """
    return [x for x in graph.segments.values() if x.number not in graph.copy_depths]


def get_segments_with_two_or_more_copies(graph):
    return [x for x in graph.segments.values() if
            x.number in graph.copy_depths and len(graph.copy_depths[x.number]) > 1]


def get_error_for_multiple_segments_and_depths(graph, segment_numbers, copy_depths):
    """
    For the given segments, this function assesses how well the given copy depths match up.
    The maximum error for any segment is what's returned at the end.
    """
    max_error = 0.0
    for i, num in enumerate(segment_numbers):
        segment_depth = graph.segments[num].depth
        depth_sum = sum(copy_depths[i])
        max_error = max(max_error, get_error(depth_sum, segment_depth))
    return max_error


def assign_copy_depths_where_needed(graph, segment_numbers, new_depths, error_margin):
    """
    For the given segments, this function assigns the corresponding copy depths, scaled to fit
    the segment.  If a segment already has copy depths, it is skipped (i.e. this function only
    write new copy depths, doesn't overwrite existing ones).
    It will only create copy depths if doing so is within the allowed error margin.
    """
    success = False
    for i, num in enumerate(segment_numbers):
        if num not in graph.copy_depths:
            new_copy_depths, error = scale_copy_depths(graph.segments[num].depth, new_depths[i])
            if error <= error_margin:
                graph.copy_depths[num] = new_copy_depths
                success = True
    return success


def get_error(source, target):
    """
    Returns the relative error from trying to assign the source value to the target value.
    E.g. if source = 1.6 and target = 2.0, the error is 0.2
    """
    if target > 0.0:
        return abs(source - target) / target
    else:
        return float('inf')


def shuffle_into_bins(items, bins, targets):
    """
    Shuffle items into bins in all possible arrangements that satisfy these conditions:
      1) All bins must have at least one item.
      2) Any bins with a specified target must have exactly that number of items.
    """
    arrangements = []

    # If there are items not yet in a bin, place the first item in each possible bin and call this
    # function recursively.
    if items:

        # If there are only enough items to fill the empty bins, then we will only put the next
        # item in an empty bin (because putting it in a non-empty bin would prevent us from filling
        # all bins).
        empty_bin_count = sum(1 for x in bins if not x)
        only_put_in_empty = len(items) <= empty_bin_count

        for i, _ in enumerate(bins):

            # Don't put an item in a bin if that bin is already at capacity.
            if targets[i] and len(bins[i]) >= targets[i]:
                continue

            if only_put_in_empty and bins[i]:
                continue

            bins_copy = [list(x) for x in bins]
            bins_copy[i].append(items[0])
            arrangements += shuffle_into_bins(items[1:], bins_copy, targets)

    # If all items are in a bin, all bins have at least one item and any bins with a target have
    # the appropriate amount, then add the arrangement to the results.
    elif all(x for x in bins) and \
            all([not target or target == len(bins[i]) for i, target in enumerate(targets)]):
        arrangements.append(bins)
    return arrangements
