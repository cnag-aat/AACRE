"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module has functions for finding graph paths connecting two nodes given a consensus read
sequence.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import sys
from collections import defaultdict
from .misc import weighted_average, reverse_complement, get_num_agreement
from . import settings

try:
    from .cpp_wrappers import fully_global_alignment, path_alignment
except AttributeError as e:
    sys.exit('Error when importing C++ library: ' + str(e) + '\n'
             'Have you successfully built the library file using make?')


class TooManyPaths(Exception):
    pass


def get_best_paths_for_seq(graph, start_seg, end_seg, target_length, sequence, scoring_scheme,
                           expected_scaled_score):
    """
    Given a sequence and target length, this function finds the best paths from the start
    segment to the end segment.
    """
    assert graph.overlap == 0

    # Limit the path search to lengths near the target.
    min_length = min(int(round(target_length * settings.MIN_RELATIVE_PATH_LENGTH)),
                     target_length - settings.RELATIVE_PATH_LENGTH_BUFFER_SIZE)
    max_length = max(int(round(target_length * settings.MAX_RELATIVE_PATH_LENGTH)),
                     target_length + settings.RELATIVE_PATH_LENGTH_BUFFER_SIZE)

    # If there are few enough possible paths, we just try aligning to them all.
    try:
        paths = all_paths(graph, start_seg, end_seg, min_length, max_length)
        progressive_path_search = False

    # If there are too many paths to try exhaustively, we use a progressive approach to find
    # the best path.
    except TooManyPaths:
        progressive_path_search = True
        paths = progressive_path_find(graph, start_seg, end_seg, min_length, max_length,
                                      sequence, scoring_scheme, expected_scaled_score)

    # Sort by length discrepancy from the target so the closest length matches come first.
    paths = sorted(paths, key=lambda x: abs(target_length - graph.get_bridge_path_length(x)))

    # We now align the consensus to each of the possible paths.
    paths_and_scores = []
    for path in paths:
        path_len = graph.get_bridge_path_length(path)
        length_discrepancy = abs(path_len - target_length)

        # If there is a consensus sequence, then we actually do an alignment against the path.
        if sequence:
            path_seq = graph.get_path_sequence(path)
            alignment_result = fully_global_alignment(sequence, path_seq, scoring_scheme,
                                                      True, 1000)
            if not alignment_result:
                continue

            seqan_parts = alignment_result.split(',', 9)
            raw_score = int(seqan_parts[6])
            scaled_score = float(seqan_parts[7])

        # If there isn't a consensus sequence (i.e. the start and end overlap), then each
        # path is only scored on how well its length agrees with the target length.
        else:
            raw_score = get_num_agreement(path_len, target_length) * 100.0
            scaled_score = 100.0

        paths_and_scores.append((path, raw_score, length_discrepancy, scaled_score))

    # Sort the paths from highest to lowest quality.
    paths_and_scores = sorted(paths_and_scores, key=lambda x: (-x[1], x[2], -x[3]))

    # Don't bother keeping paths which are much worse than the best.
    if paths_and_scores:
        best_scaled_score = paths_and_scores[0][3]
        min_scaled_score = best_scaled_score * 0.95
        paths_and_scores = [x for x in paths_and_scores if x[3] >= min_scaled_score]

    return paths_and_scores, progressive_path_search


def all_paths(graph, start, end, min_length, max_length):
    """
    Returns a list of all paths which connect the starting segment to the ending segment and
    are within the length bounds. The start and end segments are not themselves included in the
    paths. Returns an empty list if no paths exist.
    Loops in the graph (especially loops of short segments which don't add much to the path
    length) can result in very large numbers of potential paths in complex areas. To somewhat
    manage this, we exclude paths which include too many copies of a segment. 'Too many copies'
    is defined as double the copy depth count or the double the depth over start/end depth.
    """
    if start not in graph.forward_links:
        return []

    start_seg = graph.segments[abs(start)]
    end_seg = graph.segments[abs(end)]
    start_end_depth = weighted_average(start_seg.depth, end_seg.depth,
                                       start_seg.get_length(), end_seg.get_length())
    working_paths = [[x] for x in graph.forward_links[start]]
    final_paths = []
    while working_paths:
        new_working_paths = []
        for working_path in working_paths:
            last_seg = working_path[-1]
            if last_seg == end:
                potential_result = working_path[:-1]
                if graph.get_path_length(potential_result) >= min_length:
                    final_paths.append(potential_result)
                    if len(final_paths) > settings.ALL_PATH_SEARCH_MAX_FINAL_PATHS:
                        raise TooManyPaths
            elif graph.get_path_length(working_path) <= max_length and \
                    last_seg in graph.forward_links:
                for next_seg in graph.forward_links[last_seg]:
                    max_allowed_count = graph.max_path_segment_count(next_seg, start_end_depth)
                    count_so_far = working_path.count(next_seg) + working_path.count(-next_seg)
                    if count_so_far < max_allowed_count:
                        new_working_paths.append(working_path + [next_seg])

        # If the number of working paths is too high, we give up.
        if len(working_paths) > settings.ALL_PATH_SEARCH_MAX_WORKING_PATHS:
            raise TooManyPaths
        working_paths = new_working_paths

    return final_paths


def progressive_path_find(graph, start, end, min_length, max_length, sequence, scoring_scheme,
                          expected_scaled_score):
    """
    This function is called when all_paths fails due to too many paths. It searches for paths by
    extended outward from both the start and end, making paths where the two searches meet. When
    the number of working paths gets too high, it is culled by performing alignments with the
    in-progress paths.
    """
    reverse_sequence = reverse_complement(sequence)

    # This set will collect all of the final paths produced. Paths are stored here as tuples
    # (because they are hashable for a set).
    final_paths = set()

    # We will work with a list of forward paths and a list of reverse paths. We set them up now
    # along with dictionaries to allow easy access to all working paths which end in a particular
    # segment.
    forward_working_paths = [[start]]
    reverse_working_paths = [[-end]]

    # Knowing the start/end depth lets us put some limits on how many times a segment can be in a
    # path, which we use to avoid going through loops forever.
    start_seg = graph.segments[abs(start)]
    end_seg = graph.segments[abs(end)]
    start_end_depth = weighted_average(start_seg.depth, end_seg.depth,
                                       start_seg.get_length(), end_seg.get_length())

    # If one of the two directions gets clogged, then only the other direction will be advanced.
    # If both directions get clogged, then the culling score fraction will be increased (brought
    # closer to 1.0) such that path culling is more aggressive.
    forward_clogged = False
    reverse_clogged = False

    while True:
        if not forward_clogged:
            shortest_reverse_path = min(graph.get_path_length(x[1:]) for x in reverse_working_paths)
            reverse_paths_dict = build_path_dictionary(reverse_working_paths)
            forward_working_paths = advance_paths(forward_working_paths, reverse_paths_dict,
                                                  shortest_reverse_path, final_paths, False,
                                                  sequence, scoring_scheme, expected_scaled_score,
                                                  graph, start_end_depth, max_length,
                                                  settings.PROGRESSIVE_PATH_SEARCH_SCORE_FRACTION)
            if not forward_working_paths:
                break
            elif len(forward_working_paths) > settings.PROGRESSIVE_PATH_SEARCH_MAX_WORKING_PATHS:
                forward_clogged = True

        if not reverse_clogged:
            shortest_forward_path = min(graph.get_path_length(x[1:]) for x in forward_working_paths)
            forward_paths_dict = build_path_dictionary(forward_working_paths)
            reverse_working_paths = advance_paths(reverse_working_paths, forward_paths_dict,
                                                  shortest_forward_path, final_paths, True,
                                                  reverse_sequence, scoring_scheme,
                                                  expected_scaled_score, graph, start_end_depth,
                                                  max_length,
                                                  settings.PROGRESSIVE_PATH_SEARCH_SCORE_FRACTION)
            if not reverse_working_paths:
                break
            elif len(reverse_working_paths) > settings.PROGRESSIVE_PATH_SEARCH_MAX_WORKING_PATHS:
                reverse_clogged = True

        # If both paths have clogged up, then we give up!
        if forward_clogged and reverse_clogged:
            return []

    # Trim the start/end segments, filter for appropriate length and return the final paths!
    final_paths = [list(x)[1:-1] for x in final_paths]
    return [x for x in final_paths if min_length <= graph.get_path_length(x) <= max_length]


def build_path_dictionary(path_list):
    """
    Constructs a dictionary where the key is the furthest segment in the path and the value is a
    list of paths which end in that segment. The paths are reversed because they'll be used by the
    paths coming from the opposite direction.
    """
    path_dict = defaultdict(list)
    for path in path_list:
        r_path = reverse_path(path)
        path_dict[r_path[0]].append(r_path)
    return path_dict


def reverse_path(path):
    """
    Reverses the order and sign of path.
    """
    return [-x for x in path[::-1]]


def advance_paths(working_paths, opposite_paths_dict, shortest_opposite_path,
                  final_paths, flip_new_final_paths, sequence, scoring_scheme,
                  expected_scaled_score, graph, start_end_depth, total_max_length,
                  cull_score_fraction):
    """
    This function takes the working paths for one direction and extends them until there are too
    many or there are no more.
    """
    # For this function, the longest we'll allow paths to get is the the max length minus how far
    # the other side has gotten.
    max_length = total_max_length - shortest_opposite_path

    while True:
        # If the working paths have run out or grown too large, then we're finished with this
        # round of advancing.
        if not 0 < len(working_paths) <= settings.PROGRESSIVE_PATH_SEARCH_MAX_WORKING_PATHS:
            break

        shortest_path_len = min(graph.get_path_length(x) for x in working_paths)

        # Extend the shortest working path(s) by adding downstream segments.
        new_working_paths = []
        for path in working_paths:
            path_len = graph.get_path_length(path)

            # If this path isn't the shortest path, we don't deal with it this time.
            if path_len > shortest_path_len:
                new_working_paths.append(path)

            # If it is the shortest path and has downstream segments...
            elif path[-1] in graph.forward_links:
                downstream_segments = graph.forward_links[path[-1]]
                for next_seg in downstream_segments:

                    # Make sure we haven't already used this segment too many times in the path.
                    max_allowed_count = graph.max_path_segment_count(next_seg, start_end_depth)
                    count_so_far = path.count(next_seg) + path.count(-next_seg)
                    if count_so_far < max_allowed_count:

                        # If the next segment is in the dictionary of the opposite direction's
                        # paths, that means we've found a path through to the other side!
                        if next_seg in opposite_paths_dict:
                            for final_part in opposite_paths_dict[next_seg]:
                                final_path = path + final_part
                                if flip_new_final_paths:
                                    final_path = reverse_path(final_path)
                                final_paths.add(tuple(final_path))

                        # Finally, extend the path if doing so won't make it too long.
                        if graph.get_path_length(path[1:] + [next_seg]) <= max_length:
                            new_working_paths.append(path + [next_seg])

        working_paths = new_working_paths

    # If we've exceeded the allowable working count, cull the paths down to size now.
    if len(working_paths) > settings.PROGRESSIVE_PATH_SEARCH_MAX_WORKING_PATHS:
        working_paths = cull_paths(graph, working_paths, sequence, scoring_scheme,
                                   expected_scaled_score, cull_score_fraction)

    return working_paths


def cull_paths(graph, paths, sequence, scoring_scheme, expected_scaled_score, cull_score_fraction):
    """
    Returns a reduced list of paths - the ones which best align to the given sequence.
    """
    # It's possible that all of the working paths share quite a bit in common at their
    # start. We can therefore find the common starting sequence and align to that once,
    # and then only do separate alignments for the remainder of the paths, saving some time.
    common_start = []
    smallest_seg_count = min(len(x) for x in paths)
    for i in range(smallest_seg_count):
        potential_common_seg = paths[0][i]
        for path in paths:
            if path[i] != potential_common_seg:
                break
        else:
            common_start.append(potential_common_seg)
            continue
        break

    # Align the consensus sequence to the common start of the paths. We exclude the first segment
    # (which is the start segment and not part of the consensus) and back up a little bit so our
    # different alignments to follow won't start right at a difference. I.e. it's better to begin
    # with a bit of common sequence.
    common_path_seq = graph.get_path_sequence(common_start[1:])[:-100]
    path_align_start = len(common_path_seq)
    if common_path_seq:
        alignment_result = path_alignment(common_path_seq, sequence, scoring_scheme, True, 1000)
        seq_align_start = int(alignment_result.split(',', 6)[5])
    else:
        seq_align_start = 0

    scored_paths = []
    shortest_len = min(graph.get_path_length(x[1:]) for x in paths)
    seq_after_common_path = sequence[seq_align_start:]
    for path in paths:
        path_seq_after_common_path = \
            graph.get_path_sequence(path[1:])[path_align_start:shortest_len]
        alignment_result = path_alignment(path_seq_after_common_path, seq_after_common_path,
                                          scoring_scheme, True, 500)
        if alignment_result:
            scaled_score = float(alignment_result.split(',', 8)[7])
            scored_paths.append((path, scaled_score))

    scored_paths = sorted(scored_paths, key=lambda x: x[1], reverse=True)
    if not scored_paths:
        return []

    # If our path finding has taken a wrong turn (i.e. all of the paths are wrong), then we want
    # to give up to save time. To check for this we see if the best one falls well below our
    # expectation and isn't that much better than the worst one.
    best_score = scored_paths[0][1]
    worst_score = scored_paths[-1][1]
    if best_score < 0.9 * expected_scaled_score and best_score * 0.95 < worst_score:
        return []

    # Now that each path is scored we keep the ones that are closest in score to the best one.
    surviving_paths = list(x for x in scored_paths
                           if x[1] >= best_score * cull_score_fraction)

    # If any of the surviving paths end in the same segment but have different scores, only keep
    # the ones with the top score. This is because the segments with a lower score will have the
    # same future paths, and so will always be lower. Doing this also helps to ensure that future
    # passes through this function will have a larger common start and will therefore go faster.
    surviving_paths_by_terminal_seg = {}
    for path in surviving_paths:
        terminal_seg = path[0][-1]
        score = path[1]
        if terminal_seg not in surviving_paths_by_terminal_seg:  # First
            surviving_paths_by_terminal_seg[terminal_seg] = [path]
        else:  # We've seen this terminal segment already
            current_best_score = surviving_paths_by_terminal_seg[terminal_seg][0][1]
            if score > current_best_score:
                surviving_paths_by_terminal_seg[terminal_seg] = [path]  # Replace the list
            elif score == current_best_score:
                surviving_paths_by_terminal_seg[terminal_seg].append(path)  # Add to the list
            else:  # score < current_best_score
                pass
    surviving_paths = []
    for paths_with_same_terminal_seg in surviving_paths_by_terminal_seg.values():
        surviving_paths += [x[0] for x in paths_with_same_terminal_seg]

    return surviving_paths
