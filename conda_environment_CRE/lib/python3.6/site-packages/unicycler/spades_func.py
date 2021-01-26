"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module contains functions relating to SPAdes assembly.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import os
import subprocess
import gzip
import shutil
import statistics
from .misc import round_to_nearest_odd, get_compression_type, int_to_str, quit_with_error,\
    strip_read_extensions, bold, dim, print_table, get_left_arrow, float_to_str
from .assembly_graph import AssemblyGraph
from . import log


class BadFastq(Exception):
    pass


def get_best_spades_graph(short1, short2, short_unpaired, out_dir, read_depth_filter, verbosity,
                          spades_path, threads, keep, kmer_count, min_k_frac, max_k_frac, kmers,
                          no_spades_correct, expected_linear_seqs, spades_tmp_dir):
    """
    This function tries a SPAdes assembly at different k-mers and returns the best.
    'The best' is defined as the smallest dead-end count after low-depth filtering.  If multiple
    graphs have the same dead-end count (e.g. zero!) then the highest kmer is used.
    """
    spades_dir = os.path.join(out_dir, 'spades_assembly')
    if not os.path.exists(spades_dir):
        os.makedirs(spades_dir)

    # SPAdes can possibly crash if given too many threads, so limit it to 32.
    threads = min(threads, 32)

    # Make sure that the FASTQ files look good.
    using_paired_reads = bool(short1) and bool(short2)
    using_unpaired_reads = bool(short_unpaired)
    if using_paired_reads:
        count_1, count_2 = 0, 0
        try:
            count_1 = get_read_count(short1)
        except BadFastq:
            quit_with_error('this read file is not a properly formatted FASTQ: ' + short1)
        try:
            count_2 = get_read_count(short2)
        except BadFastq:
            quit_with_error('this read file is not a properly formatted FASTQ: ' + short2)
        if count_1 != count_2:
            quit_with_error('the paired read input files have an unequal number of reads')
    if using_unpaired_reads:
        try:
            get_read_count(short_unpaired)
        except BadFastq:
            quit_with_error('this read file is not properly formatted as FASTQ: ' + short_unpaired)

    if no_spades_correct:
        reads = (short1, short2, short_unpaired)
    else:
        reads = spades_read_correction(short1, short2, short_unpaired, spades_dir, threads,
                                       spades_path, keep, spades_tmp_dir)
    if kmers is not None:
        kmer_range = kmers
    else:
        kmer_range = get_kmer_range(short1, short2, short_unpaired, spades_dir, kmer_count,
                                    min_k_frac, max_k_frac, spades_path)
    assem_dir = os.path.join(spades_dir, 'assembly')

    log.log_section_header('SPAdes assemblies')
    log.log_explanation('Unicycler now uses SPAdes to assemble the short reads. It scores the '
                        'assembly graph for each k-mer using the number of contigs (fewer is '
                        'better) and the number of dead ends (fewer is better). The score '
                        'function is 1/(c*(d+2)), where c is the contig count and d is the '
                        'dead end count.')

    # Conduct a SPAdes assembly for each k-mer and score them to choose the best.
    if verbosity > 1:
        spades_results_table = [['K-mer', 'Contigs', 'Links', 'Total length', 'N50',
                                 'Longest contig', 'Dead ends', 'Score']]
    else:
        spades_results_table = [['K-mer', 'Contigs', 'Dead ends', 'Score']]
    best_score = 0.0
    best_kmer = 0
    best_graph_filename = ''

    graph_files, insert_size_mean, insert_size_deviation = \
        spades_assembly(reads, assem_dir, kmer_range, threads, spades_path, spades_tmp_dir)

    existing_graph_files = [x for x in graph_files if x is not None]
    if not existing_graph_files:
        quit_with_error('SPAdes failed to produce assemblies. '
                        'See spades_assembly/assembly/spades.log for more info.')
    median_segment_count = statistics.median(count_segments_in_spades_fastg(x)
                                             for x in existing_graph_files)

    for graph_file, kmer in zip(graph_files, kmer_range):
        table_line = [int_to_str(kmer)]

        if graph_file is None:
            table_line += [''] * (7 if verbosity > 1 else 2)
            table_line.append('failed')
            spades_results_table.append(table_line)
            continue

        assembly_graph = AssemblyGraph(graph_file, kmer, paths_file=None,
                                       insert_size_mean=insert_size_mean,
                                       insert_size_deviation=insert_size_deviation)

        # If this graph has way too many segments, then we will just skip it because very complex
        # graphs take forever to clean up.
        # TO DO: I can remove this awkward hack if I make the graph cleaning more efficient.
        if len(assembly_graph.segments) > 4 * median_segment_count:
            table_line += [''] * (6 if verbosity > 1 else 2)
            table_line.append('too complex')
            spades_results_table.append(table_line)
            continue

        log.log('\nCleaning k{} graph'.format(kmer), 2)
        assembly_graph.clean(read_depth_filter)
        clean_graph_filename = os.path.join(spades_dir, ('k%03d' % kmer) + '_assembly_graph.gfa')
        assembly_graph.save_to_gfa(clean_graph_filename, verbosity=2)

        segment_count = len(assembly_graph.segments)
        dead_ends = assembly_graph.total_dead_end_count()

        # If the user is expecting some linear sequences, then the dead end count can be adjusted
        # down so expected dead ends don't penalise this k-mer.
        adjusted_dead_ends = max(0, dead_ends - (2 * expected_linear_seqs))
        if segment_count == 0:
            score = 0.0
        else:
            score = 1.0 / (segment_count * (adjusted_dead_ends + 2))

        # Prepare the table line for this k-mer graph.
        table_line += [int_to_str(segment_count)]
        if verbosity > 1:
            n50, shortest, lower_quartile, median, upper_quartile, longest = \
                assembly_graph.get_contig_stats()
            table_line += [int_to_str(assembly_graph.get_total_link_count()),
                           int_to_str(assembly_graph.get_total_length()),
                           int_to_str(n50), int_to_str(longest)]
        table_line += [int_to_str(dead_ends), '{:.2e}'.format(score)]
        spades_results_table.append(table_line)

        if score > best_score:
            best_kmer = kmer
            best_score = score
            best_graph_filename = graph_file

    log.log('', 2)

    if not best_kmer:
        quit_with_error('none of the SPAdes graphs were suitable for scaffolding in Unicycler')

    # If the best k-mer is the top k-mer, then SPAdes has already done the repeat resolution and
    # we can just grab it now. Easy! If the best k-mer was a different k-mer size, then we need
    # to run SPAdes again to get that repeat resolution.
    if best_kmer != kmer_range[-1]:
        new_kmer_range = [x for x in kmer_range if x <= best_kmer]
        graph_file, insert_size_mean, insert_size_deviation = \
            spades_assembly(reads, assem_dir, new_kmer_range, threads, spades_path, spades_tmp_dir,
                            just_last=True)
        best_graph_filename = graph_file
    paths_file = os.path.join(assem_dir, 'contigs.paths')
    if os.path.isfile(paths_file):
        copied_paths_file = os.path.join(spades_dir,
                                         'k' + ('%03d' % best_kmer) + '_contigs.paths')
        shutil.copyfile(paths_file, copied_paths_file)
    else:
        paths_file = None

    # Now we can load and clean the graph again, this time giving it the SPAdes contig paths.
    assembly_graph = AssemblyGraph(best_graph_filename, best_kmer, paths_file=paths_file,
                                   insert_size_mean=insert_size_mean,
                                   insert_size_deviation=insert_size_deviation)
    assembly_graph.clean(read_depth_filter)
    clean_graph_filename = os.path.join(spades_dir, 'k' + str(best_kmer) + '_assembly_graph.gfa')
    assembly_graph.save_to_gfa(clean_graph_filename, verbosity=2)

    if best_score == 0.0:
        quit_with_error('none of the SPAdes assemblies produced assembled sequence')

    # Print the SPAdes result table, highlighting the best k-mer in green.
    log.log_section_header('SPAdes assembly graph summary', 2)
    best_kmer_row = [x[0] for x in spades_results_table].index(int_to_str(best_kmer))
    print_table(spades_results_table, alignments='RRRRRRRR', indent=0,
                row_colour={best_kmer_row: 'green'},
                row_extra_text={best_kmer_row: ' ' + get_left_arrow() + 'best'})

    # Clean up.
    if keep < 3 and os.path.isdir(spades_dir):
        log.log('\nDeleting ' + spades_dir + '/')
        shutil.rmtree(spades_dir, ignore_errors=True)
    if keep < 3 and spades_tmp_dir is not None and os.path.isdir(spades_tmp_dir):
        log.log('Deleting ' + spades_tmp_dir + '/')
        shutil.rmtree(spades_tmp_dir, ignore_errors=True)

    return assembly_graph


def spades_read_correction(short1, short2, unpaired, spades_dir, threads, spades_path, keep,
                           spades_tmp_dir):
    """
    This runs SPAdes with the --only-error-correction option.
    """
    log.log_section_header('SPAdes read error correction')
    log.log_explanation('Unicycler uses the SPAdes read error correction module to reduce the '
                        'number of errors in the short read before SPAdes assembly. This can make '
                        'the assembly faster and simplify the assembly graph structure.')

    using_paired_reads = bool(short1) and bool(short2)
    using_unpaired_reads = bool(unpaired)

    # If the corrected reads already exist, then we just use them and proceed.
    corrected_1 = os.path.join(spades_dir, 'corrected_1.fastq.gz')
    corrected_2 = os.path.join(spades_dir, 'corrected_2.fastq.gz')
    corrected_u = os.path.join(spades_dir, 'corrected_u.fastq.gz')
    corrected_1_exists = os.path.isfile(corrected_1)
    corrected_2_exists = os.path.isfile(corrected_2)
    corrected_u_exists = os.path.isfile(corrected_u)

    reads_already_exist = True
    if using_paired_reads and (not corrected_1_exists or not corrected_2_exists):
        reads_already_exist = False
    if using_unpaired_reads and not corrected_u_exists:
        reads_already_exist = False

    if reads_already_exist:
        log.log('Corrected reads already exist. Will use these reads instead of running SPAdes '
                'error correction:')
        if using_paired_reads:
            log.log('  ' + corrected_1)
            log.log('  ' + corrected_2)
        if using_unpaired_reads:
            log.log('  ' + corrected_u)
        if not using_paired_reads:
            corrected_1, corrected_2 = None, None
        if not using_unpaired_reads:
            corrected_u = None
        return corrected_1, corrected_2, corrected_u

    # If the corrected reads don't exist, then we run SPAdes in error correction only mode.
    read_correction_dir = os.path.join(spades_dir, 'read_correction')
    command = [spades_path]
    if using_paired_reads:
        assert bool(short1) and bool(short2)
        command += ['-1', short1, '-2', short2]
    if using_unpaired_reads:
        command += ['-s', unpaired]
    command += ['-o', read_correction_dir, '--threads', str(threads), '--only-error-correction']
    if spades_tmp_dir is not None:
        command += ['--tmp-dir', spades_tmp_dir]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    while process.poll() is None:
        spades_output = process.stdout.readline().rstrip().decode()
        if 'Command line:' in spades_output:
            spades_output = ' '.join(spades_output.split()).replace('Command line: ', '')
            log.log('Command: ' + bold(spades_output))
            log.log('', 2)
        elif spades_output:
            log.log(dim(spades_output), 2)

    spades_error = process.stderr.readline().strip().decode()
    return_code = process.returncode
    if spades_error:
        quit_with_error('SPAdes encountered an error: ' + spades_error)
    if return_code != 0:
        quit_with_error('SPAdes crashed! Please view spades.log for more information.')

    # Read error correction should be done now, so copy the corrected read files to a more
    # permanent location.
    if using_paired_reads:
        short1_no_extension = strip_read_extensions(short1)
        short2_no_extension = strip_read_extensions(short2)

        # Try to deal with particularly strange names for paired reads.
        if short1_no_extension in short2_no_extension or short2_no_extension in short1_no_extension:
            for i in range(min(len(short1), len(short2))):
                short1_no_extension = os.path.basename(short1)[:i]
                short2_no_extension = os.path.basename(short2)[:i]
                if short1_no_extension != short2_no_extension:
                    break
        if short1_no_extension in short2_no_extension or short2_no_extension in short1_no_extension:
            quit_with_error('could not process read file names')
    else:
        short1_no_extension, short2_no_extension = '', ''

    corrected_dir = os.path.join(read_correction_dir, 'corrected')
    files = os.listdir(corrected_dir)
    for spades_file in files:
        file_path = os.path.join(corrected_dir, spades_file)
        if using_paired_reads and short1_no_extension in spades_file:
            shutil.move(file_path, corrected_1)
        elif using_paired_reads and short2_no_extension in spades_file:
            shutil.move(file_path, corrected_2)
        elif using_unpaired_reads and '_unpaired' in spades_file:
            shutil.move(file_path, corrected_u)
        elif using_unpaired_reads and not using_paired_reads and \
                spades_file.endswith('.cor.fastq.gz'):
            shutil.move(file_path, corrected_u)

    corrected_1_exists = os.path.isfile(corrected_1)
    corrected_2_exists = os.path.isfile(corrected_2)
    corrected_u_exists = os.path.isfile(corrected_u)

    if using_paired_reads and (not corrected_1_exists or not corrected_2_exists):
        quit_with_error('SPAdes read error correction failed. '
                        'See spades_assembly/read_correction/spades.log for more info.')
    if using_unpaired_reads and not corrected_u_exists:
        quit_with_error('SPAdes read error correction failed. '
                        'See spades_assembly/read_correction/spades.log for more info.')
    if keep < 3:
        shutil.rmtree(read_correction_dir, ignore_errors=True)

    if not using_paired_reads:
        corrected_1 = None
        corrected_2 = None
    if not using_unpaired_reads:
        corrected_u = None

    log.log('')
    log.log('Corrected reads:')
    if using_paired_reads:
        log.log('  ' + corrected_1)
        log.log('  ' + corrected_2)
    if using_unpaired_reads:
        log.log('  ' + corrected_u)

    return corrected_1, corrected_2, corrected_u


def spades_assembly(read_files, out_dir, kmers, threads, spades_path, spades_tmp_dir,
                    just_last=False):
    """
    This runs a SPAdes assembly, possibly continuing from a previous assembly.
    """
    short1 = read_files[0]
    short2 = read_files[1]
    unpaired = read_files[2]

    using_paired_reads = short1 is not None and short2 is not None and \
        os.path.isfile(short1) and os.path.isfile(short2)
    using_unpaired_reads = unpaired is not None and os.path.isfile(unpaired)

    kmer_string = ','.join([str(x) for x in kmers])
    command = [spades_path, '-o', out_dir, '-k', kmer_string, '--threads', str(threads)]
    if just_last:
        command += ['--restart-from', 'k' + str(kmers[-1])]
    else:
        if using_paired_reads:
            command += ['--only-assembler', '-1', short1, '-2', short2]
        if using_unpaired_reads:
            command += ['--only-assembler', '-s', unpaired]
    if spades_tmp_dir is not None:
        command += ['--tmp-dir', spades_tmp_dir]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    insert_size_mean = None
    insert_size_deviation = None
    while process.poll() is None:
        spades_output = process.stdout.readline().rstrip().decode()
        if spades_output:
            # Some SPAdes output lines use tabs where spaces would look better. Fix those up here
            # for aesthetics.
            if spades_output.startswith('Command line:') or \
                    spades_output.startswith('Restored from Command line:'):
                spades_output = ' '.join(spades_output.split())

            if spades_output.startswith('Command line:'):
                spades_output = spades_output.replace('Command line: ', '')
                log.log('Command: ' + bold(spades_output), 2)
                log.log('', 2)
            elif 'Running assembler: K' in spades_output:
                log.log(spades_output, 2)
            elif spades_output:
                log.log(dim(spades_output), 2)

        try:
            insert_size_mean = float(spades_output.split('Insert size = ')[-1].split(',')[0])
            insert_size_deviation = float(spades_output.split('deviation = ')[-1].split(',')[0])
        except ValueError:
            pass

    # If we couldn't get the insert size from the SPAdes output (e.g. it was an unpaired-reads-only
    # assembly), we'll use the read length instead.
    if insert_size_mean is None or insert_size_deviation is None:
        read_lengths = get_read_lengths(short1) + get_read_lengths(short2) + \
                       get_read_lengths(unpaired)
        insert_size_mean = statistics.mean(read_lengths)
        insert_size_deviation = max(statistics.stdev(read_lengths), 1.0)

    log.log('', 2)
    log.log('Insert size mean: ' + float_to_str(insert_size_mean, 1) + ' bp', 2)
    log.log('Insert size stdev: ' + float_to_str(insert_size_deviation, 1) + ' bp', 2)
    log.log('', 2)

    spades_error = process.stderr.readline().strip().decode()
    if spades_error:
        quit_with_error('SPAdes encountered an error: ' + spades_error)

    if just_last:
        graph_file = os.path.join(out_dir, 'K' + str(kmers[-1]), 'assembly_graph.fastg')
        return graph_file, insert_size_mean, insert_size_deviation
    else:
        graph_files = []
        for kmer in kmers:
            graph_file = os.path.join(out_dir, 'K' + str(kmer), 'assembly_graph.fastg')
            if os.path.isfile(graph_file):
                parent_dir = os.path.dirname(out_dir)
                copied_graph_file = os.path.join(parent_dir,
                                                 ('k%03d' % kmer) + '_assembly_graph.fastg')
                shutil.copyfile(graph_file, copied_graph_file)
                graph_files.append(copied_graph_file)
            else:
                graph_files.append(None)
        return graph_files, insert_size_mean, insert_size_deviation


def get_max_spades_kmer(spades_path):
    """
    SPAdes usually has a maximum k-mer size of 127, but this can be changed when compiling SPAdes,
    so this function checks the help text to see what it is.
    https://github.com/ablab/spades/issues/40
    """
    try:
        process = subprocess.Popen([spades_path, '--help'], stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        out, err = process.communicate()
        all_output = out.decode() + err.decode()
        all_output = all_output.replace('\n', ' ')
        all_output = ' '.join(all_output.split())
        max_kmer = all_output.split('must be odd and less than ')[1].split(')')[0]
        return int(max_kmer) - 1
    except (IndexError, ValueError):
        return 127


def get_kmer_range(reads_1_filename, reads_2_filename, unpaired_reads_filename, spades_dir,
                   kmer_count, min_kmer_frac, max_kmer_frac, spades_path):
    """
    Uses the read lengths to determine the k-mer range to be used in the SPAdes assembly.
    """
    log.log_section_header('Choosing k-mer range for assembly')
    log.log_explanation('Unicycler chooses a k-mer range for SPAdes based on the length of the '
                        'input reads. It uses a wide range of many k-mer sizes to maximise the '
                        'chance of finding an ideal assembly.')

    # If the k-mer range file already exists, we use its values and proceed.
    kmer_range_filename = os.path.join(spades_dir, 'kmer_range')
    if os.path.isfile(kmer_range_filename):
        with open(kmer_range_filename, 'rt') as kmer_range_file:
            kmer_range = kmer_range_file.readline().strip().split(', ')
        if kmer_range:
            try:
                kmer_range = [int(x) for x in kmer_range]
                log.log('K-mer range already exists:')
                log.log('  ' + kmer_range_filename)
                log.log('\nWill use this existing range:')
                log.log('  ' + ', '.join([str(x) for x in kmer_range]))
                return kmer_range
            except ValueError:
                pass

    max_spades_kmer = get_max_spades_kmer(spades_path)
    log.log('SPAdes maximum k-mer: {}'.format(max_spades_kmer))
    if max_spades_kmer != 127:
        log.log('    (unusual value, probably indicates custom SPAdes compilation)')

    # If the code got here, then the k-mer range doesn't already exist and we'll create one by
    # examining the read lengths.
    read_lengths = get_read_lengths(reads_1_filename) + get_read_lengths(reads_2_filename) + \
        get_read_lengths(unpaired_reads_filename)
    read_lengths = sorted(read_lengths)
    median_read_length = read_lengths[len(read_lengths) // 2 - 1]
    max_kmer = round_to_nearest_odd(max_kmer_frac * median_read_length)
    if max_kmer > max_spades_kmer:
        max_kmer = max_spades_kmer
    starting_kmer = round_to_nearest_odd(min_kmer_frac * max_kmer / max_kmer_frac)
    if starting_kmer < 11:
        starting_kmer = 11

    if kmer_count == 1:
        kmer_range = [max_kmer]
    elif kmer_count == 2:
        kmer_range = [starting_kmer, max_kmer]
    else:
        # Create the k-mer range from a non-linear function that spaces out the early k-mers more
        # and makes the later k-mers (which are most likely to be the good, used ones) closer
        # together.
        kmer_range = []
        for x in [x / (kmer_count - 1) for x in range(kmer_count)]:
            kmer_range.append((max_kmer - starting_kmer) * (2 - 2 / (x + 1)) + starting_kmer)
        kmer_range = sorted(list(set([round_to_nearest_odd(x) for x in kmer_range])))

    kmer_range_str = ', '.join([str(x) for x in kmer_range])

    log.log('Median read length: ' + str(median_read_length))
    log.log('K-mer range: ' + kmer_range_str)

    kmer_range_file = open(kmer_range_filename, 'w')
    kmer_range_file.write(kmer_range_str)
    kmer_range_file.close()
    return kmer_range


def get_read_lengths(reads_filename):
    """
    Returns a list of the read lengths for the given read file.
    """
    if reads_filename is None:
        return []
    if get_compression_type(reads_filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    with open_func(reads_filename, 'rb') as reads:
        read_lengths = []
        i = 0
        for line in reads:
            if i % 4 == 1:
                read_lengths.append(len(line.strip()))
            i += 1
    return read_lengths


def get_read_count(reads_filename):
    """
    Returns the number of reads in the given file.
    """
    if reads_filename is None:
        return 0
    if get_compression_type(reads_filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    with open_func(reads_filename, 'rb') as reads:
        read_count = 0
        i = 0
        for line in reads:
            if i % 4 == 0:
                try:
                    assert line.startswith(b'@')
                except AssertionError:
                    raise BadFastq
                read_count += 1
            i += 1
    return read_count


def count_segments_in_spades_fastg(fastg_file):
    seq_count = 0
    with open(fastg_file, 'rt') as fastg:
        for line in fastg:
            if line.startswith('>'):
                seq_count += 1
    return seq_count // 2
