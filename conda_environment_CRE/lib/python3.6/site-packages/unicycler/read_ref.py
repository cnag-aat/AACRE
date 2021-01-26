"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module contains classes for reads and reference sequences, and related functions.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import random
import gzip
import os
import math
from .misc import quit_with_error, get_nice_header, get_compression_type, get_sequence_file_type,\
    strip_read_extensions, print_table, float_to_str, range_is_contained, range_overlap_size, \
    simplify_ranges, add_line_breaks_to_sequence
from . import settings
from . import log


def load_references(fasta_filename, contamination=False, section_header='Loading references',
                    show_progress=True):
    """
    This function loads in sequences from a FASTA file and returns a list of Reference objects.
    """
    references = []
    total_bases = 0
    if section_header:
        log.log_section_header(section_header)
    try:
        if get_sequence_file_type(fasta_filename) != 'FASTA':
            quit_with_error(fasta_filename + ' is not in FASTA format')
    except ValueError:
        quit_with_error(fasta_filename + ' is not in FASTA format')

    if get_compression_type(fasta_filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open

    with open_func(fasta_filename, 'rt') as fasta:
        num_refs = sum(1 for line in fasta if line.startswith('>'))
    if not num_refs:
        quit_with_error('There are no references sequences in ' + fasta_filename)
    if show_progress:
        log.log_progress_line(0, num_refs)

    fasta_file = open_func(fasta_filename, 'rt')
    name = ''
    sequence = ''
    last_progress = 0.0
    step = settings.LOADING_REFERENCES_PROGRESS_STEP
    for line in fasta_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):  # Header line = start of new contig
            if name:
                if contamination:
                    name = 'CONTAMINATION_' + name
                references.append(Reference(name, sequence))
                total_bases += len(sequence)
                progress = 100.0 * len(references) / num_refs
                progress_rounded_down = math.floor(progress / step) * step
                if progress == 100.0 or progress_rounded_down > last_progress:
                    if show_progress:
                        log.log_progress_line(len(references), num_refs, total_bases)
                    last_progress = progress_rounded_down
                sequence = ''
            name = get_nice_header(line[1:])
        else:
            sequence += line
    fasta_file.close()
    if name:
        if contamination:
            name = 'CONTAMINATION_' + name
        references.append(Reference(name, sequence))
        total_bases += len(sequence)
        if show_progress:
            log.log_progress_line(len(references), num_refs, total_bases)
    if show_progress:
        log.log_progress_line(len(references), len(references), total_bases, end_newline=True)

    return references


def load_long_reads(filename, silent=False, section_header='Loading reads', output_dir=None):
    """
    This function loads in long reads from a FASTQ file and returns a dictionary where key = read
    name and value = Read object. It also returns a list of read names, in the order they are in
    the file.
    """
    # Read files can be either FASTA or FASTQ and optionally gzipped.
    try:
        file_type = get_sequence_file_type(filename)
    except ValueError:
        file_type = ''
        quit_with_error(filename + ' is not in either FASTA or FASTQ format')
    if get_compression_type(filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open

    if not silent:
        log.log_section_header(section_header)

    read_dict = {}
    read_names = []
    total_bases = 0
    last_progress = 0.0
    step = settings.LOADING_READS_PROGRESS_STEP
    duplicate_read_names_found = False

    if file_type == 'FASTQ':
        with open_func(filename, 'rt') as fastq:
            num_reads = sum(1 for _ in fastq) // 4
    else:  # file_type == 'FASTA'
        with open_func(filename, 'rt') as fasta:
            num_reads = sum(1 for line in fasta if line.startswith('>'))
    if not num_reads:
        quit_with_error('There are no read sequences in ' + filename)
    if not silent:
        log.log_progress_line(0, num_reads)

    if file_type == 'FASTQ':
        with open_func(filename, 'rt') as fastq:
            for line in fastq:
                stripped_line = line.strip()
                if len(stripped_line) == 0:
                    continue
                if not stripped_line.startswith('@'):
                    continue
                original_name = stripped_line[1:].split()[0]
                sequence = next(fastq).strip()
                _ = next(fastq)
                qualities = next(fastq).strip()

                # Don't allow duplicate read names, so add a trailing number when they occur.
                name = original_name
                duplicate_name_number = 1
                while name in read_dict:
                    duplicate_read_names_found = True
                    duplicate_name_number += 1
                    name = original_name + '_' + str(duplicate_name_number)

                read_dict[name] = Read(name, sequence, qualities)
                read_names.append(name)
                total_bases += len(sequence)
                progress = 100.0 * len(read_dict) / num_reads
                progress_rounded_down = math.floor(progress / step) * step
                if progress == 100.0 or progress_rounded_down > last_progress:
                    if not silent:
                        log.log_progress_line(len(read_dict), num_reads, total_bases)
                    last_progress = progress_rounded_down

    else:  # file_type == 'FASTA'
        with open_func(filename, 'rt') as fasta:
            name = ''
            sequence = ''
            last_progress = 0.0
            for line in fasta:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):  # Header line = start of new contig
                    if name:
                        read_dict[name] = Read(name, sequence, None)
                        read_names.append(name)
                        total_bases += len(sequence)
                        progress = 100.0 * len(read_dict) / num_reads
                        progress_rounded_down = math.floor(progress / step) * step
                        if progress == 100.0 or progress_rounded_down > last_progress:
                            if not silent:
                                log.log_progress_line(len(read_dict), num_reads, total_bases)
                            last_progress = progress_rounded_down
                        sequence = ''

                    # Don't allow duplicate read names, so add a trailing number when they occur.
                    original_name = get_nice_header(line[1:])
                    name = original_name
                    duplicate_name_number = 1
                    while name in read_dict:
                        duplicate_read_names_found = True
                        duplicate_name_number += 1
                        name = original_name + '_' + str(duplicate_name_number)

                else:
                    sequence += line
            if name:
                read_dict[name] = Read(name, sequence, None)
                read_names.append(name)
                total_bases += len(sequence)
                if not silent:
                    log.log_progress_line(len(read_dict), num_reads, total_bases)

    if not silent:
        log.log_progress_line(len(read_dict), len(read_dict), total_bases, end_newline=True)

    # If there were duplicate read names, then we save the reads back out to file with their fixed
    # names.
    if duplicate_read_names_found:
        no_dup_filename = strip_read_extensions(filename) + '_no_duplicates'
        if file_type == 'FASTQ':
            no_dup_filename += '.fastq.gz'
        else:  # file_type == 'FASTA'
            no_dup_filename += '.fasta.gz'

        # If an output directory was provided, we put the no duplicate read file there.
        if output_dir is not None:
            no_dup_filename = os.path.join(output_dir, no_dup_filename)

        # If an output directory wasn't provided, we put the no duplicate read file in the same
        # directory as the normal read file.
        else:
            no_dup_filename = os.path.join(os.path.dirname(os.path.abspath(filename)),
                                           no_dup_filename)

        if not silent:
            log.log('\nDuplicate read names found. Saving duplicate-free file:')
            log.log(no_dup_filename)
        with gzip.open(no_dup_filename, 'wb') as f:
            for read_name in read_names:
                read = read_dict[read_name]
                if file_type == 'FASTQ':
                    f.write(read.get_fastq().encode())
                else:  # file_type == 'FASTA'
                    f.write(read.get_fasta().encode())

    else:
        no_dup_filename = filename

    return read_dict, read_names, no_dup_filename


class Reference(object):
    """
    This class holds a reference sequence: just a name and a nucleotide sequence.
    """

    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence.upper()

        # If the reference name also happens to be a number, store it as an int.
        try:
            self.number = int(name)
        except ValueError:
            self.number = 0

    def __repr__(self):
        return self.name + ' (' + str(len(self.sequence)) + ' bp)'

    def get_length(self):
        """
        Returns the sequence length.
        """
        return len(self.sequence)


class Read(object):
    """
    This class holds a long read, e.g. from PacBio or Oxford Nanopore.
    """

    def __init__(self, name, sequence, qualities):
        self.name = name
        self.sequence = sequence.upper()

        if qualities:
            self.qualities = qualities

        # If no qualities are given, then they are all set to '+', the Phred+33 score for 10% error.
        else:
            self.qualities = '+' * len(self.sequence)

        self.alignments = []

    def __repr__(self):
        return self.name + ' (' + str(len(self.sequence)) + ' bp)'

    def get_length(self):
        """
        Returns the sequence length.
        """
        return len(self.sequence)

    def remove_conflicting_alignments(self, allowed_overlap):
        """
        This function removes alignments from the read which are likely to be spurious or
        redundant.
        """
        self.alignments = sorted(self.alignments, reverse=True,
                                 key=lambda x: (x.raw_score, random.random()))
        kept_alignments = []
        kept_alignment_ranges = []
        for alignment in self.alignments:
            this_range = alignment.read_start_end_positive_strand()

            # Don't keep alignments for which their part of the read is already aligned.
            if range_is_contained(this_range, kept_alignment_ranges):
                continue

            # Don't keep alignments which overlap too much with existing alignments.
            if range_overlap_size(this_range, kept_alignment_ranges) > allowed_overlap:
                continue

            # Don't keep alignments that seem to be very similar to an already kept alignment.
            keep_alignment = True
            for kept_alignment in kept_alignments:
                if kept_alignment.is_very_similar(alignment):
                    keep_alignment = False
                    break

            if keep_alignment:
                kept_alignments.append(alignment)
                kept_alignment_ranges = simplify_ranges(kept_alignment_ranges + [this_range])

        kept_alignments = sorted(kept_alignments,
                                 key=lambda x: x.read_start_end_positive_strand()[0])
        self.alignments = kept_alignments

    def remove_low_score_alignments(self, low_score_threshold):
        """
        This function removes alignments with identity below the cutoff.
        """
        self.alignments = [x for x in self.alignments
                           if x.scaled_score is not None and x.scaled_score >= low_score_threshold]

    def remove_short_alignments(self, min_align_length):
        """
        This function removes alignments with identity below the cutoff.
        """
        self.alignments = [x for x in self.alignments
                           if x.get_aligned_ref_length() >= min_align_length]

    def get_fastq(self):
        """
        Returns a string for the read in FASTQ format. It contains four lines and ends in a line
        break.
        """
        return '@' + self.name + '\n' + \
               self.sequence + '\n' + \
               '+\n' + \
               self.qualities + '\n'

    def get_fasta(self):
        """
        Returns a string for the read in FASTA format (ending in a line break).
        """
        return '>' + self.name + '\n' + add_line_breaks_to_sequence(self.sequence, 70)

    def get_fraction_aligned(self):
        """
        This function returns the fraction of the read which is covered by any of the read's
        alignments.
        """
        if len(self.sequence) == 0:
            return 0.0
        read_ranges = [x.read_start_end_positive_strand()
                       for x in self.alignments]
        read_ranges = simplify_ranges(read_ranges)
        aligned_length = sum([x[1] - x[0] for x in read_ranges])
        return aligned_length / len(self.sequence)

    def get_reference_bases_aligned(self):
        """
        This function returns the number of bases aligned with respect to the reference.
        """
        return sum([x.get_aligned_ref_length() for x in self.alignments])

    def has_one_contained_alignment(self):
        """
        Returns true if this read aligned entirely within a reference (i.e. no read end gaps).
        """
        return len(self.alignments) == 1 and \
            self.alignments[0].read_start_pos == 0 and \
            self.alignments[0].read_end_gap == 0

    def mostly_aligns_to_contamination(self):
        """
        Returns true if 50% or more of the alignments are to contaminant sequences.
        """
        if len(self.sequence) == 0:
            return False
        if not self.alignments:
            return False
        contamination_alignment_length = sum(x.get_aligned_read_length() for x in self.alignments
                                             if x.ref.name.startswith('CONTAMINATION_'))
        good_alignment_length = sum(x.get_aligned_read_length() for x in self.alignments
                                    if not x.ref.name.startswith('CONTAMINATION_'))
        return contamination_alignment_length >= good_alignment_length

    def aligns_to_multiple_single_copy_segments(self, single_copy_segment_names):
        return sum(x.ref.name in single_copy_segment_names for x in self.alignments) > 1

    def get_alignment_table(self):
        alignment_table = [['Ref name', 'Ref start', 'Ref end', 'Read start', 'Read end', 'Strand',
                            'Raw score', 'Scaled score', 'Identity']]
        for alignment in self.alignments:
            read_start, read_end = alignment.read_start_end_positive_strand()
            strand = '-' if alignment.rev_comp else '+'
            ref_name = alignment.ref.name
            if ref_name.startswith('CONTAMINATION'):
                ref_name = 'CONTAM'
            alignment_row = [ref_name, str(alignment.ref_start_pos), str(alignment.ref_end_pos),
                             str(read_start), str(read_end), strand]
            if alignment.scaled_score is not None:
                alignment_row += [str(alignment.raw_score), float_to_str(alignment.scaled_score, 2)]
            else:
                alignment_row += ['', '']
            if alignment.percent_identity is not None:
                alignment_row.append(float_to_str(alignment.percent_identity, 2) + '%')
            else:
                alignment_row.append('')
            alignment_table.append(alignment_row)
        return print_table(alignment_table, alignments='RRRRRRRRR', return_str=True,
                           header_format=None, col_separation=2, indent=2)


def get_read_nickname_dict(read_names):
    """
    Read names can be quite long, so for the sake of output brevity, this function tries to come
    up with some shorter nicknames for the reads.
    """
    max_read_name_len = max(len(name) for name in read_names)
    for nickname_length in range(1, max_read_name_len):
        nicknames = set()
        for name in read_names:
            nickname = name[:nickname_length]
            if nickname in nicknames:
                break
            nicknames.add(nickname)
        else:
            return {name: name[:nickname_length] for name in read_names}

    # If we couldn't find a length for shorter nicknames, then the nicknames are just the full
    # names. Oh well.
    return {name: name for name in read_names}
