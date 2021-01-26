"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module contains functions relating to producing the VCF with the final Unicycler assembly.

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
import shutil
from .misc import dim
from . import log


class CannotMakeVcf(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


def make_vcf(vcf_filename, args, assembly_file, insert_size_1st, insert_size_99th):
    log.log_section_header('Building VCF')

    if insert_size_1st is None or insert_size_99th is None:
        raise CannotMakeVcf('no insert size range')

    vcf_dir = os.path.join(args.out, 'vcf')
    if not os.path.exists(vcf_dir):
        os.makedirs(vcf_dir)

    # To avoid issues with paths that contain spaces, we will move into the temporary directory
    # to run these commands.
    starting_dir = os.getcwd()
    os.chdir(vcf_dir)

    input_fasta = os.path.join(vcf_dir, 'assembly.fasta')
    using_paired_reads = bool(args.short1) and bool(args.short2)
    using_unpaired_reads = bool(args.unpaired)
    sam_filename = 'alignments.sam'
    bam_filename = 'alignments.bam'

    shutil.copyfile(assembly_file, input_fasta)

    # Prepare the FASTA for Bowtie2 alignment.
    bowtie2_build_command = [args.bowtie2_build_path, input_fasta, input_fasta]
    log.log(dim('  ' + ' '.join(bowtie2_build_command)), 2)
    try:
        subprocess.check_output(bowtie2_build_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotMakeVcf('bowtie2-build encountered an error:\n' + e.output.decode())
    if not any(x.endswith('.bt2') for x in os.listdir(vcf_dir)):
        raise CannotMakeVcf('bowtie2-build failed to build an index')

    # Perform the alignment with Bowtie2.
    bowtie2_command = [args.bowtie2_path, '--local', '--very-sensitive-local',
                       '--threads', str(args.threads), '-I', str(insert_size_1st),
                       '-X', str(insert_size_99th), '-x', input_fasta, '-S', sam_filename]
    if using_paired_reads:
        bowtie2_command += ['-1', args.short1, '-2', args.short2]
    if using_unpaired_reads:
        bowtie2_command += ['-U', args.unpaired]
    log.log(dim('  ' + ' '.join(bowtie2_command)), 2)
    try:
        subprocess.check_output(bowtie2_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotMakeVcf('Bowtie2 encountered an error:\n' + e.output.decode())

    # Sort the alignments.
    samtools_sort_command = [args.samtools_path, 'sort', '-@', str(args.threads),
                             '-o', bam_filename, '-O', 'bam', '-T', 'temp', sam_filename]
    log.log(dim('  ' + ' '.join(samtools_sort_command)), 2)
    try:
        subprocess.check_output(samtools_sort_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotMakeVcf('Samtools encountered an error:\n' + e.output.decode())

    # Index the alignments.
    samtools_index_command = [args.samtools_path, 'index', bam_filename]
    log.log(dim('  ' + ' '.join(samtools_index_command)), 2)
    try:
        subprocess.check_output(samtools_index_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotMakeVcf('Samtools encountered an error:\n' + e.output.decode())

    # Index the FASTA.
    faidx_command = [args.samtools_path, 'faidx', input_fasta]
    log.log(dim('  ' + ' '.join(faidx_command)), 2)
    try:
        subprocess.check_output(faidx_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotMakeVcf('samtools faidx encountered an error:\n' + e.output.decode())

    mpileup_command = args.samtools_path + ' mpileup -v -u -f ' + input_fasta + ' ' + \
        bam_filename + ' | ' + args.bcftools_path + ' call -c -v -o ' + vcf_filename
    log.log(mpileup_command)
    mpileup = subprocess.Popen(mpileup_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               shell=True)
    out, err = mpileup.communicate()
    if mpileup.returncode != 0:
        raise CannotMakeVcf('mpileup/bcftools encountered an error: ' + err.decode())

    if not os.path.isfile(vcf_filename):
        raise CannotMakeVcf('failed to make VCF file')

    os.chdir(starting_dir)
    if args.keep < 3 and os.path.exists(vcf_dir):
        shutil.rmtree(vcf_dir, ignore_errors=True)
