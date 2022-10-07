from datetime import datetime
import re
import os

date = datetime.now().strftime('%Y%m%d.%H%M%S')

rule trim_galore:
  input:
    read1 = "illumina.1.fastq.gz",
    read2 = "illumina.2.fastq.gz",
  output:
    trim1 = "illumina.trimed.1.fastq.gz",
    trim2 = "illumina.trimed.2.fastq.gz",
  params:
    outdir = "illumina_trim",
    opts = "--gzip -q 20 --paired --retain_unpaired",
  threads: 4
  conda:
    "../envs/trim_galore0.6.7.yaml"
  shell:
    "mkdir -p {params.outdir};"
    "cd {params.outdir}; "
    "trim_galore -j {threads} {params.opts} {input.read1} {input.read2} ;"
    "b=`basename {input.read1} .1.fastq.gz`;"
    "ln -s $b.1_val_1.fq.gz {output.trim1};"
    "ln -s $b.2_val_2.fq.gz {output.trim2};"

rule Trim_ONT_qcat:
  input:
    read = "ont.fastq.gz"
  output:
    trim = "ont.trim.fastq.gz"   
  params:
    options = "-b ont_trim/tmp/trimmed.fastq",
    workdir = "ont_trim",
    tempdir = "ont_trim/tmp",
    tempfile = "ont_trim/tmp/trimmed.fastq"
  threads: 4
  conda:
    "../envs/qcat1.1.0.yaml"
  shell:
    "mkdir -p {params.tempdir} ;"
    "cd {params.workdir}; "
    "zcat {input.read} | qcat {params.options} --detect-middle --trim ; "
    "gzip -c {params.tempfile} > {output.trim} ; "
    "rm -r {params.tempdir};"

rule concat_reads:
  input:
    fastqs = ["reads1.fastq.gz", "reads2.fastq.gz"]
  output:
    final_fastq = "all_reads.fastq.gz"
  conda:
    '../envs/ass_base.yaml'
  threads: 2
  shell:
    "zcat {input.fastqs} | pigz -p {threads} -c  > {output.final_fastq};"

rule nanoplot:
  input:
    fastq = "reads.ont.fastq.gz"
  output:
    stats = "nanostats_out/NanoStats.txt",
    yield_len = "nanostats_out/Yield_By_Length.png",
    read_len = "nanostats_out/WeightedHistogramReadlength.png"
  params:
    outdir = "nanostats_out/",
  log:
    "{dir}logs/" + str(date) + ".NanoStats.out",
    "{dir}logs/" + str(date) + ".NanoStats.err"
  benchmark:
    "{dir}logs/" + str(date) + ".NanoStats.benchmark.txt"
  threads: 4
  conda:
    '../envs/nanoplot1.40.0'
  shell:
    "mkdir -p {params.outdir};" 
    "cd {params.outdir}; "
    "NanoPlot -t {threads} --plots dot --fastq {input.fastq} -o .; "

rule Kraken2:
  input:
    read = "reads.fastq.gz",
    database = "/scratch/groups/assembly/shared/databases/kraken2/bacteria_db/",
    kmers = "/scratch/groups/assembly/shared/databases/kraken2/bacteria_db/database100mers.kmer_distrib"
  output:
    report = "kraken2.report",
    readsout = "kraken2.seqs.out",
    abundance =  "bracken_abundance.txt",
    tophits= "bracken_abundance.tophits.txt"
  params:
    additional = "",
    prefix = "assembly",
    scripts_dir = "../scripts"
  threads: 4
  conda:
    '../envs/kraken2.1.2.yaml'
  shell:
     "kraken2 --threads {threads} --db {input.database}  --use-names --report {output.report} {params.additional} {input.read} > {output.readsout}; "
     "est_abundance.py -i {output.report} -k {input.kmers} -l S -t 10 -o {output.abundance}; "
     "{params.scripts_dir}bracken-top-hits.v01.pl -f  {output.abundance} > {output.tophits}; "

rule Verify_Species:
  input: 
    illumina_tophit = "illumina.bracken_abundance.tophits.txt",
    nanopore_tophit = "ont.bracken_abundance.tophits.txt",
  output:
    verification = "species_verification.out"
  params:
    workdir = "01_preprocessing",
    metadata_species = "\"Eschericcia coli\"",
    scripts_dir = "../scripts"
  threads: 1
  shell:
    "cd {params.workdir}; "
    "{params.scripts_dir}verify_species.v01.pl -i {input.illumina_tophit} -n {input.nanopore_tophit} -l {params.metadata_species};"