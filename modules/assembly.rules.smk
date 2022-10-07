from datetime import datetime
import re
import os

date = datetime.now().strftime('%Y%m%d.%H%M%S')

rule unicycler:
  input:
    read1 = "reads.illumina.1.fastq.gz",
    read2 =  "reads.illumina.2.fastq.gz",
    nanopore ="reads.ont.fastq.gz",
  output:
    fasta = "02_assembly/unicycler_out/assembly.fasta",
    plot = "plots/assembly.png",
  params:
    assemblydir = "02_assembly",
    start_genes = "start_genes.fasta",
  threads: 12
  conda:
    "../envs/unicycler0.5.0.yaml"
  envmodules:
    'Mesa/21.1.7-GCCcore-11.2.0'
  shell:    
    "export HDF5_USE_FILE_LOCKING=FALSE;"
    "XDG_RUNTIME_DIR=/scratch/groups/assembly/shared/projects/BALSALOBRE/BALSALOBRE_02_03/assembly/xdg_runtime_dir;"
    "mkdir -p {params.assemblydir};"
    "cd {params.assemblydir}; "
    "unicycler -1 {input.read1} -2 {input.read2}  -l {input.nanopore} -o unicycler_out -t {threads} --start_genes {params.start_genes}; "
    "mkdir -p plots ; "
    "Bandage image unicycler_out/assembly.gfa {output.plot}; "

rule produce_final_assembly:
  input:
    fasta = "02_assembly/unicycler_out/assembly.fasta",
    plot = "plots/assembly.png",
  output:
    fasta = "assembly.fasta",
    table1 = "tables/assembly.tbl",
    table2 = "tables/scaffolds.tbl",
    nseries = "tables/assembly.nseries.txt"
  params:
    assemblydir = "02_assembly",
    scripts_dir = "scripts/",
    sample = "sample",
    pipeline_version = "v1"
  conda:
    "../envs/ass_base.yaml"
  shell:  
    "cd {params.assemblydir};"
    "mkdir -p tables; "
    "perl {params.scripts_dir}unicycler_postprocessing.v02.pl -f {input.fasta} -s {params.sample} -p {params.pipeline_version} -a Unicycler_v0.5.0;"
    "{params.scripts_dir}fastalength {output.fasta} | {params.scripts_dir}Nseries.pl > {output.nseries}; "
