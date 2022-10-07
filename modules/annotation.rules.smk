from datetime import datetime
import re
import os

date = datetime.now().strftime('%Y%m%d.%H%M%S')

rule MLST:
  input:
    assembly = "assembly.fasta"
  output: 
    mlst = "tables/mlst.tbl"
  params:
    scripts_dir = "../scripts/",
    sample = "sample",
    pipeline_version = "v1"
  threads: 1
  conda:
    "../envs/mlst2.22.1.yaml"
  shell:
    "mlst {input.assembly} | {params.scripts_dir}process_mlst.pl {params.sample}_{params.pipeline_version} > {output.mlst};"

rule Prokka:
  input:
    assembly = "assembly.fasta"
  output:
    gff = "03_Annotation/prokka.gff3",
    faa = "03_Annotation/prokka.faa"
  params:
    annotation_dir = "03_Annotation/",
    prokka_out = "03_Annotation/prokka/",
    genus = "Eschericcia",
    species = "coli",
    prefix = "Ecoli"
  threads: 4
  conda:
    "../envs/prokka1.14.6.yaml"
  shell:
    "prokka --addgenes --addmrna --cpus {threads} --norrna --outdir {params.prokka_out} --kingdom Bacteria --genus {params.genus} --species {params.species} {input.assembly};"
    "ln -s {params.prokka_out}/*.gff {output.gff};"
    "ln -s {params.prokka_out}/*.faa {output.faa};"

rule RGI:
  input:
    assembly = "assembly.fasta",
    prokka_gff = "03_Annotation/prokka.gff3",
    prokka_faa = "03_Annotation/prokka.faa"
  output:
    rgi_out = "rgi.txt",
    prokka_tbl = "annot.tbl",
    rgi_tbl = "rgi.tbl",
    annotation_4jb = "annot.4jb.gff3",
    proteins = "annotation.pep.fa",
  params:
    annotation_dir = "03_Annotation/",
    prefix = "Ecoli",
    scripts_dir = "../scripts/",
  threads: 4
  conda:
    "../envs/RGI6.0.0.yaml"
  shell:
    "cd `dirname {output.rgi_out}`;"
    "echo \"Running RGI. Database version:\";"
    "rgi database -v;"
    "rgi main -n {threads} -t contig -i {input.assembly} -o rgi;"
    "cd {params.annotation_dir};"
    "{params.scripts_dir}join_annotations.pl -prokka {input.prokka_gff} -pep {input.prokka_faa} -rgi {output.rgi_out} -pref {params.prefix};"
    "cut -f 5,12 {output.prokka_tbl} | grep -v protein_sequence |  {params.scripts_dir}TblToFasta > {output.proteins};"

rule AMRFinder:
  input:
    annotation =  "annot.tbl",
    annotation_gff = "annot.4jb.gff3",
    proteins = "annotation.pep.fa",
  output:
    amrfinder =  "amrfinder.tsv",
    amrfinder_gff = "amrfinder.4jb.gff3"
  params:
    scripts_dir = "../scripts/",
    additional = ""
  conda:
    "../envs/amrfinder3.10.42.yaml"    
  threads: 1
  shell:
    "amrfinder {params.additional} --threads {threads} -p {input.proteins} | {params.scripts_dir}change_header_amrfinder.py > {output.amrfinder};"
    "{params.scripts_dir}parseAMRFinder4jb.pl -annot {input.annotation_gff} -amrf {output.amrfinder} > {output.amrfinder_gff};"

rule resfinder:
  input:
    assembly = "assembly.fasta",
    db = "databases/resfinder_db/"
  output:
    gff = "resfinder_results.gff3",
    tbl =  "ResFinder_results_tab.txt",
    pheno = "pheno_table.txt"
  params:
    resfinder_dir = "resfinder_out/",
    score = "0.8",
    scripts_dir = "../scripts/"
  conda:
    "../envs/resfinder4.1.11.yaml"
  threads: 1
  shell: 
    "mkdir -p {params.resfinder_dir};"
   # "run_resfinder.py -ifa {input.assembly} -o {params.resfinder_dir}  -t {params.score} -acq;"
    "run_resfinder.py -ifa {input.assembly} -o {params.resfinder_dir} -db_res {input.db} -t {params.score} -acq;"
    "{params.scripts_dir}parse_resfinder_output.pl {params.resfinder_dir}/resfinder_blast/tmp > {output.gff};"

rule IS_Finder:
  input:
    assembly = "assembly.fasta",
    db = "databases/ISFinder/ISfinder_genes.fasta"
  output:
    ISF_out = "ISFinder_results.gff3"
  params:
    IS_dir = "03_Annotation"
  conda:
    "../envs/BLAST2.13.0.yaml"
  threads: 4
  shell:
    "TMPDIR={params.IS_dir}/tmp_isfinder;"
    "mkdir -p $TMPDIR;"
    "cd $TMPDIR;"
    "ln -s {input.assembly} assembly.fa;"
    "makeblastdb -in assembly.fa -dbtype nucl -parse_seqids;"
    "blastn -db assembly.fa -query {input.db} -out {params.IS_dir}/ISFinder.out -outfmt 7 -num_threads {threads};"
    "cd {params.IS_dir};"
    "grep -v \'#\' ISFinder.out | sort -k2,2 | perl -ane \'chomp; $c++; if ($F[8]<$F[9]){{print \"$F[1]\\tISFinder\\tIS_element\\t$F[8]\\t$F[9]\\t$F[11]\\t+\\t.\\tID=IS.$c;Target=$F[0]:$F[6]-$F[7];Name=$F[0];Evalue=$F[10]\\n\";}}else {{print \"$F[1]\\tISFinder\\tIS_element\\t$F[9]\\t$F[8]\\t$F[11]\\t-\\t.\\tID=IS.$c;Target=$F[0]:$F[6]-$F[7];Name=$F[0];Evalue=$F[10]\\n\";}}\' > {output.ISF_out};"
    "rm -r $TMPDIR;"

rule centrifuge:
  input:
    assembly = "assembly.fasta",
    species_verification = "species_verification.out",
    scaffolds_table = "02_assembly/unicycler_out/assembly.fasta",
  output:
    out = "centrifuge.out",
    report =  "centrifuge.report",
    table = "centrifuge.tbl"
  params:
    centrifuge_db = "/scratch/devel/jgomez/centrifuge_cre_database/bacteria/bact",
    centrifuge_k = 5,
    scripts_dir = "../scripts/"
  conda:
    "../envs/centrifuge1.0.4.yaml"
  threads: 8
  shell:
    "centrifuge -k {params.centrifuge_k} -x {params.centrifuge_db} -f {input.assembly} -p {threads} -S {output.out} --report-file {output.report};"
    "echo \"scaffold\tcircular\tcentrifuge_sp\tseqID\" > {output.table};"
    "{params.scripts_dir}combine_centrifuge.pl {output.out} {output.report} {input.species_verification} {input.scaffolds_table} | cut -f 1,2,3,4 >> {output.table};"