#Authors: Jessica Gomez-Garrido, CNAG-CRG; Fernando Cruz, CNAG-CRG.
#Contact email: jessica.gomez@cnag.crg.eu
#Date:20190704
#Snakemake pipeline to assemble and annotate bacterial genomes from ONT and Illumina genomic data

from datetime import datetime

date = datetime.now().strftime('%Y%m%d.%H%M%S')

scripts_dir = config["Inputs"]["scripts_dir"]
shell.prefix("export PATH=" + scripts_dir + ":$PATH;")
logs_dir = config["Parameters"]["logs_dir"]
if not os.path.exists(logs_dir):
    os.makedirs(logs_dir)
base_dir= config["Parameters"]["basedir"]

# file wildcards 
illumina_files = config["Wildcards"]["ILLUMINA_fastqs"]
if illumina_files == None:
  illumina_files = ""
  print ("WARNING no Illumina fastqs have been given")

ont_files = config["Wildcards"]["ONT_fastqs"]
if ont_files == None:
  ont_files = ""
  print ("WARNING no ONT fastqs have been given")

#Global variables definition
sample = config["Parameters"]["sample"]
pipeline_version = "v" + str(config["Parameters"]["version"])
species= config["Parameters"]["metadata_species"].split(" ")
databases = scripts_dir + "../databases/"

#Rule for outputs
rule all:
  input:
    expand(config["Outputs"]["ILLUMINA_trim"] + "{file}.1_val_1.fq.gz", file=illumina_files.split(',')),
    expand(config["Outputs"]["ILLUMINA_trim"] + "{file}.2_val_2.fq.gz", file=illumina_files.split(',')),
    config["Outputs"]["ILLUMINA_cat"] + sample + ".1.fastq.gz",
    config["Outputs"]["ILLUMINA_cat"] + sample + ".2.fastq.gz",
    config["Outputs"]["ILLUMINA_bracken"] + sample + ".illumina.kraken2.report",
    config["Outputs"]["ILLUMINA_bracken"] + sample + ".illumina.bracken_abundance.txt",
    config["Outputs"]["ILLUMINA_bracken"] + sample + ".illumina.bracken_abundance.tophits",
    expand(config["Outputs"]["ONT_trim"] + "{ontfile}.trim.fastq.gz", ontfile=ont_files.split(',')),
    expand(config["Outputs"]["nanostats_dir"] +"{ontfile}.HistogramReadlength.png", ontfile=ont_files.split(',')),
    expand(config["Outputs"]["nanostats_dir"] +"{ontfile}.LengthvsQualityScatterPlot_dot.png", ontfile=ont_files.split(',')),
    expand(config["Outputs"]["nanostats_dir"] +"{ontfile}.LogTransformed_HistogramReadlength.png", ontfile=ont_files.split(',')),
    expand(config["Outputs"]["nanostats_dir"] +"{ontfile}.NanoPlot-report.html", ontfile=ont_files.split(',')),
    expand(config["Outputs"]["nanostats_dir"] +"{ontfile}.NanoStats.txt", ontfile=ont_files.split(',')),
    config["Outputs"]["ONT_cat"] + sample + ".ont.fastq.gz",
    config["Outputs"]["ONT_bracken"] + sample + ".ont.kraken2.report",
    config["Outputs"]["ONT_bracken"] + sample + ".ont.bracken_abundance.txt",
    config["Outputs"]["ONT_bracken"] + sample + ".ont.bracken_abundance.tophits",
    config["Outputs"]["preprocessing_dir"] + sample + ".species_verification.out",
    config["Outputs"]["assembly"],
    config["Outputs"]["assembly_dir"] + "plots/" + sample + "." + pipeline_version + ".assembly.png",
    config["Outputs"]["assembly_dir"] + "tables/" + sample + "." + pipeline_version +  ".assembly.tbl",
    config["Outputs"]["assembly_dir"] + "tables/" + sample + "." + pipeline_version +  ".scaffolds.tbl",
    config["Outputs"]["assembly_dir"] + "tables/" + sample + "." + pipeline_version +  ".assembly.nseries.txt",
    config["Outputs"]["assembly_dir"] + "tables/" + sample + "." + pipeline_version +  ".mlst.tbl",
    config["Outputs"]["centrifuge_dir"] + sample + "." + pipeline_version + ".out",
    config["Outputs"]["centrifuge_dir"] + sample + "." + pipeline_version + ".report",
    config["Outputs"]["centrifuge_dir"] + sample + "." + pipeline_version + ".tbl",
    config["Outputs"]["annotation_dir"] + sample + pipeline_version +  ".annot.tbl",
    config["Outputs"]["annotation_dir"] + sample + pipeline_version + ".rgi.tbl",
    config["Outputs"]["annotation_dir"] + sample + pipeline_version + ".annot.4jb.gff3",
    config["Outputs"]["annotation_dir"] + sample + pipeline_version + ".pep.fa",
    config["Outputs"]["annotation_dir"] + sample + pipeline_version + ".amrfinder.tsv",
    config["Outputs"]["annotation_dir"] + sample + pipeline_version + ".amrfinder.4jb.gff3",
    config["Outputs"]["resfinder_gff"],
    config["ResFinder"]["dir"] + "ResFinder_results_tab.txt",
    config["ResFinder"]["dir"] + "pheno_table.txt",
    config["Outputs"]["ISFinder_out"]
  log:
    logs_dir + str(date) + ".rule_all.log"

#Pipeline rules
rule Trim_Galore:
  input:
    read1 = config["Inputs"]["ILLUMINA_reads"] + "{file}.1.fastq.gz",
    read2 = config["Inputs"]["ILLUMINA_reads"] + "{file}.2.fastq.gz"
  output:
    trim1 = config["Outputs"]["ILLUMINA_trim"] + "{file}.1_val_1.fq.gz",
    trim2 = config["Outputs"]["ILLUMINA_trim"] + "{file}.2_val_2.fq.gz",
  params:
    outdir = directory(config["Outputs"]["ILLUMINA_trim"]),
    quality = config["Trim_Galore"]["quality_threshold"],
    length = config["Trim_Galore"]["illumina_RL"]
  log: 
    logs_dir + str(date) + ".{file}.trim_galore.log"
  threads: config["Trim_Galore"]["Trim_Illumina_cores"]
  shell:
    "mkdir -p {params.outdir};"
    "cd {params.outdir}; "
    "trim_galore --gzip -q {params.quality} --length {params.length} --paired {input.read1} {input.read2} ;"

rule Concat_Illumina:
  input:
    trim1 = lambda wildcards: expand(rules.Trim_Galore.output.trim1, file=illumina_files.split(',') ),
    trim2 = lambda wildcards: expand(rules.Trim_Galore.output.trim2,  file=illumina_files.split(','))
  output:
    concat1 = config["Outputs"]["ILLUMINA_cat"] + sample + ".1.fastq.gz",
    concat2 = config["Outputs"]["ILLUMINA_cat"] + sample + ".2.fastq.gz"
  params:
    outdir =config["Outputs"]["ILLUMINA_cat"] 
  log:
    logs_dir + str(date) + ".illumina_concatenation.log"
  threads: config["Trim_Galore"]["Concat_Illumina_cores"]
  shell:
    "mkdir -p {params.outdir}; "
    "zcat {input.trim1} | pigz -p {threads} -c  > {output.concat1}; "
    "zcat {input.trim2} | pigz -p {threads} -c  > {output.concat2}; "

rule Bracken_Illumina:
  input:
    read1 = rules.Concat_Illumina.output.concat1,
    read2 = rules.Concat_Illumina.output.concat2
  output:
    report = config["Outputs"]["ILLUMINA_bracken"] + sample + ".illumina.kraken2.report",
    abundance =  config["Outputs"]["ILLUMINA_bracken"] + sample + ".illumina.bracken_abundance.txt",
    tophits =  config["Outputs"]["ILLUMINA_bracken"] + sample + ".illumina.bracken_abundance.tophits"
  params:
    database = config["Bracken"]["kraken_db"],
    kmers = config["Bracken"]["kraken_kmers"]
  log:
    logs_dir + str(date) + ".illumina_bracken.log"
  threads: config["Bracken"]["brackenIlluminaCores"]
  shell:
   # Running kraken2 and bracken simultaneously using both illumina fastq
    "kraken2 --threads {threads} --db {params.database}  <(zcat {input.read1} {input.read2}) --use-names --report {output.report}; "
    "wait; "
    "est_abundance.py -i {output.report} -k {params.kmers} -l S -t 10 -o {output.abundance}; "
    "bracken-top-hits.v01.pl -f  {output.abundance} > {output.tophits}; "

rule Trim_ONT_qcat:
  input:
    read = config["Inputs"]["ONT_reads"] + "{ontfile}.1.fastq.gz"
  output:
    trim = config["Outputs"]["ONT_trim"] + "{ontfile}.trim.fastq.gz"   
  params:
    barcode = lambda wildcards: config["ONT_Barcodes"][wildcards.ontfile],
    workdir = config["Outputs"]["ONT_trim"],
    tempdir = config["Outputs"]["ONT_trim"] + "temp_{ontfile}",
  log:
    logs_dir + str(date) + ".{ontfile}.qcat.log"   
  threads: config["Parameters"]["QCAT_cores"]
  shell:
    "mkdir -p {params.tempdir} ;"
    "cd {params.workdir}; "
    "zcat {input.read} | qcat -b {params.tempdir} --detect-middle --trim ; "
    "gzip -c {params.tempdir}/barcode{params.barcode}.fastq > {output.trim} ; "
    "rm -r {params.tempdir};"

rule NanoStats_ONT:
  input:
    trim = lambda wildcards: expand(rules.Trim_ONT_qcat.output.trim, ontfile=ont_files.split(',') )
  output:
    histogram = config["Outputs"]["nanostats_dir"] +  "{ontfile}.HistogramReadlength.png", 
    lengthPlot = config["Outputs"]["nanostats_dir"] + "{ontfile}.LengthvsQualityScatterPlot_dot.png", 
    logPlot = config["Outputs"]["nanostats_dir"] + "{ontfile}.LogTransformed_HistogramReadlength.png",
    html = config["Outputs"]["nanostats_dir"] +  "{ontfile}.NanoPlot-report.html",
    stats = config["Outputs"]["nanostats_dir"] + "{ontfile}.NanoStats.txt", 
  params:
    outdir = config["Outputs"]["nanostats_dir"]
  log:
    logs_dir + str(date) + "{ontfile}.NanoPlot.log"
  threads: config["Parameters"]["nanostat_cores"]
  shell:
    "source ~/init_shell.sh;"
    "conda deactivate;"
    "conda activate {scripts_dir}/../conda_environment_CRE2;"
    "mkdir -p {params.outdir}; "
    "cd {params.outdir}; "
    "NanoPlot -t {threads} --plots dot -p {wildcards.ontfile}. --fastq {input.trim} -o .; "   

rule Concat_ONT:
  input:
    trim = lambda wildcards: expand(rules.Trim_ONT_qcat.output.trim, ontfile=ont_files.split(',') )
  output:
    concat = config["Outputs"]["ONT_cat"] + sample + ".ont.fastq.gz"
  params:
    outdir = config["Outputs"]["ONT_cat"]
  log:
    logs_dir + str(date) + ".ont_concatenation.log"
  threads: config["Parameters"]["Concat_ONT_cores"]
  shell:
    "mkdir -p {params.outdir}; "
    "zcat {input.trim} | pigz -p {threads} -c  > {output.concat}; "

rule Bracken_ONT:
  input:
    read = rules.Concat_ONT.output.concat
  output:
    report = config["Outputs"]["ONT_bracken"] + sample + ".ont.kraken2.report",
    abundance =  config["Outputs"]["ONT_bracken"] + sample + ".ont.bracken_abundance.txt",
    tophits =  config["Outputs"]["ONT_bracken"] + sample + ".ont.bracken_abundance.tophits"
  params:
    database = config["Bracken"]["kraken_db"],
    kmers = config["Bracken"]["kraken_kmers"]
  log:
    logs_dir + str(date) + ".ont_bracken.log"
  threads: config["Bracken"]["brackenONTCores"]
  shell:
   # Running kraken2 and bracken simultaneously using both illumina fastq
     "kraken2 --threads {threads} --db {params.database}  <( zcat {input.read} ) --use-names --report {output.report}; "
     "wait; "
     "est_abundance.py -i {output.report} -k {params.kmers} -l S -t 10 -o {output.abundance}; "
     "bracken-top-hits.v01.pl -f {output.abundance} > {output.tophits}; "

rule Verify_Species:
  input: 
    illumina_tophit = rules.Bracken_Illumina.output.tophits,
    nanopore_tophit = rules.Bracken_ONT.output.tophits
  output:
    verification = config["Outputs"]["preprocessing_dir"] + sample + ".species_verification.out"
  params:
    workdir = config["Outputs"]["preprocessing_dir"],
    metadata_species = "\"" + config["Parameters"]["metadata_species"] + "\"",
  log:
    logs_dir + str(date) + ".verify_species.log"
  threads: 1
  shell:
    "cd {params.workdir}; "
    "verify_species.v01.pl -i {input.illumina_tophit} -n {input.nanopore_tophit} -l {params.metadata_species};"

rule Assembly:
  input:
    read1 = rules.Concat_Illumina.output.concat1,
    read2 = rules.Concat_Illumina.output.concat2,
    nanopore = rules.Concat_ONT.output.concat
  output:
    fasta = config["Outputs"]["assembly"],
    plot = config["Outputs"]["assembly_dir"] + "plots/" + sample + "." + pipeline_version + ".assembly.png",
    table1 = config["Outputs"]["assembly_dir"] + "tables/" + sample + "." + pipeline_version +  ".assembly.tbl",
    table2 = config["Outputs"]["assembly_dir"] + "tables/" + sample + "." + pipeline_version +  ".scaffolds.tbl",
    nseries = config["Outputs"]["assembly_dir"] + "tables/" + sample + "." + pipeline_version +  ".assembly.nseries.txt"
  params:
    assemblydir = config["Outputs"]["assembly_dir"],
    pilon_path = config["Assembly"]["pilon_path"],
  log:
    logs_dir + str(date) + ".assembly.log"
  threads: config["Assembly"]["assemblyCores"]
  shell:    
    "export HDF5_USE_FILE_LOCKING=FALSE;"
    "cd {params.assemblydir}; "
    "unicycler --pilon_path {params.pilon_path} -1 {input.read1} -2 {input.read2}  -l {input.nanopore} -o unicycler_out -t {threads} --start_genes {databases}/unicycler_gene_data/start_genes.fasta; "
    
    #convert assembly graph into PNG file
    "mkdir -p plots ; "
    "Bandage image unicycler_out/assembly.gfa {output.plot}; "

    #Create assembly table for database
    # get UniCycler version by capturing the standard error
    "univer=$(unicycler --version 2>&1 | sed 's/ /_/g') ; "
    "mkdir -p tables; "
    "unicycler_postprocessing.v02.pl -f unicycler_out/assembly.fasta -s {sample} -p {pipeline_version} -a $univer;"
    "fastalength {output.fasta} | Nseries.pl > {output.nseries}; "

rule MLST:
  input:
    assembly = rules.Assembly.output.fasta
  output: 
    mlst = config["Outputs"]["assembly_dir"] + "tables/" + sample + "." + pipeline_version +  ".mlst.tbl"
  params:
  log:
    logs_dir + str(date) + ".mlst.out"
  threads: 1
  shell:
    "mlst {input.assembly} | process_mlst.pl {sample}_{pipeline_version} > {output.mlst};"

dirCentrifuge = config["Outputs"]["centrifuge_dir"]
if not os.path.exists(dirCentrifuge):
    os.makedirs(dirCentrifuge)

rule centrifuge:
  input:
    assembly = rules.Assembly.output.fasta,
    species_verification = rules.Verify_Species.output.verification,
    scaffolds_table = rules.Assembly.output.table2
  output:
    out = dirCentrifuge + sample + "." + pipeline_version + ".out",
    report = dirCentrifuge + sample + "." + pipeline_version + ".report",
    table = dirCentrifuge + sample + "." + pipeline_version + ".tbl"
  params:
    centrifuge_db = config["Centrifuge"]["db"],
    centrifuge_k = config["Centrifuge"]["k"]
  log:
    logs_dir + str(date) + ".centrifuge.out"
  threads: config["Centrifuge"]["centrifugeCores"]
  shell:
    "centrifuge -k {params.centrifuge_k} -x {params.centrifuge_db} -f {input.assembly} -p {threads} -S {output.out} --report-file {output.report};"
    "echo \"scaffold\tcircular\tcentrifuge_sp\tseqID\" > {output.table};"
    "combine_centrifuge.pl {output.out} {output.report} {input.species_verification} {input.scaffolds_table} | cut -f 1,2,3,4 >> {output.table};"

rule Annotation:
  input:
    assembly = rules.Assembly.output.fasta
  output:
    rgi_out = config["Outputs"]["rgi_outdir"] + "rgi.txt",
    prokka_tbl = config["Outputs"]["annotation_dir"] + sample + pipeline_version +  ".annot.tbl",
    rgi_tbl = config["Outputs"]["annotation_dir"] + sample + pipeline_version + ".rgi.tbl",
    annotation_4jb = config["Outputs"]["annotation_dir"] + sample + pipeline_version + ".annot.4jb.gff3"
  params:
    prokka_out = config["Outputs"]["prokka_outdir"],
    genus = species[0],
    species = species[1],
    prefix = sample + pipeline_version
  log: 
    logs_dir + str(date) + ".annotation.out"
  threads: config["Annotation"]["annotCores"]
  shell:
    "prokka --addgenes --addmrna --cpus {threads} --norrna --outdir {params.prokka_out} --kingdom Bacteria --genus {params.genus} --species {params.species} {input.assembly};"
    "cd `dirname {output.rgi_out}`;"
    "source ~/init_shell.sh;"
    "conda deactivate;"
    "conda activate {scripts_dir}/../conda_environment_CRE3;"
    "echo \"Running RGI. Database version:\";"
    "rgi database -v;"
    "rgi main -n {threads} -t contig -i {input.assembly} -o rgi;"
    "cd `dirname {output.annotation_4jb}`;"
    "ln -s {params.prokka_out}/*.gff prokka.gff3;"
    "ln -s {params.prokka_out}/*.faa prokka.faa;"
    "join_annotations.pl -prokka prokka.gff3 -pep prokka.faa -rgi {output.rgi_out} -pref {params.prefix};"

additional_amrf_opts = None

if species[0] == "Escherichia" or species[0] == "Shigella":
  additional_amrf_opts = "-O Escherichia --plus "
elif species[0] == "Klebsiella":
  additional_amrf_opts = "-O Klebsiella --plus "
elif species[0] == "Salmonella":
  additional_amrf_opts = "-O Salmonella --plus "
elif species[0] == "Staphylococcus":
  additional_amrf_opts = "-O Staphylococcus --plus "
elif species[0] == "Campylobacter":
  additional_amrf_opts = "-O Campylobacter"
elif config["Parameters"]["LIMS_species"] == "Enterococcus faecalis":
  additional_amrf_opts = "-O Enterococcus_faecalis"
elif config["Parameters"]["LIMS_species"] == "Enterococcus faecium":
  additional_amrf_opts = "-O Enterococcus_faecium"
elif config["Parameters"]["LIMS_species"] == "Vibrio cholerae":
  additional_amrf_opts = "-O Vibrio_cholerae --plus "

if additional_amrf_opts == None:
    additional_amrf_opts = ""

rule AMRFinder:
  input:
    annotation = rules.Annotation.output.prokka_tbl,
    annotation_gff = rules.Annotation.output.annotation_4jb
  output:
    proteins = config["Outputs"]["annotation_dir"] + sample + pipeline_version + ".pep.fa",
    amrfinder =  config["Outputs"]["annotation_dir"] + sample + pipeline_version + ".amrfinder.tsv",
    amrfinder_gff =  config["Outputs"]["annotation_dir"] + sample + pipeline_version + ".amrfinder.4jb.gff3",
  params:
  log: 
    logs_dir + str(date) + ".amrf.out"
  threads: 1
  shell:
    "cut -f 5,12 {input.annotation} | grep -v protein_sequence |  TblToFasta > {output.proteins};"
    "amrfinder {additional_amrf_opts} --threads {threads} -p {output.proteins} | change_header_amrfinder.py > {output.amrfinder};"
    "parseAMRFinder4jb.pl -annot {input.annotation_gff} -amrf {output.amrfinder} > {output.amrfinder_gff};"

rule resfinder:
  input:
    assembly = rules.Assembly.output.fasta,
    db = config["ResFinder"]["db"]
  output:
    gff = config["Outputs"]["resfinder_gff"],
    tbl = config["ResFinder"]["dir"] + "ResFinder_results_tab.txt",
    pheno = config["ResFinder"]["dir"] + "pheno_table.txt"
  params:
    resfinder_dir = config["ResFinder"]["dir"],
    path = config["ResFinder"]["path"],
    score = config["ResFinder"]["score"]
  log:
    logs_dir + str(date) + ".resfinder.out"
  threads: 1
  shell: 
    "mkdir -p {params.resfinder_dir};"
    "export PATH={params.path}:$PATH;"
    "run_resfinder.py -ifa {input.assembly} -o {params.resfinder_dir} -db_res {input.db} -t {params.score} -acq;"
    "parse_resfinder_output.pl {params.resfinder_dir}/resfinder_blast/tmp > {output.gff};"

rule IS_Finder:
  input:
    assembly = rules.Assembly.output.fasta,
    db = config["ISFinder"]["IS_db"]
  output:
    ISF_out = config["Outputs"]["ISFinder_out"] 
  params:
    IS_dir = os.path.dirname(config["Outputs"]["ISFinder_out"])
  log: 
    logs_dir + str(date) + ".ISFinder.out"
  threads: config["ISFinder"]["IS_cores"]
  shell:
    "cd $TMPDIR;"
    "ln -s {input.assembly} assembly.fa;"
    "makeblastdb -in assembly.fa -dbtype nucl -parse_seqids;"
    "blastn -db assembly.fa -query {input.db} -out {params.IS_dir}/ISFinder.out -outfmt 7 -num_threads {threads};"
    "cd {params.IS_dir};"
    "grep -v \'#\' ISFinder.out | sort -k2,2 | perl -ane \'chomp; $c++; if ($F[8]<$F[9]){{print \"$F[1]\\tISFinder\\tIS_element\\t$F[8]\\t$F[9]\\t$F[11]\\t+\\t.\\tID=IS.$c;Target=$F[0]:$F[6]-$F[7];Name=$F[0];Evalue=$F[10]\\n\";}}else {{print \"$F[1]\\tISFinder\\tIS_element\\t$F[9]\\t$F[8]\\t$F[11]\\t-\\t.\\tID=IS.$c;Target=$F[0]:$F[6]-$F[7];Name=$F[0];Evalue=$F[10]\\n\";}}\' > {output.ISF_out};"
