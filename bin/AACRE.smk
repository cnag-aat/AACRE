#Authors: Jessica Gomez-Garrido, CNAG-CRG; Fernando Cruz, CNAG-CRG.
#Contact email: jessica.gomez@cnag.crg.eu
#Date:20190704
#Snakemake pipeline to assemble and annotate bacterial genomes from ONT and Illumina genomic data

from datetime import datetime

date = datetime.now().strftime('%Y%m%d.%H%M%S')

scripts_dir = config["Inputs"]["scripts_dir"]
logs_dir = config["Parameters"]["logs_dir"]
if not os.path.exists(logs_dir):
    os.makedirs(logs_dir)
base_dir= config["Parameters"]["basedir"]

# file wildcards 
illumina_files = config["Wildcards"]["ILLUMINA_fastqs"].split(',')
if illumina_files == None:
  illumina_files = ""
  print ("WARNING no Illumina fastqs have been given")

ont_files = config["Wildcards"]["ONT_fastqs"].split(',')
if ont_files == None:
  ont_files = ""
  print ("WARNING no ONT fastqs have been given")

#Global variables definition
sample = config["Parameters"]["sample"]
pipeline_version = "v" + str(config["Parameters"]["version"])
species= config["Parameters"]["metadata_species"].split(" ")
databases = scripts_dir + "../databases/"

fastqs={}
extensions = ["1.fastq.gz", "2.fastq.gz"]
for i in extensions:
  fastqs["illumina." + i] = []
  for file in illumina_files:
    fastqs["illumina." + i].append(config["Outputs"]["ILLUMINA_trim"] + file + ".trimmed." + i)
  if not os.path.exists(config["Outputs"]["ILLUMINA_cat"] + "/logs"):
    os.makedirs(config["Outputs"]["ILLUMINA_cat"] + "/logs")

extensions = ["fastq.gz"]
for i in extensions:
  fastqs["ont."+i] = []
  for file in ont_files:
    fastqs["ont."+i].append(config["Outputs"]["ONT_trim"] +  file + ".trim." + i)
  if not os.path.exists(config["Outputs"]["ONT_cat"] + "/logs"):
    os.makedirs(config["Outputs"]["ONT_cat"] + "/logs")

kraken_ins = {}
kraken_ins[sample + ".ont"] = config["Outputs"]["ONT_cat"] + "reads.ont.fastq.gz"
kraken_ins[sample + ".illumina"] = [config["Outputs"]["ILLUMINA_cat"] + "reads.illumina.1.fastq.gz", config["Outputs"]["ILLUMINA_cat"] + "reads.illumina.2.fastq.gz"]
  
#Modules
module preprocess_workflow:
  snakefile: "../modules/process_reads.rules.smk"
module assembly_workflow:
  snakefile: "../modules/assembly.rules.smk"
module annotation_workflow:
  snakefile: "../modules/annotation.rules.smk"

#Rule for outputs
rule all:
  input:
    config["Outputs"]["ILLUMINA_cat"] +  "reads.illumina.1.fastq.gz",
    config["Outputs"]["ILLUMINA_cat"] + "reads.illumina.2.fastq.gz",
    config["Outputs"]["ONT_cat"] + "reads.ont.fastq.gz",
    config["Outputs"]["assembly"],
    config["Outputs"]["assembly_dir"] + "tables/" + sample + "." + pipeline_version +  ".mlst.tbl",
    expand(config["Outputs"]["nanostats_dir"] + "NanoStats.txt"),
    config["Outputs"]["bracken"] + sample + ".illumina.kraken2.report.txt",
    config["Outputs"]["bracken"] + sample + ".ont.kraken2.report.txt",
    config["Outputs"]["assembly_dir"] + "plots/" + sample + "." + pipeline_version + ".assembly.png",
    config["Outputs"]["preprocessing_dir"] + sample + ".species_verification.out",
    config["Outputs"]["annotation_dir"] + "prokka.gff3",
    config["Outputs"]["annotation_dir"] + sample + pipeline_version + ".annot.4jb.gff3",
    config["Outputs"]["annotation_dir"] + sample + pipeline_version + ".rgi.tbl",
    config["Outputs"]["annotation_dir"] + sample + pipeline_version + ".amrfinder.4jb.gff3",
    config["Outputs"]["centrifuge_dir"] + sample + "." + pipeline_version + ".out",
    config["Outputs"]["resfinder_gff"],
    config["Outputs"]["ISFinder_out"] 
  log:
    logs_dir + str(date) + ".rule_all.j%j.out",
    logs_dir + str(date) + ".rule_all.j%j.err",

#Additional options

dirCentrifuge = config["Outputs"]["centrifuge_dir"]
if not os.path.exists(dirCentrifuge):
    os.makedirs(dirCentrifuge)

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

#Pipeline rules
use rule trim_galore from preprocess_workflow with:
  input:
    read1 = config["Inputs"]["ILLUMINA_reads"] + "{file}.1.fastq.gz",
    read2 = config["Inputs"]["ILLUMINA_reads"] + "{file}.2.fastq.gz"
  output:
    trim1 = config["Outputs"]["ILLUMINA_trim"] + "{file}.trimmed.1.fastq.gz",
    trim2 = config["Outputs"]["ILLUMINA_trim"] + "{file}.trimmed.2.fastq.gz",
    report1 = report(config["Outputs"]["ILLUMINA_trim"] + "{file}.1.fastq.gz_trimming_report.txt",
                caption = "../report/trimgalore1.rst",
                category = "Process reads",
                subcategory = "Illumina"),
    report2 = report(config["Outputs"]["ILLUMINA_trim"] + "{file}.2.fastq.gz_trimming_report.txt",
                caption = "../report/trimgalore2.rst",
                category = "Process reads",
                subcategory = "Illumina")
  params:
    outdir = directory(config["Outputs"]["ILLUMINA_trim"]),
    opts = config["Trim_Galore"]["options"]
  log: 
    logs_dir + str(date) + ".{file}.j%j.trim_galore.out",
    logs_dir + str(date) + ".{file}.j%j.trim_galore.err",
  benchmark:
    logs_dir + str(date) + ".{file}.trim_galore.benchmark.txt",
  conda:
    "../envs/trim_galore0.6.7.yaml"
  threads: config["Trim_Galore"]["Trim_Illumina_cores"]

use rule Trim_ONT_qcat from preprocess_workflow with:
  input:
    read = config["Inputs"]["ONT_reads"] + "{ontfile}.fastq.gz"
  output:
    trim = config["Outputs"]["ONT_trim"] + "{ontfile}.trim.fastq.gz"   
  params:
    options = lambda wildcards: "-b " + config["Outputs"]["ONT_trim"] + wildcards.ontfile + "/tmp/barcode" + config["ONT_Barcodes"][wildcards.ontfile] + ".fastq" if config["ONT_Barcodes"][wildcards.ontfile] != "" \
              else  " -o " + config["Outputs"]["ONT_trim"] + wildcards.ontfile + "/tmp/trimmed.fastq",
    workdir = config["Outputs"]["ONT_trim"] + "{ontfile}/",
    tempdir = config["Outputs"]["ONT_trim"] + "{ontfile}/tmp",
    tempfile = lambda wildcards:  config["Outputs"]["ONT_trim"] + wildcards.ontfile + "/tmp/barcode" + config["ONT_Barcodes"][wildcards.ontfile] + ".fastq" if config["ONT_Barcodes"][wildcards.ontfile] != "" \
              else config["Outputs"]["ONT_trim"] + wildcards.ontfile + "/tmp/trimmed.fastq",
  log:
    logs_dir + str(date) + ".{ontfile}.QCAT.j%j.out",
    logs_dir + str(date) + ".{ontfile}.QCAT.j%j.err",   
  benchmark:
    logs_dir + str(date) + ".{ontfile}.QCAT.benchmark.txt",
  conda:
    "../envs/qcat1.1.0.yaml"
  threads: config["Parameters"]["QCAT_cores"]

use rule concat_reads from preprocess_workflow with:
  input:
    fastqs = lambda wildcards: fastqs[wildcards.ext]
  output:
    final_fastq = "{dir}reads.{ext}"
  log:
    "{dir}logs/" + str(date) + ".j%j.concat.{ext}.out",
    "{dir}logs/" + str(date) + ".j%j.concat.{ext}.err"
  benchmark:
    "{dir}logs/" + str(date) + ".concat.benchmark.{ext}.txt"
  conda:
    '../envs/ass_base.yaml'
  threads: config["Parameters"]["Concat_cores"] 

use rule Kraken2 from preprocess_workflow with:
  input:
    read = lambda wildcards: kraken_ins[wildcards.base],
    database = config["Bracken"]["kraken_db"],
    kmers = config["Bracken"]["kraken_kmers"]
  output:
    report = report (config["Outputs"]["bracken"] + "{base}.kraken2.report.txt",
             caption="../report/kraken.rst",
             category = "Process reads",
             subcategory = "Kraken reports"),
    abundance =  config["Outputs"]["bracken"] + "{base}.bracken_abundance.txt",
    tophits =  config["Outputs"]["bracken"] + "{base}.bracken_abundance.tophits.txt",
    readsout = config["Outputs"]["bracken"] + "{base}.kraken2.seqs.out",
  params:
    additional = lambda wildcards: config["Bracken"]["additional_opts"] + " --paired " \
                 if re.search("illumina", wildcards.base) \ 
                 else config["Bracken"]["additional_opts"],
    prefix = lambda wildcards: os.path.basename(wildcards.base),
    scripts_dir = scripts_dir
  log:
    logs_dir + str(date) + ".j%j.{base}.bracken.out",
    logs_dir + str(date) + ".j%j.{base}.bracken.err",
  benchmark:
    logs_dir + str(date) + ".{base}.bracken.benchmark.txt",
  conda:
    '../envs/kraken2.1.2.yaml'
  threads: config["Bracken"]["brackenCores"]

use rule nanoplot from preprocess_workflow with:
  input:
    fastq = config["Outputs"]["ONT_cat"] + "reads.ont.fastq.gz",
  output:
    stats = report(config["Outputs"]["nanostats_dir"] + "NanoStats.txt",
            caption="../report/nanostats.rst",
            category = "Process reads"),
    yield_len = report(config["Outputs"]["nanostats_dir"] + "Yield_By_Length.png",
                caption="../report/nanostats.rst",
                category = "Process reads"),
    read_len = report(config["Outputs"]["nanostats_dir"] + "WeightedHistogramReadlength.png",
               caption= "../report/nanostats.rst",
               category = "Process reads"),
  params:
    outdir = config["Outputs"]["nanostats_dir"]
  log:
    logs_dir + str(date) + ".j%j.NanoPlot.out",
    logs_dir + str(date) + ".j%j.NanoPlot.err",
  benchmark:
    logs_dir + str(date) + ".NanoPlot.benchmark.txt",
  conda:
    '../envs/nanoplot1.40.0.yaml'
  threads: config["Parameters"]["nanostat_cores"]

use rule Verify_Species from preprocess_workflow with:
  input: 
    illumina_tophit = config["Outputs"]["bracken"] + sample + ".illumina.bracken_abundance.tophits.txt",
    nanopore_tophit = config["Outputs"]["bracken"] + sample + ".ont.bracken_abundance.tophits.txt",
  output:
    verification = report(config["Outputs"]["preprocessing_dir"] + sample + ".species_verification.out",
                   caption="../report/verify_sp.rst",
                   category = "Process reads")
  params:
    workdir = config["Outputs"]["preprocessing_dir"],
    metadata_species = "\"" + config["Parameters"]["metadata_species"] + "\"",
    scripts_dir = scripts_dir
  log:
    logs_dir + str(date) + ".j%j.verify_species.out",
    logs_dir + str(date) + ".j%j.verify_species.err",
  benchmark:
    logs_dir + str(date) + ".verify_species.benchmark.txt",
  threads: 1

use rule unicycler from assembly_workflow with:
  input:
    read1 = config["Outputs"]["ILLUMINA_cat"] +  "reads.illumina.1.fastq.gz",
    read2 = config["Outputs"]["ILLUMINA_cat"] + "reads.illumina.2.fastq.gz",
    nanopore = config["Outputs"]["ONT_cat"] + "reads.ont.fastq.gz",
  output:
     fasta = config["Outputs"]["assembly_dir"] + "unicycler_out/assembly.fasta",
     plot = report(config["Outputs"]["assembly_dir"] + "plots/" + sample + "." + pipeline_version + ".assembly.png",
                  caption="../report/assembly.rst",
                  category = "Assembly")
  params:
    assemblydir = config["Outputs"]["assembly_dir"],
    start_genes =  databases + "/unicycler_gene_data/start_genes.fasta",
  log:
    logs_dir + str(date) + ".j%j.assembly.out",
    logs_dir + str(date) + ".j%j.assembly.err",
  benchmark:
    logs_dir + str(date) + ".assembly.benchmark.txt"
  conda:
    "../envs/unicycler0.5.0.yaml"
  threads: config["Assembly"]["assemblyCores"]

use rule produce_final_assembly from assembly_workflow with:
  input:
    fasta = config["Outputs"]["assembly_dir"] + "unicycler_out/assembly.fasta",
    plot = config["Outputs"]["assembly_dir"] + "plots/" + sample + "." + pipeline_version + ".assembly.png",
  output:
    fasta = config["Outputs"]["assembly"],
    table1 = report(config["Outputs"]["assembly_dir"] + "tables/" + sample + "." + pipeline_version +  ".assembly.tbl",
                    caption="../report/assembly.rst",
                    category = "Assembly"),
    table2 = report (config["Outputs"]["assembly_dir"] + "tables/" + sample + "." + pipeline_version +  ".scaffolds.tbl",
                    caption="../report/assembly.rst",
                    category = "Assembly"),
    nseries = config["Outputs"]["assembly_dir"] + "tables/" + sample + "." + pipeline_version +  ".assembly.nseries.txt"
  params:
    assemblydir = config["Outputs"]["assembly_dir"],
    scripts_dir = scripts_dir,
    sample = sample,
    pipeline_version = pipeline_version
  log:
    logs_dir + str(date) + ".j%j.finalize_assembly.out",
    logs_dir + str(date) + ".j%j.finalize_assembly.err",
  benchmark:
    logs_dir + str(date) + ".finalize_assembly.benchmark.txt"
  conda:
    "../envs/ass_base.yaml"

use rule MLST from annotation_workflow with:
  input:
    assembly = config["Outputs"]["assembly"],
  output: 
    mlst = report (config["Outputs"]["assembly_dir"] + "tables/" + sample + "." + pipeline_version +  ".mlst.tbl",
                  caption="../report/mlst.rst",
                  category = "Assembly")
  params:
    scripts_dir = scripts_dir,
    sample = sample,
    pipeline_version = pipeline_version
  log:
    logs_dir + str(date) + ".j%j.mlst.out",
    logs_dir + str(date) + ".j%j.mlst.err"
  benchmark:
    logs_dir + str(date) + ".mlst.benchmark.txt"
  conda:
    "../envs/mlst2.22.1.yaml"
  threads: 1

use rule centrifuge from annotation_workflow with:
  input:
    assembly =  config["Outputs"]["assembly"],
    species_verification = config["Outputs"]["preprocessing_dir"] + sample + ".species_verification.out",
    scaffolds_table = config["Outputs"]["assembly_dir"] + "tables/" + sample + "." + pipeline_version +  ".scaffolds.tbl",
  output:
    out = dirCentrifuge + sample + "." + pipeline_version + ".out",
    report = report (dirCentrifuge + sample + "." + pipeline_version + ".report",
                    caption="../report/centrifuge.rst",
                    category = "Assembly"),
    table = dirCentrifuge + sample + "." + pipeline_version + ".tbl"
  params:
    centrifuge_db = config["Centrifuge"]["db"],
    centrifuge_k = config["Centrifuge"]["k"],
    scripts_dir = scripts_dir
  log:
    logs_dir + str(date) + ".j%j.centrifuge.out",
    logs_dir + str(date) + ".j%j.centrifuge.err"
  benchmark:
    logs_dir + str(date) + ".centrifuge.benchmark.txt"
  conda:
    "../envs/centrifuge1.0.4.yaml"
  threads: config["Centrifuge"]["centrifugeCores"]

use rule Prokka from annotation_workflow with:
  input:
    assembly = config["Outputs"]["assembly"],
  output:
    gff = config["Outputs"]["annotation_dir"] + "prokka.gff3",
    faa = config["Outputs"]["annotation_dir"] + "prokka.faa"
  params:
    annotation_dir = config["Outputs"]["annotation_dir"],
    prokka_out = config["Outputs"]["prokka_outdir"],
    genus = species[0],
    species = species[1],
    prefix = sample + pipeline_version
  log:
    logs_dir + str(date) + ".j%j.prokka.out",
    logs_dir + str(date) + ".j%j.prokka.err"
  benchmark:
    logs_dir + str(date) + ".prokka.benchmark.txt"
  conda:
    "../envs/prokka1.14.6.yaml"
  threads: config["Annotation"]["annotCores"]

use rule RGI from annotation_workflow with:
  input:
    assembly = config["Outputs"]["assembly"],
    prokka_gff = config["Outputs"]["annotation_dir"] + "prokka.gff3",
    prokka_faa = config["Outputs"]["annotation_dir"] + "prokka.faa"
  output:
    rgi_out = config["Outputs"]["rgi_outdir"] + "rgi.txt",
    prokka_tbl = config["Outputs"]["annotation_dir"] + sample + pipeline_version +  ".annot.tbl",
    rgi_tbl = report (config["Outputs"]["annotation_dir"] + sample + pipeline_version + ".rgi.tbl",
                     caption="../report/annotation.rst",
                     category = "Annotation"),
    annotation_4jb = config["Outputs"]["annotation_dir"] + sample + pipeline_version + ".annot.4jb.gff3",
    proteins = config["Outputs"]["annotation_dir"] + sample + pipeline_version + ".pep.fa",
  params:
    annotation_dir = config["Outputs"]["annotation_dir"],
    scripts_dir = scripts_dir,
    prefix = sample + pipeline_version
  log:
    logs_dir + str(date) + ".j%j.rgi.out",
    logs_dir + str(date) + ".j%j.rgi.err"
  benchmark:
    logs_dir + str(date) + ".rgi.benchmark.txt"
  conda:
    "../envs/RGI6.6.0.yaml"
  threads: config["Annotation"]["annotCores"]

use rule AMRFinder from annotation_workflow with:
  input:
    annotation = config["Outputs"]["annotation_dir"] + sample + pipeline_version +  ".annot.tbl",
    annotation_gff = config["Outputs"]["annotation_dir"] + sample + pipeline_version + ".annot.4jb.gff3",
    proteins = config["Outputs"]["annotation_dir"] + sample + pipeline_version + ".pep.fa",
  output:
    amrfinder =  report (config["Outputs"]["annotation_dir"] + sample + pipeline_version + ".amrfinder.tsv",
                        caption="../report/amrfinder.rst",
                        category = "Annotation"),
    amrfinder_gff =  config["Outputs"]["annotation_dir"] + sample + pipeline_version + ".amrfinder.4jb.gff3",
  params:
    scripts_dir = scripts_dir,
    additional = additional_amrf_opts
  log: 
    logs_dir + str(date) + ".j%j.amrfinder.out",
    logs_dir + str(date) + ".j%j.amrfinder.err",
  benchmark:
    logs_dir + str(date) + ".amrfinder.benchmark.txt",
  conda:
    "../envs/amrfinder3.10.42.yaml"    
  threads: 1

use rule resfinder from annotation_workflow with:
  input:
    assembly = config["Outputs"]["assembly"],
    db = config["ResFinder"]["db"]
  output:
    gff = config["Outputs"]["resfinder_gff"],
    tbl = report (config["ResFinder"]["dir"] + "ResFinder_results_tab.txt",
                  caption="../report/annotation.rst",
                  category = "Annotation"),
    pheno = config["ResFinder"]["dir"] + "pheno_table.txt"
  params:
    resfinder_dir = config["ResFinder"]["dir"],
    score = config["ResFinder"]["score"],
    scripts_dir = scripts_dir
  log:
    logs_dir + str(date) + ".j%j.resfinder.out",
    logs_dir + str(date) + ".j%j.resfinder.err"
  benchmark:
    logs_dir + str(date) + ".resfinder.benchmark.txt"
  conda:
    "../envs/resfinder4.1.11.yaml"
  threads: 1

use rule IS_Finder from annotation_workflow with:
  input:
    assembly = config["Outputs"]["assembly"],
    db = config["ISFinder"]["IS_db"]
  output:
    ISF_out = config["Outputs"]["ISFinder_out"],
  params:
    IS_dir = os.path.dirname(config["Outputs"]["ISFinder_out"])
  log: 
    logs_dir + str(date) + ".j%j.ISFinder.out",
    logs_dir + str(date) + ".j%j.ISFinder.err"
  benchmark:
    logs_dir + str(date) + ".ISFinder.benchmark.txt"
  conda:
    "../envs/BLAST2.13.0.yaml"
  threads: config["ISFinder"]["IS_cores"]
