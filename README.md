# AACRE
The AACRE pipeline (Assembly and Annotation of Carbapenem-Resistant Enterobacteriaceae) is a Snakemake-based workflow that was developed to automatically process, assemble and annotate Carbapenem-Resistant enterobactericiae that have been sequenced with Illumina and ONT (Oxford Nanopore Technologies). It takes as input fastq files and takes care of pre-processing the reads, performing the hybrid assembly and annotating the assembled genomes with special attention on genes involved in conferring resistance to antiobiotics. 
All the needed software has been installed under conda environments. To install:

1- Clone this directory

2- Install the non-provided databases:

-- Bacteria centrifuge database (https://ccb.jhu.edu/software/centrifuge/manual.shtml#custom-database),
-- Resfinder database with the following command from the location where you want to place it:
```bash
  git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git db_resfinder
```

3- If you are using a HPC cluster and want to specify resources for each job, modify the AACRE_pipeline.spec accordingly or adapt it to your situation (do not specify number of threads in the .spec file, this can be done directly in the configuration file of the run).

4- Before running the snakemake pipeline you need to create a config file by using the bin/CRE_get_config.py script. Afterwards, you can launch the pipeline specifying the path to the configuration file. The conda environments will be installed the first time you launch the pipeline. To launch the pipeline you need to have snakemake in your path (it has been tested with snakemake version 6.3.0). An example of a basic run:
```bash
  bin/CRE_get_config.py --sample AH0325 --configFile AH0325.config.json --illumina-reads directory_with_illumina_reads --ont-reads directory_with_ont_reads --sp "Klebsiella pneumoniae" --resfinder-db path/to/resfinder_db
  ##Run it in a single machine:
  snakemake --snakefile bin/AACRE.smk  --configfile AH0325.config.json --is --use-conda --use-envmodules  -np  ##np means "dry-run" it's useful to first check the commands that will be run, remove the -np when you're sure that you want to run them

  ##Example of a run in a lustre cluster:
snakemake  --notemp -j 999 --snakefile AACRE.smk --cluster "python3 sbatch.py {dependencies}" --configfile AH0325.config.json  --cluster-config AACRE_pipeline.spec --use-conda --use-envmodules --is -np
```

5- The first time you run the pipeline, after the AMRFinder conda environment has been created, you will need to download the amrfinder database. With the amrfinder loaded, please do:
```bash
armfinder -u 
```

6- After a successful run has been completed, you can create a report with details on the process and the results. Please do:
```bash
snakemake --snakefile bin/AACRE.smk  --configfile AH0325.config.json --is --use-conda --use-envmodules --report myreport.zip
```
