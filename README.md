# AACRE
The AACRE pipeline (Assembly and Annotation of Carbapenem-Resistant Enterobacteriaceae) is a Snakemake-based workflow that was developed to automatically process, assemble and annotate Carbapenem-Resistant enterobactericiae that have been sequenced with Illumina and ONT (Oxford Nanopore Technologies). It takes as input fastq files and takes care of pre-processing the reads, performing the hybrid assembly and annotating the assembled genomes with special attention on genes involved in conferring resistance to antiobiotics. 
All the needed software has been installed under conda environments. To install:

1- Clone this directory

2- Install the three conda environments:
```bash
  conda env create --file conda_environment_CRE.yml
  conda env create --file conda_environment_CRE2.yml
  conda env create --file conda_environment_CRE3.yml
```

3- The Resfinder software needs to be downloaded from (https://cge.cbs.dtu.dk/services/ResFinder/). Please, download and install the software and database following the developer's instructions.
Place the resfinder directory in the base 'AACRE' directory and the database in AACRE/databases/resfinder_db or specify the correct path when running CRE_get_config.py (options --resfinder-db and --resfinder-path, do CRE_get_config.py -h for more details).

4- If you are using a HPC cluster and want to specify resources for each job, modify the AACRE_pipeline.spec accordingly or adapt it to your situation (do not specify number of threads in the .spec file, this can be done directly in the configuration file of the run).

5- To launch the pipeline, you need to first activate the "conda_environment_CRE". Only two steps use any of the other environments and the pipeline itself loads the correct one during the process. Here a example of a basic run:
```bash
  conda activate PATH/AACRE/conda_environment_CRE
  CRE_get_config.py --sample AH0325 --configFile AH0325.config.json --illumina-reads directory_with_illumina_reads --ont-reads directory_with_ont_reads --sp "Klebsiella pneumoniae"
  ##Run it in a single machine:
  snakemake --snakefile AACRE.smk  --configfile AH0325.config.json --is -np  ##np means "dry-run" it's useful to first check the commands that will be run, remove the -np when you're sure that you want to run them

  ##Example of a run in a lustre cluster:
snakemake  --notemp -j 999 --snakefile AACRE.smk --cluster "python3 sbatch.py {dependencies}" --configfile AH0325.config.json  --cluster-config AACRE_pipeline.spec --is -np
```
