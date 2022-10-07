#!/usr/bin/env python3
import os
import json
import argparse
import sys
import re

#Author: Jessica Gomez-Garrido, CNAG-CRG.
#Contact email: jessica.gomez@cnag.crg.eu
#Date:20190704

ont_barcodes = []
class CreateConfigurationFile(object):
    """Class which manages Configuration file Manager"""
      
    def __init__(self):
        """Class constructor"""
        #GENERAL PARAMETERS
        self.version = 2                                                                          #Pipeline run version
        self.sample = None                                                                        #Sample barcode 
        self.configFile = None                                                                    #Name of the json configuration file to be created.
        self.basedir = None                                                                       #Base directory for the pipeline run
        self.logs_dir = "logs/"                                                                   #Directory to keep all the log files
        self.metadata_species = None                                                              #Species listed in the metadata
        self.Qcat_cores = 4                                                                       #Number of threads to run the QCAT step
        self.Concat_cores = 4                                                                     #Number of threads to run the concat reads step
        self.nanostat_cores = 2                                                                   #Number of threads to tun the nanostats step

        #INPUT PARAMETERS
        self.scripts_dir = os.path.dirname(sys.argv[0]) + "/../scripts/"                          #Directory with the different scripts for the pipeline
        self.ILLUMINA_reads = None                                                                #Directory where the illumina fastqs are stored
        self.ONT_reads = None                                                                     #Directory where the ont fastqs are stored     

        #OUTPUT PARAMETERS
        self.preprocessing_dir =  "01_preprocessing/"                                             #Base directory of the preprocessing step
        self.ILLUMINA_trim = self.preprocessing_dir + "illumina/trim/"                            #Base directory of the illumina trimming step
        self.ILLUMINA_cat = self.preprocessing_dir + "illumina/out/"                              #Out directory of the illumina trimming step
        self.ILLUMINA_qc = self.preprocessing_dir + "illumina/qc/"                                #Base directory of the illumina qc step
        self.bracken = self.preprocessing_dir + "bracken/"                                        #Directory for the Bracken step of the reads
        self.ONT_trim = self.preprocessing_dir + "ont/trim/"                                      #Base directory of the ont trimming step
        self.ONT_cat = self.preprocessing_dir + "ont/out/"                                        #Out directory of the ont trimming step
        self.ONT_qc =  self.preprocessing_dir + "ont/qc/"                                         #Base directory of the ont qc step
        self.nanostats_dir = self.preprocessing_dir + "ont/nanostats/"                            #Base directory for the nanostats step
        self.assembly_dir = "02_Assembly/"                                                        #Base directory of the assembly step
        self.annotation_dir = "03_Annotation/"                                                    #Base directory of the annotation step
        self.prokka_outdir = self.annotation_dir + "prokka/"                                      #Directory to store the prokka outputs
        self.rgi_outdir = self.annotation_dir + "RGI/"                                            #Directory to store the RGI outputs
        self.assembly = None                                                                      #Assembly file 
        self.centrifuge_dir = self.assembly_dir + "centrifuge/"                                   #Directory to store Centrifuge output
        self.ISFinder = None                                                                      #File to keep the ISFinder results
        self.resfinder = None                                                                     #File to keep the resfinder results         

        #WILDCARDS
        self.ILLUMINA_fastqs = None                                                               #List with basename of the illumina fastqs
        self.ONT_fastqs = None                                                                    #List with basename of the ONT fastqs

        #TRIMGALORE PARAMETERS
        self.trim_galore_opts = "--gzip -q 10 --paired --retain_unpaired"
        self.Trim_Illumina_cores = 4                                                              #Number of threads to run the trim Illumina step

        #BRACKEN PARAMETERS
        self.brackenCores = 4                                                                     #Number of threads to run the bracken step 
        self.brackenDB = os.path.dirname(sys.argv[0]) + "/../databases/minikraken2_v1_8GB"           #Database to run Bracken
        self.brackenKmers = self.brackenDB + "/database200mers.kmer_distrib"                      #Kmer distribution to run Bracken
        self.additional_kraken2_opts = ""

        #ASSEMBLY PARAMETERS
        self.assemblyCores = 16                                                                    #Number of threads needed to run the assembly step

        #CENTRIFUGE PARAMETERS
        self.centrifugeCores = 12                                                                  #Number of threads needed to run the centrifuge step
        self.centrifuge_db ="/scratch/devel/jgomez/centrifuge_cre_database/bacteria/bact"         #Database to run centrifuge, CHANGE THIS LINE accordingly to the location of your db
        self.centrifuge_k = 5                                                                     #Centrifuge: report upto <int> distinct, primary assignments for each read or pair

        #ANNOTATION PARAMETERS
        self.annotCores = 4                                                                       #Default number of threads for the annotation step

        #ISFINDER PARAMETERS
        self.IS_db = os.path.dirname(sys.argv[0])+ "/../databases/ISFinder/ISfinder_genes.fasta"      #Database to run ISFinder
        self.IS_cores = 8                                                                         #Default number of threads for the annotation step

        #RESFINDER PARAMETERS
        self.resfinder_db = os.path.dirname(sys.argv[0])+ "/../databases/resfinder_db/"              #Location of the database to run resfinder
        self.resfinder_score = 0.8
        self.resfinder_dir = None                                                                 #Directory to run resfinder
###
        #DICTIONARIES
        self.allParameters = {}
        self.generalParameters = {}
        self.inputParameters = {}
        self.outputParameters = {}
        self.wildcardParameters = {}
        self.trimgaloreParameters = {}
        self.qcatParameters = {}
        self.brackenParameters = {}
        self.assemblyParameters = {}
        self.centrifugeParameters = {}
        self.annotationParameters = {}
        self.ISFinderParameters = {}
        self.ResFinderParameters = {}

####

    def register_parameter(self, parser):
        """Register all parameters with the given
        argparse parser"""
        self.register_general(parser)
        self.register_input(parser)
        self.register_output(parser)
        self.register_wildcards(parser)
        self.register_trimgalore(parser)
        self.register_bracken(parser)
        self.register_assembly(parser)
        self.register_annotation(parser)
        self.register_centrifuge(parser)
        self.register_ISFinder(parser)
        self.register_ResFinder(parser)

###

    def register_general(self, parser):
        """Register all general parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        general_group = parser.add_argument_group('General Parameters')
        general_group.add_argument('--version', dest="version", metavar="version", type=int, default=self.version, help='Pipeline run version. Default %s' % self.version)
        general_group.add_argument('--sample', dest="sample", metavar="sample", help='Sample barcode. Default %s' % self.sample)
        general_group.add_argument('--configFile', dest="configFile", metavar="configFile", help='Configuration JSON to be generated. Default %s' % self.configFile)
        general_group.add_argument('--basedir', dest="basedir", metavar="basedir", help='Base directory for the pipeline run. Default %s' % self.basedir)
        general_group.add_argument('--logs-dir', dest="logs_dir", metavar="logs_dir", help='Directory to keep all the log files. Default %s' % self.logs_dir)
        general_group.add_argument('--sp', dest="metadata_species", metavar="metadata_species", help='Species listed in the metadata. Default %s' % self.metadata_species)
        general_group.add_argument('--qcat-cores', type = int, dest="Qcat_cores", metavar="QCAT_cores", default=self.Qcat_cores, help='Number of threads to run the QCAT ONT step. Default %s' % self.Qcat_cores)
        general_group.add_argument('--cat-cores', type = int, dest="Concat_cores", metavar="Concat_cores", default=self.Concat_cores, help='Number of threads to run the concat reads step. Default %s' % self.Concat_cores)
        general_group.add_argument('--nanostat-cores', type = int, dest="nanostat_cores", metavar="nanostat_cores", default=self.nanostat_cores, help='Number of threads to run the nanostat step. Default %s' % self.nanostat_cores)

    def register_input(self, parser):
        """Register all input parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        input_group = parser.add_argument_group('Inputs')
        input_group.add_argument('--scripts-dir', dest="scripts_dir", help='Directory with the different scripts for the pipeline. Default %s' % self.scripts_dir)
        input_group.add_argument('--illumina-reads', dest="ILLUMINA_reads", help='Directory where the illumina fastqs are stored. Default %s' % self.ILLUMINA_reads)
        input_group.add_argument('--ont-reads', dest="ONT_reads", help='Directory where the ont fastqs are stored. Default %s' % self.ONT_reads)

    def register_output(self, parser):
        """Register all output parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        output_group = parser.add_argument_group('Outputs')
        output_group.add_argument('--preprocessing-dir', dest="preprocessing_dir", help='Base directory of the preprocessing step. Default %s' % self.preprocessing_dir)
        output_group.add_argument('--illumina-trim', dest="ILLUMINA_trim", help='Base directory of the illumina trimming step. Default %s' % self.ILLUMINA_trim)
        output_group.add_argument('--illumina-cat', dest="ILLUMINA_cat", help='Out directory of the illumina trimming step. Default %s' % self.ILLUMINA_cat)
        output_group.add_argument('--illumina-qc', dest="ILLUMINA_qc", help='Base directory of the illumina qc step. Default %s' % self.ILLUMINA_qc)
        output_group.add_argument('--bracken-dir', dest="bracken", help='Directory for the Bracken step of the reads. Default %s' % self.bracken)
        output_group.add_argument('--ont-trim', dest="ONT_trim", help='Base directory of the ONT trimming step. Default %s' % self.ONT_trim)
        output_group.add_argument('--ont-cat', dest="ONT_cat", help='Out directory of the ONT trimming step. Default %s' % self.ONT_cat)
        output_group.add_argument('--ont-qc', dest="ONT_qc", help='Base directory of the ONT qc step. Default %s' % self.ONT_qc)
        output_group.add_argument('--nanostats-dir', dest="nanostats_dir", help='Base directory of the nanostats step. Default %s' % self.nanostats_dir)
        output_group.add_argument('--assembly-dir', dest="assembly_dir", help='Base directory of the assembly step. Default %s' % self.assembly_dir)
        output_group.add_argument('-g', dest="assembly", help='Assembly file (It will be given through the pipeline once it is complete). Default %s' % self.assembly)
        output_group.add_argument('--centrifuge-dir', dest="centrifuge_dir", help='Base directory of the centrifuge step. Default %s' % self.centrifuge_dir)
        output_group.add_argument('--annotation-dir', dest="annotation_dir", help='Base directory of the annotation step. Default %s' % self.annotation_dir)
        output_group.add_argument('--prokka-outdir', dest="prokka_outdir", help='Directory to store the prokka outputs. Default %s' % self.prokka_outdir)
        output_group.add_argument('--rgi-outdir', dest="rgi_outdir", help='Directory to store the RGI outputs. Default %s' % self.rgi_outdir)
        output_group.add_argument('--ISFinder', dest="ISFinder", help='File to keep te ISFinder results. Default %s' % self.ISFinder)
        output_group.add_argument('--resfinder', dest="resfinder", help='File to keep te resfinder results. Default %s' %self.resfinder)

    def register_wildcards(self, parser):
        """Register all wildcards parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        wildcards_group = parser.add_argument_group('Wildcards')
        wildcards_group.add_argument('--ILLUMINA-fastqs', dest="ILLUMINA_fastqs", metavar="ILLUMINA_fastqs", help='List with basename of the illumina fastqs. Default %s' % self.ILLUMINA_fastqs)
        wildcards_group.add_argument('--ONT-fastqs', dest="ONT_fastqs", metavar="ONT_fastqs", help='List with basename of the ONT fastqs. Default %s' % self.ONT_fastqs)

    def register_trimgalore(self, parser):
        """Register all trimgalore parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        trimgalore_group = parser.add_argument_group('Trim_Galore')
        trimgalore_group.add_argument('--trim-galore-opts', dest="trim_galore_opts", metavar="trim_galore_opts", default=self.trim_galore_opts, help='Optional parameters for the rule trim_galore. Default %s' % self.trim_galore_opts)
        trimgalore_group.add_argument('--trim-Illumina-cores', type = int, dest="Trim_Illumina_cores", metavar="Trim_Illumina_cores", default=self.Trim_Illumina_cores, help='Number of threads to run the Illumina trimming step. Default %s' % self.Trim_Illumina_cores)

    def register_bracken(self, parser):
        """Register all Bracken parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        bracken_group = parser.add_argument_group('Bracken')
        bracken_group.add_argument('--brackenCores', dest="brackenCores", metavar="brackenCores", type=int, default=self.brackenCores, help='Number of threads to run the bracken step with the  reads . Default %s' % self.brackenCores)
        bracken_group.add_argument('--brackenDB', dest="brackenDB", metavar="brackenDB", default=self.brackenDB, help='Database to run Bracken. Default %s' % self.brackenDB)
        bracken_group.add_argument('--brackenKmers', dest="brackenKmers", metavar="brackenKmers", default=self.brackenKmers, help='Kmer distribution to run Bracken. Default %s' % self.brackenKmers)
        bracken_group.add_argument('--kraken2-opts', dest="additional_kraken2_opts", metavar="additional_kraken2_opts", default=self.additional_kraken2_opts, help='Optional parameters for the rule Kraken2. Default %s' % self.additional_kraken2_opts)
       
    def register_assembly(self, parser):
        """Register all assembly parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        assembly_group = parser.add_argument_group('Assembly')
        assembly_group.add_argument('--assemblyCores', dest="assemblyCores", metavar="assemblyCores", type=int, default=self.assemblyCores, help='Number of threads needed to run the assembly step. Default %s' % self.assemblyCores)

    def register_centrifuge(self, parser):
        """Register all centrifuge parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        centrifuge_group = parser.add_argument_group('Centrifuge')
        centrifuge_group.add_argument('--centrifugeCores', dest="centrifugeCores", metavar="centrifugeCores", type=int, default=self.centrifugeCores, help='Number of threads needed to run the centrifuge step. Default %s' % self.centrifugeCores)
        centrifuge_group.add_argument('--centrifuge_db', dest="centrifuge_db", metavar="centrifuge_db", default=self.centrifuge_db, help=' Default %s' % self.centrifuge_db)
        centrifuge_group.add_argument('--centrifuge_k', dest="centrifuge_k", metavar="centrifuge_k", type=int, default=self.centrifuge_k, help='report upto <int> distinct, primary assignments for each read or pair. Default %s' % self.centrifuge_k)

    def register_annotation(self, parser):
        """Register all annotation parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        annotation_group = parser.add_argument_group('Annotation')
        annotation_group.add_argument('--annotCores', dest="annotCores", metavar="annotCores", type=int, default=self.annotCores, help='Number of threads needed to run the annotation step. Default %s' % self.annotCores)

    def register_ISFinder(self, parser):
        """Register all ISFinder parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        ISFinder_group = parser.add_argument_group('ISFinder')
        ISFinder_group.add_argument('--IS-cores', dest="IS_cores", metavar="IS_cores", type=int, default=self.IS_cores, help='Number of threads needed to run the ISFinder step. Default %s' % self.IS_cores)
        ISFinder_group.add_argument('--IS-db', dest="IS_db", metavar="IS_db", default=self.IS_db, help='Database for ISFinder stepe. Default %s' % self.IS_db)

    def register_ResFinder(self, parser):
        """Register all Resfinder parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        ResFinder_group = parser.add_argument_group('Resfinder')
        ResFinder_group.add_argument('--resfinder-score', dest="resfinder_score", metavar="resfinder_score", type=float, default=self.resfinder_score, help='Number of threads needed to run the Resfinder step. Default %s' % self.resfinder_score)
        ResFinder_group.add_argument('--resfinder-db', dest="resfinder_db", metavar="resfinder_db", default=self.resfinder_db, help='Database for resfinder step. Default %s' % self.resfinder_db)
        ResFinder_group.add_argument('--resfinder-dir', dest="resfinder_dir", metavar="resfinder_dir", default=self.resfinder_dir, help='Directory to run the resfinder step. Default %s' % self.resfinder_dir)

####

    def check_parameters(self,args):
        """Check parameters consistency
            
        args -- set of parsed arguments"""

        working_dir = os.getcwd() + "/"

        if args.sample == None:
            print ("No sample specified")
            parser.print_help()
            sys.exit(-1)

        if args.basedir:
            args.basedir = os.path.abspath(args.basedir) + "/"
        else: 
            args.basedir = working_dir + "v" + str(args.version) + "/" + args.sample + "/"

        if args.configFile != None:
            args.configFile = os.path.abspath(args.configFile) 
        else:
            print ("You need to specify the configfile")
            parser.print_help()
            sys.exit(-1)

        if args.scripts_dir:
            args.scripts_dir = os.path.abspath(args.scripts_dir) + "/"
        else:
            args.scripts_dir = os.path.abspath(self.scripts_dir) + "/"
        if not os.path.exists(args.scripts_dir):
            print (args.scripts_dir + " not found")

        if args.ILLUMINA_reads:
            args.ILLUMINA_reads = os.path.abspath(args.ILLUMINA_reads) + "/"
        else:
            args.ILLUMINA_reads = working_dir + "reads/illumina/" + args.sample + "/"
        if not os.path.exists(args.ILLUMINA_reads):
            print (args.ILLUMINA_reads + " not found")

        if args.ONT_reads:
            args.ONT_reads = os.path.abspath(args.ONT_reads) + "/"
        else:
            args.ONT_reads =  working_dir + "reads/ont/" + args.sample + "/"
        if not os.path.exists(args.ONT_reads):
            print (args.ONT_reads + " not found")

        if args.logs_dir:
            args.logs_dir = os.path.abspath(args.logs_dir) + "/"
        else:
            args.logs_dir = args.basedir + self.logs_dir 

        if args.preprocessing_dir:
            args.preprocessing_dir = os.path.abspath(args.preprocessing_dir) + "/"
        else:
            args.preprocessing_dir = args.basedir + self.preprocessing_dir 

        if args.ILLUMINA_trim:
            args.ILLUMINA_trim = os.path.abspath(args.ILLUMINA_trim) + "/"
        else:
            args.ILLUMINA_trim = args.basedir + self.ILLUMINA_trim 

        if args.ILLUMINA_cat:
            args.ILLUMINA_cat = os.path.abspath(args.ILLUMINA_cat) 
        else:
            args.ILLUMINA_cat = args.basedir + self.ILLUMINA_cat 

        if args.ILLUMINA_qc:
            args.ILLUMINA_qc = os.path.abspath(args.ILLUMINA_qc) 
        else:
            args.ILLUMINA_qc = args.basedir+ self.ILLUMINA_qc 

        if args.bracken:
            args.bracken = os.path.abspath(args.bracken) 
        else:
            args.bracken =args.basedir + self.bracken

        if args.ONT_trim:
            args.ONT_trim = os.path.abspath(args.ONT_trim) 
        else:
            args.ONT_trim = args.basedir +  self.ONT_trim 

        if args.ONT_cat:
            args.ONT_cat = os.path.abspath(args.ONT_cat) 
        else:
            args.ONT_cat = args.basedir + self.ONT_cat 

        if args.ONT_qc:
            args.ONT_qc = os.path.abspath(args.ONT_qc) + "/"
        else:
            args.ONT_qc = args.basedir +  self.ONT_qc 

        if args.nanostats_dir:
            args.nanostats_dir = os.path.abspath(args.nanostats_dir) + "/"
        else:
            args.nanostats_dir = args.basedir + self.nanostats_dir

        if args.assembly_dir:
            args.assembly_dir = os.path.abspath(args.assembly_dir) + "/"
        else:
            args.assembly_dir = args.basedir + self.assembly_dir 

        if args.centrifuge_dir:
            args.centrifuge_dir = os.path.abspath(args.centrifuge_dir) + "/"
        else:
            args.centrifuge_dir = args.basedir + self.centrifuge_dir 

        if args.annotation_dir:
            args.annotation_dir = os.path.abspath(args.annotation_dir) + "/"
        else:
            args.annotation_dir = args.basedir + self.annotation_dir 

        if args.prokka_outdir:
            args.prokka_outdir = os.path.abspath(args.prokka_outdir) + "/"
        else:
            args.prokka_outdir =args.basedir + self.prokka_outdir 

        if args.rgi_outdir:
            args.rgi_outdir = os.path.abspath(args.rgi_outdir) + "/"
        else:
            args.rgi_outdir =args.basedir + self.rgi_outdir 

        if args.brackenDB:
            args.brackenDB = os.path.abspath(args.brackenDB) 

        if args.brackenKmers:
            args.brackenKmers = os.path.abspath(args.brackenKmers) 

        if args.assembly:
            args.assembly = os.path.abspath(args.assembly) 
        else:
            args.assembly =  args.assembly_dir + args.sample + ".v" + str(args.version) + ".assembly.fasta" 

        if args.ISFinder:
            args.ISFinder = os.path.abspath(args.ISFinder)
        else:
            args.ISFinder = args.annotation_dir + args.sample + "v" + str(args.version) + ".ISFinder_results.gff3"

        if args.resfinder:
            args.resfinder = os.path.abspath(args.resfinder)
        else:
            args.resfinder = args.annotation_dir + args.sample + "v" + str(args.version) + ".resfinder_results.gff3"

        if args.resfinder_dir:
            args.resfinder_dir = os.path.abspath(args.resfinder_dir)
        else:
            args.resfinder_dir = args.annotation_dir + "resfinder/"

        if args.resfinder_db:
            args.resfinder_db = os.path.abspath(args.resfinder_db)

        ##Checking general parameters
        if args.metadata_species==None:
            print ("No metadata species given")

        ##Assign wildcards
        if args.ILLUMINA_fastqs == None:
            for r, d, f in os.walk(args.ILLUMINA_reads):
                for file in f:
                    if re.search('1.fastq.gz', file):
                        a = file.replace('.1.fastq.gz','')
                        if args.ILLUMINA_fastqs == None:
                            args.ILLUMINA_fastqs = a
                        else:
                            args.ILLUMINA_fastqs += "," + a

        if args.ONT_fastqs == None:
            for r, d, f in os.walk(args.ONT_reads):
                for file in f:
                    if re.search('.fastq.gz', file):
                        a = file.replace('.fastq.gz','')
                        ont_barcodes.append(a)
                        if args.ONT_fastqs == None:
                            args.ONT_fastqs = a
                        else:
                            args.ONT_fastqs += "," + a               

        
###

###

    def storeGeneralParameters(self,args):
        """Updates general parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.generalParameters["configFile"] = args.configFile
        self.generalParameters["version"] = args.version
        self.generalParameters["basedir"] = args.basedir
        self.generalParameters["logs_dir"] = args.logs_dir
        self.generalParameters["sample"] = args.sample
        self.generalParameters["metadata_species"] = args.metadata_species
        self.generalParameters["QCAT_cores"] = args.Qcat_cores
        self.generalParameters["Concat_cores"] = args.Concat_cores
        self.generalParameters["nanostat_cores"] = args.nanostat_cores
        self.allParameters["Parameters"] = self.generalParameters

    def storeInputParameters(self,args):
        """Updates input parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.inputParameters["scripts_dir"] = args.scripts_dir
        self.inputParameters["ILLUMINA_reads"] = args.ILLUMINA_reads
        self.inputParameters["ONT_reads"] = args.ONT_reads
        self.allParameters ["Inputs"] = self.inputParameters

    def storeOutputParameters(self,args):
        """Updates output parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.outputParameters["preprocessing_dir"] = args.preprocessing_dir
        self.outputParameters["ILLUMINA_trim"] = args.ILLUMINA_trim
        self.outputParameters["ILLUMINA_cat"] = args.ILLUMINA_cat
        self.outputParameters["ILLUMINA_qc"] = args.ILLUMINA_qc
        self.outputParameters["bracken"] = args.bracken
        self.outputParameters["ONT_trim"] = args.ONT_trim
        self.outputParameters["ONT_cat"] = args.ONT_cat
        self.outputParameters["ONT_qc"] = args.ONT_qc
        self.outputParameters["nanostats_dir"] = args.nanostats_dir
        self.outputParameters["assembly_dir"] = args.assembly_dir
        self.outputParameters["centrifuge_dir"] = args.centrifuge_dir
        self.outputParameters["annotation_dir"] = args.annotation_dir
        self.outputParameters["prokka_outdir"] = args.prokka_outdir 
        self.outputParameters["rgi_outdir"] = args.rgi_outdir 
        self.outputParameters["assembly"] = args.assembly
        self.outputParameters["ISFinder_out"] = args.ISFinder
        self.outputParameters["resfinder_gff"] = args.resfinder
        self.allParameters ["Outputs"] = self.outputParameters

    def storeWildcardParameters(self,args):
        """Updates wildcard parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.wildcardParameters["ILLUMINA_fastqs"] = args.ILLUMINA_fastqs
        self.wildcardParameters["ONT_fastqs"] = args.ONT_fastqs
        self.allParameters ["Wildcards"] = self.wildcardParameters

    def storeTrimgaloreParameters(self,args):
        """Updates the Trim_Galore parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.trimgaloreParameters["options"] = args.trim_galore_opts
        self.trimgaloreParameters["Trim_Illumina_cores"] = args.Trim_Illumina_cores
        self.allParameters ["Trim_Galore"] = self.trimgaloreParameters

    def storeQcatParameters(self,args):
        """Updates qcat parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        for l in ont_barcodes:
          p = l.split('.')
          if len(p) > 1:
            b = p[1].replace('NB','')
            self.qcatParameters[l] = b
          else:
            b = l
            self.qcatParameters[l] = ""
        self.allParameters ["ONT_Barcodes"] = self.qcatParameters

    def storeBrackenParameters(self,args):
        """Updates the Bracken parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.brackenParameters["brackenCores"] = args.brackenCores
        self.brackenParameters["kraken_db"] = args.brackenDB
        self.brackenParameters["kraken_kmers"] = args.brackenKmers
        self.brackenParameters["additional_opts"] = args.additional_kraken2_opts
        self.allParameters ["Bracken"] = self.brackenParameters

    def storeAssemblyParameters(self,args):
        """Updates assembly parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.assemblyParameters["assemblyCores"] = args.assemblyCores
        self.allParameters ["Assembly"] = self.assemblyParameters

    def storeCentrifugeParameters(self,args):
        """Updates centrifuge parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.centrifugeParameters["centrifugeCores"] = args.centrifugeCores
        self.centrifugeParameters["db"] = args.centrifuge_db
        self.centrifugeParameters["k"] = args.centrifuge_k
        self.allParameters ["Centrifuge"] = self.centrifugeParameters

    def storeAnnotationParameters(self,args):
        """Updates annotation parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.annotationParameters["annotCores"] = args.annotCores
        self.allParameters ["Annotation"] = self.annotationParameters

    def storeISFinderParameters(self,args):
        """Updates ISFinder parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.ISFinderParameters["IS_cores"] = args.IS_cores
        self.ISFinderParameters["IS_db"] = args.IS_db
        self.allParameters ["ISFinder"] = self.ISFinderParameters

    def storeResFinderParameters(self,args):
        """Updates resfinder parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.ResFinderParameters["score"] = args.resfinder_score
        self.ResFinderParameters["db"] = args.resfinder_db
        self.ResFinderParameters["dir"] = args.resfinder_dir
        self.allParameters ["ResFinder"] = self.ResFinderParameters

#####

#1.Create object class Configuration File
configManager = CreateConfigurationFile()

#2.Create object for argument parsinng
parser = argparse.ArgumentParser(prog="create_configuration_file",
                description="Create a configuration json file for the CRE pipeline."
                )     

#2.1 Updates arguments and parsing
configManager.register_parameter(parser)

args = parser.parse_args()

#2.2 Check Parameters
configManager.check_parameters(args)

#3. store arguments to super map structure
configManager.storeGeneralParameters(args)
configManager.storeInputParameters(args)
configManager.storeOutputParameters(args)
configManager.storeWildcardParameters(args)
configManager.storeTrimgaloreParameters(args)
configManager.storeQcatParameters(args)
configManager.storeBrackenParameters(args)
configManager.storeAssemblyParameters(args)
configManager.storeCentrifugeParameters(args)
configManager.storeAnnotationParameters(args)
configManager.storeISFinderParameters(args)
configManager.storeResFinderParameters(args)

###
#4. Store JSON file
with open(args.configFile, 'w') as of:
    json.dump(configManager.allParameters, of, indent=2)




