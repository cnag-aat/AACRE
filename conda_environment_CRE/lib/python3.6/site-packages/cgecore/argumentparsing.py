#!/usr/bin/env python3
""" THIS MODULE CONTAINS ALL THE SHARED WRAPPER FUNCTIONS """
################################################################################
#                              CGE FUNCTION MODULE                             #
################################################################################
# This script is part of the CGE Pipeline structure
import os, sys
from subprocess import Popen, PIPE
from argparse import ArgumentParser, RawDescriptionHelpFormatter, SUPPRESS

# CGE modules
from .utility import debug, open_

def get_string(string):
   """ This function checks if a path was given as string, and tries to read the
       file and return the string.
   """
   truestring = string
   if string is not None:
      if '/' in string:
         if os.path.isfile(string):
            try:
               with open_(string,'r') as f:
                  truestring = ' '.join(line.strip() for line in f)
            except: pass
      if truestring.strip() == '': truestring = None
   return truestring

def get_arguments(options):
   """ This function handles and validates the wrapper arguments. """
   # These the next couple of lines defines the header of the Help output
   parser = ArgumentParser(
      formatter_class=RawDescriptionHelpFormatter,
      usage=("""%(prog)s
--------------------------------------------------------------------------------
"""),
      description=("""
Service Wrapper
===============
This is the service wrapper script, which is a part of the CGE services.
Read the online manual for help.
A list of all published services can be found at:
cge.cbs.dtu.dk/services

"""), epilog=("""
--------------------------------------------------------------------------------
      """))

   #ADDING ARGUMENTS
   setarg = parser.add_argument
   #SERVICE SPECIFIC ARGUMENTS
   if isinstance(options, str):
      options = [[x for i,x in enumerate(line.split()) if i in [1,2]] for line in options.split('\n') if len(line)>0]
      for o in options:
         try:
            setarg(o[1], type=str, dest=o[0], default=None, help=SUPPRESS)
         except:
            None
   else:
      for o in options:
         if o[2] is True:
            # Handle negative flags
            setarg(o[0], action="store_false", dest=o[1], default=o[2],
                   help=o[3])
         elif o[2] is False:
            # Handle positive flags
            setarg(o[0], action="store_true", dest=o[1], default=o[2],
                   help=o[3])
         else:
            help_ = o[3] if o[2] is None else "%s [%s]"%(o[3], '%(default)s')
            setarg(o[0], type=str, dest=o[1], default=o[2],
                   help=help_)
   # VALIDATION OF ARGUMENTS
   args = parser.parse_args()
   debug.log("ARGS: %s"%args)
   return args

def check_file_type(files):
   """ Check whether the input files are in fasta format, reads format or
       other/mix formats.
   """
   all_are_fasta = True
   all_are_reads = True
   all_are_empty = True
   if sys.version_info < (3, 0):
      if isinstance(files, (str, unicode)): files = [files]
   else:
      if isinstance(files, str): files = [files]
   for file_ in files:
      debug.log('Checking file type: %s'%file_)
      # Check if file is empty
      if os.stat(file_).st_size == 0: continue
      else: all_are_empty = False
      with open_(file_) as f:
         fc = f.readline()[0]
         if fc != "@": all_are_reads = False
         if fc != ">": all_are_fasta = False
   if all_are_empty:   return 'empty'
   elif all_are_fasta: return 'fasta'
   elif all_are_reads: return 'fastq'
   else: return 'other'

def make_file_list(upload_path):
   """ This function returns list of files in the given dir """
   newlist = []
   for el in sorted(os.listdir(upload_path)):
      if ' ' in el:
         raise Exception('Error: Spaces are not allowed in file names!\n')
      newlist.append(os.path.normpath(upload_path+'/'+el))
   debug.log('InputFiles: %s\n'%newlist)
   return newlist
