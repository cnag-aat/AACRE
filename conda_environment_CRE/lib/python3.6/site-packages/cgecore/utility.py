#!/usr/bin/env python3
""" THIS MODULE CONTAINS ALL THE SHARED WRAPPER FUNCTIONS """
################################################################################
#                              CGE FUNCTION MODULE                             #
################################################################################
# This script is part of the CGE Pipeline structure
import sys, os, gzip, shutil, glob, re, json
from subprocess import Popen
from zipfile import ZipFile
from contextlib import closing

############# CLASSES #############
class Debug():
   """ Debug object, keeps track of all debug related matters, and provides easy
   access to logging debug text

   USAGE
      import os, sys
      debug = Debug()
      debug.setup(debug=True, logfile=None, stdout=None, stderr=None)
      debug.print_out('hello')
      debug.log_no_newline('Oh')
      debug.log_no_newline(' my')
      debug.log(' god!')
      # The graceful exit is only meant to be used with the integrated CGE platform
      debug.graceful_exit('error message')
   IMPORTING THE DEBUG OBJECT
      from cge.utility import debug
   """
   def __init__(self):
      """  """
      self.debug = False
      self.logfile = sys.stderr
      self.stdout = sys.stdout
      self.stderr = sys.stderr
      self.caught_error = None
   def setup(self, debug=None, logfile=None, stdout=None, stderr=None):
      """  """
      if debug   is not None: self.debug   = debug
      if logfile is not None: self.logfile = logfile
      if stdout  is not None: self.stdout  = stdout
      if stderr  is not None: self.stderr  = stderr
   def print_out(self, *lst):
      """ Print list of strings to the predefined stdout. """
      self.print2file(self.stdout, True, True, *lst)
   def print_err(self, *lst):
      """ Print list of strings to the predefined stdout. """
      self.print2file(self.stderr, False, True, *lst)
   def print2file(self, logfile, print2screen, addLineFeed, *lst):
      """ This function prints to the screen and logs to a file, all the strings
      given.
      # print2screen eg. True, *lst is a commaseparated list of strings
      """
      if addLineFeed:
         linefeed = '\n'
      else: linefeed = ''
      if print2screen: print(linefeed.join(str(string) for string in lst))
      try: file_instance = isinstance(logfile, file)
      except NameError as e:
         from io import IOBase
         try: file_instance = isinstance(logfile, IOBase)
         except: raise e
      if file_instance:
         logfile.write(linefeed.join(str(string) for string in lst) + linefeed)
      elif isinstance(logfile, str) and os.path.exists(logfile):
         with open_(logfile, 'a') as f:
            f.write(linefeed.join(str(string) for string in lst) + linefeed)
      elif not print2screen: # Print to screen if there is no outputfile
         print(linefeed.join(str(string) for string in lst))
   def log(self, *lst):
      """ Print list of strings to the predefined logfile if debug is set. and
      sets the caught_error message if an error is found
      """
      self.print2file(self.logfile, self.debug, True, *lst)
      if 'Error' in '\n'.join([str(x) for x in lst]):
         self.caught_error = '\n'.join([str(x) for x in lst])
   def log_no_newline(self, msg):
      """ print the message to the predefined log file without newline """
      self.print2file(self.logfile, False, False, msg)
   def graceful_exit(self, msg):
      """ This function Tries to update the MSQL database before exiting. """
      # Print stored errors to stderr
      if self.caught_error:
         self.print2file(self.stderr, False, False, self.caught_error)
      # Kill process with error message
      self.log(msg)
      sys.exit(1)

class adv_dict(dict):
   """ This class expands on the dictionary class by adding the gettree class
   method.
   """
   def get_tree(self, list_of_keys):
      """ gettree will extract the value from a nested tree
      
      INPUT
         list_of_keys: a list of keys ie. ['key1', 'key2']
      USAGE
      >>> # Access the value for key2 within the nested dictionary
      >>> adv_dict({'key1': {'key2': 'value'}}).gettree(['key1', 'key2'])
      'value'
      """
      cur_obj = self
      for key in list_of_keys:
         cur_obj = cur_obj.get(key)
         if not cur_obj: break
      return cur_obj
   def invert(self):
      ''' Return inverse mapping of dictionary with sorted values.
      USAGE
         >>> # Switch the keys and values
         >>> adv_dict({
         ...     'A': [1, 2, 3],
         ...     'B': [4, 2],
         ...     'C': [1, 4],
         ... }).invert()
         {1: ['A', 'C'], 2: ['A', 'B'], 3: ['A'], 4: ['B', 'C']}
      '''
      inv_map = {}
      for k, v in self.items():
         if sys.version_info < (3, 0):
            acceptable_v_instance = isinstance(v, (str, int, float, long))
         else:
            acceptable_v_instance = isinstance(v, (str, int, float))
         if acceptable_v_instance: v = [v]
         elif not isinstance(v, list):
            raise Exception('Error: Non supported value format! Values may only'
                            ' be numerical, strings, or lists of numbers and '
                            'strings.')
         for val in v:
            inv_map[val] = inv_map.get(val, [])
            inv_map[val].append(k)
            inv_map[val].sort()
      return inv_map

class Reg:
   """
   NAME:        Reg - Extended Regular Expression Handler
   AUTHOR:      Martin Thomsen
   DESCRIPTION:
      This class enables a simplistic usage of regular expression to get
      contained groups in a match statement. But it also allows to do some of
      the normal re call, such as findall and sub.
   DEPENDENCIES:
      re (regular expression module)
   USAGE:
      >>> import re
      >>> RegEx = Reg(pattern, flag)
      >>> if RegEx.match(string):
      >>>    RegEx.getgroup(index)
   EXAMPLE:
      >>> RegEx = Reg('[^a]*(a)[^b]*(b)[^Y]*(Y)(a)?', 'I')
      >>> if RegEx.match('aBcdefgHIJKLmnOpqrstuvwxyz'):
      ...    print(RegEx.getgroup(0), # index=0 -> full match
      ...          RegEx.getgroup(1),
      ...          RegEx.getgroup(2),
      ...          RegEx.getgroup(3),
      ...          RegEx.getgroup(4))
      ...
      ('aBcdefgHIJKLmnOpqrstuvwxy', 'a', 'B', 'y', None)
      # NIFTY SUBSTITUTION LOOP
      >>> string = 'There are {%=count%} {%=animal%} on the {%=location%}!'
      >>> # Dictionary containing place-holders and values
      ... # (make sure all placeholders in the string is included!)
      ... d = { 'count': 5, 'animal': 'cows', 'location': 'Battle Field' }
      >>> # RE Object which matches place-holders in the string
      ... tmpPH = Reg('\{\%\=(\w+)\%\}', 'I')
      >>> # substitute all placeholders
      ... while tmpPH.match(string): string = tmpPH.sub(str(d[tmpPH.getgroup(1)]), string, 1)
      ...
      >>> print(string)
      There are 5 cows on the Battle Field!
      """
   def __init__(self, pattern, *flags):
      sd = {'T': 1, 'I': 2, 'L':4,'M':8,'S':16,'U':32,'X':64}
      try:
         flag = sum([sd[f] if f in sd else int(f) for f in set(flags)]) if flags else 0
      except:
         for f in flags:
            if not isinstance(f, int) and not f in sd:
               raise Exception("Error: Unrecognised flag argument '%s' for Reg call."%f)
         flag=0
      if flag: self.re = re.compile(pattern, flag)
      else: self.re = re.compile(pattern)
      self.matches = None
   def sub(self, replace, string, count=0):
      """ returns new string where the matching cases (limited by the count) in
      the string is replaced. """
      return self.re.sub(replace, string, count)
   def find_all(self, s):
      """ Finds all matches in the string and returns them in a tuple. """
      return self.re.findall(s)
   def match(self, s):
      """ Matches the string to the stored regular expression, and stores all
      groups in mathches. Returns False on negative match. """
      self.matches = self.re.search(s)
      return self.matches
   def get_group(self, x):
      """ Returns requested subgroup. """
      return self.matches.group(x)
   def get_groups(self):
      """ Returns all subgroups. """
      return self.matches.groups()

class REGroup():
   """ Regular Expression group object

   This class simplyfies the use of groups for the Sort2Groups function.
   """
   def __init__(self, pattern, flags=''):
      self.re = Reg(pattern, flags)
      self.list = []
   def match(self, s):
      """ Matching the pattern to the input string, returns True/False and
          saves the matched string in the internal list
      """
      if self.re.match(s):
         self.list.append(s)
         return True
      else: return False

############# ITERATORS #############
def seqs_from_file(filename, exit_on_err=False, return_qual=False):
   """Extract sequences from a file
   
   Name:
      seqs_from_file
   Author(s):
      Martin C F Thomsen
   Date:
      18 Jul 2013
   Description:
      Iterator which extract sequence data from the input file
   Args:
      filename: string which contain a path to the input file
   Supported Formats:
      fasta, fastq
   
   USAGE:
   >>> import os, sys
   >>> # Create fasta test file
   >>> file_content = ('>head1 desc1\nthis_is_seq_1\n>head2 desc2\n'
                       'this_is_seq_2\n>head3 desc3\nthis_is_seq_3\n')
   >>> with open_('test.fsa', 'w') as f: f.write(file_content)
   >>> # Parse and print the fasta file
   >>> for seq, name, desc in SeqsFromFile('test.fsa'):
   ...    print ">%s %s\n%s"%(name, desc, seq)
   ...
   >head1 desc1
   this_is_seq_1
   >head2 desc2
   this_is_seq_2
   >head3 desc3
   this_is_seq_3
   """
   # VALIDATE INPUT
   if not isinstance(filename, str):
      msg = 'Filename has to be a string.'
      if exit_on_err:
         sys.stderr.write('Error: %s\n'%msg)
         sys.exit(1)
      else: raise IOError(msg)
   if not os.path.exists(filename):
      msg = 'File "%s" does not exist.'%filename
      if exit_on_err:
         sys.stderr.write('Error: %s\n'%msg)
         sys.exit(1)
      else: raise IOError(msg)
   
   # EXTRACT DATA
   with open_(filename,"rt") as f:
      query_seq_segments = []
      seq, name, desc, qual = '', '', '', ''
      add_segment = query_seq_segments.append
      for l in f:
         if len(l.strip()) == 0: continue
         #sys.stderr.write("%s\n"%line)
         fields=l.strip().split()
         if l.startswith(">"):
            # FASTA HEADER FOUND
            if query_seq_segments != []:
               # YIELD SEQUENCE AND RESET
               seq = ''.join(query_seq_segments)
               yield (seq, name, desc)
               seq, name, desc = '', '', ''
               del query_seq_segments[:]
            name = fields[0][1:]
            desc = ' '.join(fields[1:])
         
         elif l.startswith("@"):
            # FASTQ HEADER FOUND
            name = fields[0][1:]
            desc = ' '.join(fields[1:])
            try:
               # EXTRACT FASTQ SEQUENCE
               seq  = next(f).strip().split()[0]
               # SKIP SECOND HEADER LINE AND QUALITY SCORES
               l = next(f)
               qual = next(f).strip() # Qualities
            except:
               break
            else:
               # YIELD SEQUENCE AND RESET
               if return_qual:
                  yield (seq, qual, name, desc)
               else:
                  yield (seq, name, desc)
               seq, name, desc, qual = '', '', '', ''
         
         elif len(fields[0])>0:
            # EXTRACT FASTA SEQUENCE
            add_segment(fields[0])
      
      # CHECK FOR LAST FASTA SEQUENCE
      if query_seq_segments != []:
         # YIELD SEQUENCE
         seq = ''.join(query_seq_segments)
         yield (seq, name, desc)

############# FUNCTIONS #############
def open_(filename, mode=None, compresslevel=9):
   """Switch for both open() and gzip.open().
   
   Determines if the file is normal or gzipped by looking at the file
   extension.
   
   The filename argument is required; mode defaults to 'rb' for gzip and 'r'
   for normal and compresslevel defaults to 9 for gzip.
   
   >>> import gzip
   >>> from contextlib import closing
   >>> with open_(filename) as f:
   ...     f.read()
   """
   if filename[-3:] == '.gz':
      if mode is None: mode = 'rt'
      return closing(gzip.open(filename, mode, compresslevel))
   else:
      if mode is None: mode = 'r'
      return open(filename, mode)

def load_json(json_object):
   ''' Load json from file or file name '''
   content = None
   if isinstance(json_object, str) and os.path.exists(json_object):
      with open_(json_object) as f:
         try:
            content = json.load(f)
         except Exception as e:
            debug.log("Warning: Content of '%s' file is not json."%f.name)
   elif hasattr(json_object, 'read'):
      try:
         content = json.load(json_object)
      except Exception as e:
         debug.log("Warning: Content of '%s' file is not json."%json_object.name)
   else:
      debug.log("%s\nWarning: Object type invalid!"%json_object)
   return content

def sort2groups(array, gpat=['_R1','_R2']):
   """ Sort an array of strings to groups by patterns """
   groups = [REGroup(gp) for gp in gpat]
   unmatched = []
   for item in array:
      matched = False
      for m in groups:
         if m.match(item):
            matched = True
            break
      if not matched: unmatched.append(item)
   return [sorted(m.list) for m in groups], sorted(unmatched)

def sort_and_distribute(array, splits=2):
   """ Sort an array of strings to groups by alphabetically continuous
       distribution
   """
   if not isinstance(array, (list,tuple)): raise TypeError("array must be a list")
   if not isinstance(splits, int): raise TypeError("splits must be an integer")
   remaining = sorted(array)
   if sys.version_info < (3, 0):
      myrange = xrange(splits)
   else:
      myrange = range(splits)
   groups = [[] for i in myrange]
   while len(remaining) > 0:
      for i in myrange:
         if len(remaining) > 0: groups[i].append(remaining.pop(0))
   return groups

def mkpath(filepath, permissions=0o777):
   """ This function executes a mkdir command for filepath and with permissions
   (octal number with leading 0 or string only)
   # eg. mkpath("path/to/file", "0o775")
   """
   # Converting string of octal to integer, if string is given.
   if isinstance(permissions, str):
      permissions = sum([int(x)*8**i for i,x in enumerate(reversed(permissions))])
   # Creating directory
   if not os.path.exists(filepath):
      debug.log("Creating Directory %s (permissions: %s)"%(
         filepath, permissions))
      os.makedirs(filepath, permissions)
   else:
      debug.log("Warning: The directory "+ filepath +" already exists")
   return filepath

def create_zip_dir(zipfile_path, *file_list):
   """ This function creates a zipfile located in zipFilePath with the files in
   the file list
   # fileList can be both a comma separated list or an array
   """
   try:
      if isinstance(file_list, (list, tuple)): #unfolding list of list or tuple
         if len(file_list) == 1:
            if isinstance(file_list[0], (list, tuple)): file_list = file_list[0]
      #converting string to iterable list
      if isinstance(file_list, str): file_list = [file_list]
      if file_list:
         with ZipFile(zipfile_path, 'w') as zf:
            for cur_file in file_list:
               if '/' in cur_file:
                  os.chdir('/'.join(cur_file.split('/')[:-1]))
               elif '/' in zipfile_path:
                  os.chdir('/'.join(zipfile_path.split('/')[:-1]))
               zf.write(cur_file.split('/')[-1])
      else:
         debug.log('Error: No Files in list!',zipfile_path+' was not created!')
   except Exception as e:
      debug.log('Error: Could not create zip dir! argtype: '+
                 str(type(file_list)), "FileList: "+ str(file_list),
                 "Errormessage: "+ str(e))

def file_zipper(root_dir):
   """ This function will zip the files created in the runroot directory and
   subdirectories """
   # FINDING AND ZIPPING UNZIPPED FILES
   for root, dirs, files in os.walk(root_dir, topdown=False):
      if root != "":
         if root[-1] != '/': root += '/'
         for current_file in files:
            filepath = "%s/%s"%(root, current_file)
            try:
               file_size = os.path.getsize(filepath)
            except Exception as e:
               file_size = 0
               debug.log('Error: file_zipper failed to zip following file '+filepath, e)
            # Excluding small files, gzipped files and links
            if (         file_size > 50
                 and     current_file[-3:] != ".gz"
                 and not os.path.islink(filepath)
               ):
               if current_file[-4:] == ".zip":
                  # Unzip file
                  ec = Popen('unzip -qq "%s" -d %s > /dev/null 2>&1'%(filepath, root), shell=True).wait()
                  if ec > 0:
                     debug.log('Error: fileZipper failed to unzip following file %s'%filepath)
                     continue
                  else:
                     ec = Popen('rm -f "%s" > /dev/null 2>&1'%(filepath), shell=True).wait()
                     if ec > 0: debug.log('Error: fileZipper failed to delete the original zip file (%s)'%filepath)
                     filepath = filepath[:-4]
                  # Saving a gzipped version
                  with open_(filepath, 'rb') as f, open_(filepath+".gz", 'wb', 9) as gz:
                     gz.writelines(f)
                  # Deleting old (non-zipped) file
                  try: os.remove(filepath)
                  except OSError as e:
                     debug.log(("WARNING! The file %s could not be "
                                    "removed!\n%s")%(current_file, e))

def file_unzipper(directory):
   """ This function will unzip all files in the runroot directory and
   subdirectories
   """
   debug.log("Unzipping directory (%s)..."%directory)
   #FINDING AND UNZIPPING ZIPPED FILES
   for root, dirs, files in os.walk(directory, topdown=False):
      if root != "":
         orig_dir = os.getcwd()
         os.chdir(directory)
         Popen('gunzip -q -f *.gz > /dev/null 2>&1', shell=True).wait()
         Popen('unzip -qq -o "*.zip" > /dev/null 2>&1', shell=True).wait()
         Popen('rm -f *.zip > /dev/null 2>&1', shell=True).wait()
         os.chdir(orig_dir)

def move_file(src, dst):
   """ this function will simply move the file from the source path to the dest
   path given as input
   """
   # Sanity checkpoint
   src = re.sub('[^\w/\-\.\*]', '', src)
   dst = re.sub('[^\w/\-\.\*]', '', dst)
   if len(re.sub('[\W]', '', src)) < 5 or len(re.sub('[\W]', '', dst)) < 5:
      debug.log("Error: Moving file failed. Provided paths are invalid! src='%s' dst='%s'"%(src, dst))
   else:
      # Check destination
      check = False
      if dst[-1] == '/':
         if os.path.exists(dst):
            check = True # Valid Dir
         else:
            debug.log("Error: Moving file failed. Destination directory does not exist (%s)"%(dst)) #DEBUG
      elif os.path.exists(dst):
         if os.path.isdir(dst):
            check = True # Valid Dir
            dst += '/' # Add missing slash
         else:
            debug.log("Error: Moving file failed. %s exists!"%dst)
      elif os.path.exists(os.path.dirname(dst)):
         check = True # Valid file path
      else:
         debug.log("Error: Moving file failed. %s is an invalid distination!"%dst)
      if check:
         # Check source
         files = glob.glob(src)
         if len(files) != 0:
            debug.log("Moving File(s)...", "Move from %s"%src, "to %s"%dst)
            for file_ in files:
               # Check if file contains invalid symbols:
               invalid_chars = re.findall('[^\w/\-\.\*]', os.path.basename(file_))
               if invalid_chars:
                  debug.graceful_exit(("Error: File %s contains invalid "
                                      "characters %s!"
                                      )%(os.path.basename(file_), invalid_chars))
                  continue
               # Check file exists
               if os.path.isfile(file_):
                  debug.log("Moving file: %s"%file_)
                  shutil.move(file_, dst)
               else:
                  debug.log("Error: Moving file failed. %s is not a regular file!"%file_)
         else: debug.log("Error: Moving file failed. No files were found! (%s)"%src)

def copy_file(src, dst, ignore=None):
   """ this function will simply copy the file from the source path to the dest
   path given as input
   """
   # Sanity checkpoint
   src = re.sub('[^\w/\-\.\*]', '', src)
   dst = re.sub('[^\w/\-\.\*]', '', dst)
   if len(re.sub('[\W]', '', src)) < 5 or len(re.sub('[\W]', '', dst)) < 5:
      debug.log("Error: Copying file failed. Provided paths are invalid! src='%s' dst='%s'"%(src, dst))
   else:
      # Check destination
      check = False
      if dst[-1] == '/':
         if os.path.exists(dst):
            check = True # Valid Dir
         else:
            debug.log("Error: Copying file failed. Destination directory does not exist (%s)"%(dst)) #DEBUG
      elif os.path.exists(dst):
         if os.path.isdir(dst):
            check = True # Valid Dir
            dst += '/' # Add missing slash
         else:
            debug.log("Error: Copying file failed. %s exists!"%dst)
      elif os.path.exists(os.path.dirname(dst)):
         check = True # Valid file path
      else:
         debug.log("Error: Copying file failed. %s is an invalid distination!"%dst)
      if check:
         # Check source
         files = glob.glob(src)
         if ignore is not None: files = [fil for fil in files if not ignore in fil]
         if len(files) != 0:
            debug.log("Copying File(s)...", "Copy from %s"%src, "to %s"%dst) #DEBUG
            for file_ in files:
               # Check file exists
               if os.path.isfile(file_):
                  debug.log("Copying file: %s"%file_) #DEBUG
                  shutil.copy(file_, dst)
               else:
                  debug.log("Error: Copying file failed. %s is not a regular file!"%file_) #DEBUG
         else: debug.log("Error: Copying file failed. No files were found! (%s)"%src) #DEBUG

def copy_dir(src, dst):
   """ this function will simply copy the file from the source path to the dest
   path given as input
   """
   try:
      debug.log("copy dir from "+ src, "to "+ dst)
      shutil.copytree(src, dst)
   except Exception as e:
      debug.log("Error: happened while copying!\n%s\n"%e)

# Initiate Shared Global Objects
debug = Debug()
