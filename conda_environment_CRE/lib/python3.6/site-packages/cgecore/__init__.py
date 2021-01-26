#!/usr/bin/env python
"""
The CGE functions module
"""
from .utility import (adv_dict, copy_dir, copy_file, create_zip_dir, debug,
                      open_, file_unzipper, file_zipper, mkpath, move_file,
                      seqs_from_file, Reg, REGroup, sort2groups, load_json,
                      sort_and_distribute
                      )
from .cmdline import Program, proglist, cmd2list
from .argumentparsing import (check_file_type, get_arguments, get_string,
                              make_file_list
                              )

#####################
__version__ = "1.5.6"
__all__ = [
    "argumentparsing",
    "cmdline",
    "utility"
]

# Initiate Shared Objects
# debug = Debug()
# proglist = programlist_obj()
