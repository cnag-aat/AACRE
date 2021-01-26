#!/usr/bin/env python3
import os.path


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class Gramstain(dict):
    """
    """
    gram_neg_file = ("{}/gram_neg.txt"
                     .format(os.path.abspath(os.path.dirname(__file__))))
    gram_pos_file = ("{}/gram_pos.txt"
                     .format(os.path.abspath(os.path.dirname(__file__))))

    def __init__(self):
        """
        """
        self.load_gram_file(file=Gramstain.gram_neg_file, gram="-")
        self.load_gram_file(file=Gramstain.gram_pos_file, gram="+")

    def load_gram_file(self, file, gram):
        """
        """
        with open(file, "r") as fh:
            for line in fh:
                line = line.rstrip()
                if(not line):
                    continue
                if(line.startswith("#")):
                    continue
                self[line.lower()] = gram
