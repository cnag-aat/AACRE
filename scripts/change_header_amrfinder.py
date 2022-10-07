#!/usr/bin/env python3

#Author: Jessica Gomez-Garrido, CNAG-CRG.
#Contact email: jessica.gomez@cnag.crg.eu

import sys

for line in sys.stdin:
    nline = []
    line_col = line.rstrip().split('\t')
    if line_col[0] == "Protein identifier":
        for i in line_col:
            n = i.lower().replace(' ', '_')
            nline.append(n)
           # print (n)
        nline[0] = 'gene'
        fline = '\t'.join(nline)
        print (fline)
    else:
        print (line.rstrip())



    
