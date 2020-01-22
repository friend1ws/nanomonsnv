#! /usr/bin/env python

import sys, math

input_file = sys.argv[1]

with open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')

        if float(F[10]) < - math.log10(0.005): continue
        if float(F[11]) > 0.95 or float(F[11]) < 0.05: continue

        if int(F[12]) < 10: continue
        if float(F[14]) < 10: continue
        if int(F[15]) < 10: continue
        if int(F[16]) < 10: continue

        # import pdb; pdb.set_trace()        
        if float(F[28]) > 0.12: continue
        if float(F[29]) > 0.10: continue

        print('\t'.join(F))
