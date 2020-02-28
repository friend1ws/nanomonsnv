#! /usr/bin/env python

import sys, math

input_file = sys.argv[1]
filter_type = sys.argv[2]

with open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')

        if int(filter_type) >= 4:
            #
            # fisher test
            if float(F[10]) < - math.log10(0.005): continue
        
        if int(filter_type) >= 2:
            # 2
            # strandness tumor
            if float(F[11]) > 0.95 or float(F[11]) < 0.05: continue 
    
            #  max_base_qual_tumor
            # if int(F[12]) < 10: continue
            # median_base_qual_tumor
            # if float(F[14]) < 10: continue
            
            # 2
            # max_plus_base_qual_tumor
            if int(F[15]) < 10: continue
            # 2
            # max_minus_base_qual_tumor
            if int(F[16]) < 10: continue

        if int(filter_type) >= 3:
            # 3
            # import pdb; pdb.set_trace()        
            # non-matched control max error ratio
            if float(F[28]) > 0.12: continue
            # 3
            # non-matched control median error ratio
            if float(F[29]) > 0.10: continue

        print('\t'.join(F))
