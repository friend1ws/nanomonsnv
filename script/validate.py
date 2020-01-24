#! /usr/bin/env python

import sys, math, argparse
import scipy.stats as stats

from utils import check_pileup_record


def short_read_validate(var_file, illumina_file, output_file):

    chr_pos2info = {}
    with open(var_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            chr_pos2info[F[0] + '\t' +  F[1]] = '\t'.join(F)
                
    hout = open(output_file, 'w')
   
    with open(illumina_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0] + '\t' + F[1] not in chr_pos2info: continue
    
            var_info = chr_pos2info[F[0] + '\t' + F[1]].split('\t')

            # import pdb; pdb.set_trace()

            var2num_tumor, ovar2qual_tumor, depth_tumor = check_pileup_record(F[2], F[4], F[5])
            var_num_tumor = var2num_tumor[var_info[3]]
            vaf_tumor = str(round(float(var_num_tumor) / depth_tumor, 4)) if depth_tumor > 0 else "---"

            var2num_ctrl, ovar2qual_ctrl, depth_ctrl = check_pileup_record(F[2], F[7], F[8])
            var_num_ctrl = var2num_ctrl[var_info[3]]
            vaf_ctrl = str(round(float(var_num_ctrl) / depth_ctrl, 4)) if depth_ctrl > 0 else "---"
            
            oddsratio, pvalue = stats.fisher_exact([[depth_tumor - var_num_tumor, var_num_tumor], 
                                                    [depth_ctrl - var_num_ctrl, var_num_ctrl]])

            is_validated = "Ambiguous"
            if depth_tumor >= 20 and depth_ctrl >= 20:
                if vaf_tumor != "---" and float(vaf_tumor) >= 0.07 and var_num_tumor >= 3 and (vaf_ctrl == "---" or float(vaf_ctrl) < 0.02) and var_num_ctrl <= 2 and pvalue < 0.05:
                    is_validated = "True"
                else:
                    is_validated = "False"

            print(chr_pos2info[F[0] + '\t' +  F[1]] + '\t' + \
              str(depth_tumor) + '\t' + str(var_num_tumor) + '\t' + str(vaf_tumor) + '\t' + \
              str(depth_ctrl) + '\t' + str(var_num_ctrl) + '\t' + str(vaf_ctrl) + '\t' + \
              str(round(-math.log10(pvalue), 4)) + '\t' + is_validated, file = hout)          


    hout.close()


if __name__ == "__main__":

    import sys
    var_file = sys.argv[1]
    illumina_file = sys.argv[2]
    output_file = sys.argv[3]
    # proc_mpileup(input_file, output_file, min_variant_num_tumor = 5, min_variant_ratio_tumor = 0.15)
    short_read_validate(var_file, illumina_file, output_file)


