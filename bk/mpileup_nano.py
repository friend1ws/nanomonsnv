#! /usr/bin/env python

import re, sys, math, statistics
import scipy.stats as stats

from utils import check_pileup_record


def proc_mpileup(input_file, output_file, min_variant_num_tumor = 5, min_depth_tumor = 8, min_variant_ratio_tumor = 0.1, 
                 max_variant_num_ctrl = 2, min_depth_ctrl = 10, max_variant_ratio_ctrl = 0.05, min_minus_log10_fisher_pvalue = 1.30103):

    hout = open(output_file, 'w') 

    # print('\t'.join(["Chr", "Start", "End", "Ref", "Alt", "Original_Ref", "Variant_Ratio", "Depth", 
    #                  "Variant_Num", "Strand_Ratio", "Variant_Num_Info"]), file = hout)

    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            if int(F[3]) < min_depth_tumor: continue        
            # first simple check
            total_variant_num_tumor = int(F[3]) - F[4].count('.') - F[4].count(',') - F[4].count('>') - F[4].count('<') + F[4].count('+') + F[4].count('-')
            if total_variant_num_tumor < min_variant_num_tumor: continue
            if total_variant_num_tumor / int(F[3]) < min_variant_ratio_tumor: continue

            if int(F[6]) < 8: continue

            # total_variant_num_ctrl = int(F[6]) - F[7].count('.') - F[7].count(',') - F[7].count('>') - F[7].count('<') + F[7].count('+') + F[7].count('-')
            # if total_variant_num_ctrl > max_variant_num_ctrl: continue
            # if total_variant_num_ctrl / int(F[6]) > max_variant_ratio_ctrl: continue
            
            # if F[1] == "16640410":
            #     import pdb; pdb.set_trace()

            var2num_tumor, ovar2qual_tumor, depth_tumor = check_pileup_record(F[2], F[4], F[5])

            if depth_tumor < min_depth_tumor: continue

            bvar = ''
            bmis_rate_tumor = 0
            for var in ['A', 'C', 'G', 'T']:
                if var == F[2].upper(): continue
                if var2num_tumor[var] < min_variant_num_tumor: continue
                quals = sorted(ovar2qual_tumor[var.lower()] + ovar2qual_tumor[var.upper()], reverse = True)
                # if quals[0] < 20 or quals[1] < 15 or quals[2] < 15: continue
                # if max(ovar2qual_tumor[var.lower()] + ovar2qual_tumor[var.upper()]) < 15: continue
                cur_rate = float(var2num_tumor[var]) / depth_tumor
                if cur_rate > bmis_rate_tumor:
                    bmis_rate_tumor = cur_rate
                    bvar = var

            if bmis_rate_tumor < min_variant_ratio_tumor: continue
            if var2num_tumor[bvar] < min_variant_num_tumor: continue

            strandness_tumor = len(ovar2qual_tumor[bvar.lower()]) / (len(ovar2qual_tumor[bvar.lower()]) + len(ovar2qual_tumor[bvar.upper()]))
            quals = sorted(ovar2qual_tumor[bvar.lower()] + ovar2qual_tumor[bvar.upper()], reverse = True)
            if len(quals) < 1:
                import pdb; pdb.set_trace()
            max_base_qual_tumor = quals[0]
            second_max_base_qual_tumor = quals[1]
            median_base_qual_tumor = statistics.median(quals)
            max_plus_base_qual_tumor = sorted(ovar2qual_tumor[bvar.upper()], reverse = True)[0] if len(ovar2qual_tumor[bvar.upper()]) > 0 else "---"
            max_minus_base_qual_tumor = sorted(ovar2qual_tumor[bvar.lower()], reverse = True)[0] if len(ovar2qual_tumor[bvar.lower()]) > 0 else "---"

            var_info = []
            for bb in ['A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n']:
                var_info.append(','.join([str(x) for x in ovar2qual_tumor[bb]]) if len(ovar2qual_tumor[bb]) > 0 else '---')

            pos, ref, alt = F[1], F[2], bvar


            ##########
            # control check
            var2num_ctrl, ovar2qual_ctrl, depth_ctrl = check_pileup_record(F[2], F[7], F[8])
            
            bmis_rate_ctrl = float(var2num_ctrl[bvar]) / depth_ctrl

            oddsratio, pvalue = stats.fisher_exact([[depth_tumor - var2num_tumor[bvar], var2num_tumor[bvar]], 
                                                    [depth_ctrl - var2num_ctrl[bvar], var2num_ctrl[bvar]]])

            if -math.log10(pvalue) < min_minus_log10_fisher_pvalue: continue
            if bmis_rate_ctrl > max_variant_ratio_ctrl: continue
            if var2num_ctrl[bvar] > max_variant_num_ctrl: continue
            if depth_ctrl < min_depth_ctrl: continue

            print('\t'.join([F[0], pos, ref, alt, str(depth_tumor), str(var2num_tumor[bvar]), str(round(bmis_rate_tumor, 4)), 
                             str(depth_ctrl), str(var2num_ctrl[bvar]), str(round(bmis_rate_ctrl, 4)), str(round(-math.log10(pvalue), 4)),
                             str(round(strandness_tumor, 4)), str(max_base_qual_tumor), str(second_max_base_qual_tumor),
                             str(median_base_qual_tumor), str(max_plus_base_qual_tumor), str(max_minus_base_qual_tumor)] + var_info), file = hout)


    hout.close()


if __name__ == "__main__":

    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    proc_mpileup(input_file, output_file)
    
