#! /usr/bin/env python

import sys, re, statistics
import scipy.stats as stats

from utils import check_pileup_record


def proc_ctrl_pileup_record_info(ref, alt, bases, qualities):

    var2num_ctrl, ovar2qual_ctrl, depth_ctrl = check_pileup_record(ref, bases, qualities)
    depth_p_ctrl, depth_n_ctrl = 0, 0
    for var in ['A', 'C', 'G', 'T', 'N']:
        depth_p_ctrl = depth_p_ctrl + len(ovar2qual_ctrl[var])
        depth_n_ctrl = depth_n_ctrl + len(ovar2qual_ctrl[var.lower()])
   
    var_p_ctrl, var_n_ctrl = len(ovar2qual_ctrl[alt.upper()]), len(ovar2qual_ctrl[alt.lower()])

    return([depth_p_ctrl, var_p_ctrl, depth_n_ctrl, var_n_ctrl])


def short_read_validate(var_file, control_file, output_file):

    chr_pos2info = {}
    with open(var_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            chr_pos2info[F[0] + '\t' +  F[1]] = '\t'.join(F)
                
    hout = open(output_file, 'w')
   
    with open(control_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0] + '\t' + F[1] not in chr_pos2info: continue
    
            var_info = chr_pos2info[F[0] + '\t' + F[1]].split('\t')

            # import pdb; pdb.set_trace()

            depth_p_ctrl1, var_p_ctrl1, depth_n_ctrl1, var_n_ctrl1 = proc_ctrl_pileup_record_info(F[2], var_info[3], F[4], F[5])
            depth_p_ctrl2, var_p_ctrl2, depth_n_ctrl2, var_n_ctrl2 = proc_ctrl_pileup_record_info(F[2], var_info[3], F[7], F[8])
            depth_p_ctrl3, var_p_ctrl3, depth_n_ctrl3, var_n_ctrl3 = proc_ctrl_pileup_record_info(F[2], var_info[3], F[10], F[11])
            depth_p_ctrl4, var_p_ctrl4, depth_n_ctrl4, var_n_ctrl4 = proc_ctrl_pileup_record_info(F[2], var_info[3], F[13], F[14])
            depth_p_ctrl5, var_p_ctrl5, depth_n_ctrl5, var_n_ctrl5 = proc_ctrl_pileup_record_info(F[2], var_info[3], F[16], F[17])

            ctrl_info = ','.join([str(depth_p_ctrl1), str(var_p_ctrl1), str(depth_n_ctrl1), str(var_n_ctrl1)]) + ';' + \
                        ','.join([str(depth_p_ctrl2), str(var_p_ctrl2), str(depth_n_ctrl2), str(var_n_ctrl2)]) + ';' + \
                        ','.join([str(depth_p_ctrl3), str(var_p_ctrl3), str(depth_n_ctrl3), str(var_n_ctrl3)]) + ';' + \
                        ','.join([str(depth_p_ctrl4), str(var_p_ctrl4), str(depth_n_ctrl4), str(var_n_ctrl4)]) + ';' + \
                        ','.join([str(depth_p_ctrl5), str(var_p_ctrl5), str(depth_n_ctrl5), str(var_n_ctrl5)])

            vaf_ctrl_2 = float(var_p_ctrl2 + var_n_ctrl2) / (depth_p_ctrl2 + depth_n_ctrl2 + 1)
            vaf_ctrl_3 = float(var_p_ctrl3 + var_n_ctrl3) / (depth_p_ctrl3 + depth_n_ctrl3 + 1) 
            vaf_ctrl_4 = float(var_p_ctrl4 + var_n_ctrl4) / (depth_p_ctrl4 + depth_n_ctrl4 + 1)
            vaf_ctrl_5 = float(var_p_ctrl5 + var_n_ctrl5) / (depth_p_ctrl5 + depth_n_ctrl5 + 1)
    
            max_vaf = max([vaf_ctrl_2, vaf_ctrl_3, vaf_ctrl_4, vaf_ctrl_5])
            median_vaf = statistics.median([vaf_ctrl_2, vaf_ctrl_3, vaf_ctrl_4, vaf_ctrl_5])
            """
            is_error = "False"
            error_num = 0
            if var_p_ctrl2 + var_n_ctrl2 >= 2: error_num = error_num + 1
            if var_p_ctrl3 + var_n_ctrl3 >= 2: error_num = error_num + 1
            if var_p_ctrl4 + var_n_ctrl4 >= 2: error_num = error_num + 1
            if var_p_ctrl5 + var_n_ctrl5 >= 2: error_num = error_num + 1
            if error_num >= 3: is_error = "True"
            """

            print('\t'.join(var_info) + '\t' + ctrl_info + '\t' + str(round(max_vaf, 4)) + '\t' + str(round(median_vaf, 4)), file = hout)


    hout.close()


if __name__ == "__main__":

    import sys
    var_file = sys.argv[1]
    output_file = sys.argv[2]
    control_files = sys.argv[3:]
    # proc_mpileup(input_file, output_file, min_variant_num_tumor = 5, min_variant_ratio_tumor = 0.15)
    short_read_validate(var_file, control_file, output_file)


