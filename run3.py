#! /usr/bin/env python

import sys, math, subprocess, concurrent.futures, argparse
import scipy.stats as stats

from utils import check_pileup_record


def short_read_validate(pileup_line, var_info, output_file_handle):

    F = pileup_line.rstrip('\n').split('\t')

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
 
    print('\t'.join(var_info) + '\t' + \
          str(depth_tumor) + '\t' + str(var_num_tumor) + '\t' + str(vaf_tumor) + '\t' + \
          str(depth_ctrl) + '\t' + str(var_num_ctrl) + '\t' + str(vaf_ctrl) + '\t' + \
          str(round(-math.log10(pvalue), 4)) + '\t' + is_validated, file = output_file_handle)          


def add_validation_info_region(var_file, tumor_bam, control_bam, reference_genome, output_file, region, samtools_options = "-Q 0 -q 40 -B"):

    target_rname = region.split(':')[0]
    chr_pos2info = {}
    with open(var_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0] == target_rname:
                chr_pos2info[F[0] + '\t' +  F[1]] = '\t'.join(F)


    hout = open(output_file, 'w')
    samtools_option_list = samtools_options.split(' ')
    samtools_cmd = ["samtools", "mpileup", tumor_bam, control_bam, "-f", reference_genome, "-r", region] + samtools_option_list
    print(' '.join(samtools_cmd))
    proc = subprocess.Popen(samtools_cmd, stdout = subprocess.PIPE) # , stderr = subprocess.DEVNULL)

    for pileup_line in proc.stdout:
        F = pileup_line.decode().rstrip('\n').split('\t')
        if F[0] + '\t' + F[1] not in chr_pos2info: continue
        var_info = chr_pos2info[F[0] + '\t' + F[1]].split('\t')
        short_read_validate(pileup_line.decode(), var_info, hout)

    hout.close()


def main(var_file, output_file, tumor_bam, control_bam, reference_genome):

    import pysam

    tbamfile = pysam.AlignmentFile(tumor_bam)
    rname2seqlen = {}
    for i in range(tbamfile.nreferences):
        rname = tbamfile.getrname(i)
        seq_length = tbamfile.get_reference_length(rname)
        rname2seqlen[rname] = seq_length

    executor = concurrent.futures.ProcessPoolExecutor(max_workers = 8)
    futures = []
    for rname in rname2seqlen:
        region = rname + ":1-" + str(rname2seqlen[rname])
        future = executor.submit(add_validation_info_region, var_file, tumor_bam, control_bam, reference_genome, output_file + ".tmp." + rname, region)
        futures.append(future)

    for future in concurrent.futures.as_completed(futures):
        if future.exception() is not None or future.result() is not None:
            print(future.exception())
            print(future.result())
            executor.shutdown()
            sys.exit(1)


    hout = open(output_file, 'w')
    for rname in rname2seqlen:
        with open(output_file + ".tmp." + rname) as hin:
            for line in hin:
                print(line.rstrip('\n'), file = hout)


if __name__ == "__main__":

    """
    var_file = "out2.txt"
    output_file = "out3.txt"
    tumor_bam = "s3://eva-bucket-tokyo/kataoka-lab/wgs/cell-line/rawdata/COLO829/COLO829.markdup.bam"
    control_bam = "s3://eva-bucket-tokyo/kataoka-lab/wgs/cell-line/rawdata/COLO829BL/COLO829BL.markdup.bam"
    reference_genome = "/home/ubuntu/environment/workspace/seq_data/reference/GRCh37.fa"
    # proc_mpileup(input_file, output_file, min_variant_num_tumor = 5, min_variant_ratio_tumor = 0.15)
    # short_read_validate(var_file, illumina_file, output_file)
    """


    irser = argparse.ArgumentParser(description = "Validation using Illumina short read tumor and control sequencing data")
    parser.add_argument("variant_file", type = str, help = "Path to variant candidate file")
    parser.add_argument("tumor_bam", type = str, help = "Path to Illumina tumor bam file")
    parser.add_argument("control_bam", type = str, help = "Path to Illumina matched control file")
    parser.add_argument("output_file", type = str, help = "Path to the output file")
    parser.add_argument("reference", type = str, help = "Path to the reference genome")

    args = parser.parse_args()

    main(args.variant_file, args.output_file, args.tumor_bam, args.control_bam, args.reference_genome)


