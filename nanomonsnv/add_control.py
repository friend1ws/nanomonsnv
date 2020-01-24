#! /usr/bin/env python

import sys, re, argparse, statistics, subprocess, concurrent.futures
import scipy.stats as stats
import os, time

from .utils import check_pileup_record


def proc_ctrl_pileup_record_info(ref, alt, bases, qualities):

    var2num_ctrl, ovar2qual_ctrl, depth_ctrl = check_pileup_record(ref, bases, qualities)
    depth_p_ctrl, depth_n_ctrl = 0, 0
    for var in ['A', 'C', 'G', 'T', 'N']:
        depth_p_ctrl = depth_p_ctrl + len(ovar2qual_ctrl[var])
        depth_n_ctrl = depth_n_ctrl + len(ovar2qual_ctrl[var.lower()])
   
    var_p_ctrl, var_n_ctrl = len(ovar2qual_ctrl[alt.upper()]), len(ovar2qual_ctrl[alt.lower()])

    return([depth_p_ctrl, var_p_ctrl, depth_n_ctrl, var_n_ctrl])


def proc_pileup_line_eb(pileup_line, var_info, output_file_handle):

    # import pdb; pdb.set_trace()

    F = pileup_line.rstrip('\n').split('\t')
    ctrl_num = int((len(F) - 3) / 3)

    ctrl_infos = []
    vaf_ctrls = []
    for i in range(ctrl_num):
        depth_p_ctrl, var_p_ctrl, depth_n_ctrl, var_n_ctrl = proc_ctrl_pileup_record_info(F[2], var_info[3], F[3 * i + 4], F[3 * i + 5])
        ctrl_infos.append(','.join([str(depth_p_ctrl), str(var_p_ctrl), str(depth_n_ctrl), str(var_n_ctrl)]))
        vaf = float(var_p_ctrl + var_n_ctrl) / (depth_p_ctrl + depth_n_ctrl + 1)
        vaf_ctrls.append(vaf)

    max_vaf = max(vaf_ctrls)
    median_vaf = statistics.median(vaf_ctrls)
    ctrl_info = ';'.join(ctrl_infos)
 
    print('\t'.join(var_info) + '\t' + ctrl_info + '\t' + str(round(max_vaf, 4)) + '\t' + str(round(median_vaf, 4)), file = output_file_handle)



def add_control_info_region(var_file, control_bam, reference_genome, output_file, region, samtools_options = "-Q 0 -q 40 -B"):

    target_rname = region.split(':')[0]
    chr_pos2info = {}
    with open(var_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0] == target_rname:
                chr_pos2info[F[0] + '\t' +  F[1]] = '\t'.join(F)

    time.sleep(1)
    hout = open(output_file, 'w')
    samtools_option_list = samtools_options.split(' ')
    samtools_cmd = ["samtools", "mpileup"] + control_bam + ["-f", reference_genome, "-r", region] + samtools_option_list
    print(' '.join(samtools_cmd))
    with subprocess.Popen(samtools_cmd, stdout = subprocess.PIPE) as proc: # , stderr = subprocess.DEVNULL): 

        """
        if proc.returncode != 0:
            # return(Exception) 
            raise ValueError
        """
    
        for pileup_line in proc.stdout:
            F = pileup_line.decode().rstrip('\n').split('\t')
            if F[0] + '\t' + F[1] not in chr_pos2info: continue
            var_info = chr_pos2info[F[0] + '\t' + F[1]].split('\t')
            proc_pileup_line_eb(pileup_line.decode(), var_info, hout)

    hout.close()



def add_control_main(args):

    import pysam

    tbamfile = pysam.AlignmentFile(args.control_bams[0])
    rname2seqlen = {}
    for i in range(tbamfile.nreferences):
        rname = tbamfile.getrname(i)
        seq_length = tbamfile.get_reference_length(rname)
        rname2seqlen[rname] = seq_length

    if args.max_workers > 1:
        executor = concurrent.futures.ProcessPoolExecutor(max_workers = args.max_workers)
        futures = []
        for rname in rname2seqlen:
            region = rname + ":1-" + str(rname2seqlen[rname])
            future = executor.submit(add_control_info_region, args.variant_file, args.control_bams, args.reference, args.output_file + ".tmp." + rname, region)
            futures.append(future)

        for future in concurrent.futures.as_completed(futures):
            """
            print("print future.exception()")
            print(future.exception())
            print()
            print("print future.result()")
            print(future.result())
            """
            if future.exception() is not None or future.result() is not None:
                print(future.exception())
                print(future.result())
                executor.shutdown()
                sys.exit(1)

    else:
    # if args.max_workers == 1
        for rname in rname2seqlen:
            region = rname + ":1-" + str(rname2seqlen[rname])
            add_control_info_region(args.variant_file, args.control_bams, args.reference, args.output_file + ".tmp." + rname, region)
        
    hout = open(args.output_file, 'w')
    for rname in rname2seqlen:
        with open(args.output_file + ".tmp." + rname) as hin:
            for line in hin:
                print(line.rstrip('\n'), file = hout)
        os.remove(args.output_file + ".tmp." + rname)


"""
if __name__ == "__main__":


    var_file = "out.txt"
    output_file = "out2.txt"
    control_bam = ["s3://eva-bucket-tokyo/kataoka-lab/long_read_sequencing/cell-line/minimap2/HCC1954.bam",
                   "s3://eva-bucket-tokyo/kataoka-lab/long_read_sequencing/cell-line/minimap2/HCC1954BL.bam",
                   "s3://eva-bucket-tokyo/kataoka-lab/long_read_sequencing/cell-line/minimap2/H2009.bam",
                   "s3://eva-bucket-tokyo/kataoka-lab/long_read_sequencing/cell-line/minimap2/BL2009.bam"]
    reference_genome = "/home/ubuntu/environment/workspace/seq_data/reference/GRCh37.fa"

    parser = argparse.ArgumentParser(description = "Adding non-matched control infomation to the variant candidate file")
    parser.add_argument("variant_file", type = str, help = "Path to variant candidate file")
    parser.add_argument("output_file", type = str, help = "Path to the output file")
    parser.add_argument("reference", type = str, help = "Path to the reference genome")
    parser.add_argument("control_bams", type = str, nargs='+', help = "Path to non-matched control files")
    args = parser.parse_args()

    main(args.variant_file, args.output_file, args.control_bams, args.reference_genome)
"""
