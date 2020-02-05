#! /usr/bin/env python

import sys, argparse
import os, time
import pysam

target_rnames = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]


def annotate_anno(variant_file, output_file, gnomad_genome, simple_repeat):

    if gnomad_genome is not None: gnomad_genome_db = pysam.Tabixfile(gnomad_genome, encoding="utf-8")
    if simple_repeat is not None: simple_repeat_db = pysam.Tabixfile(simple_repeat, encoding="utf-8")

    header_end_flag = False
    info_start_flag = False
    info_end_flag = False

    hout = open(output_file, 'w')
    with open(variant_file, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')

            if line.startswith('#'):
                print(line, file = hout)
            else:
                header_end_flag = True

            if not header_end_flag: continue

            nanosnv_record = line
            F = line.split('\t')
            chrom = F[0].replace("chr", '')
            if gnomad_genome is not None:

                GNOMAD_GENOME = 0
                if chrom in target_rnames + ["chr" + x for x in target_rnames]:
                    for record_line in gnomad_genome_db.fetch(chrom, int(F[1]) - 1, int(F[1]) + 1):
                        record = record_line.split('\t')
                        if record[0] != chrom: continue
                        if record[1] != F[1]: continue
                        if record[3] != F[2]: continue
                        if record[4] != F[3]: continue
                        infos = record[7].split(';')
                        for info in infos:
                            if info.startswith("AF="):
                                GNOMAD_GENOME = float(info.replace("AF=", ''))
                nanosnv_record = nanosnv_record + "\t" + str(round(GNOMAD_GENOME, 4))

            if simple_repeat is not None:
            
                repeat_seq = []
                if chrom in target_rnames + ["chr" + x for x in target_rnames]:
                    for record_line in simple_repeat_db.fetch("chr"+chrom, int(F[1]) - 1, int(F[1]) + 1):
                        record = record_line.split('\t')
                        if record[0] != "chr"+chrom: continue
                        if record[1] >= F[1]: continue
                        if record[2] < F[1]: continue
                        repeat_seq.append(record[3])

                nanosnv_record = nanosnv_record + "\t" + ",".join(repeat_seq)

            print(nanosnv_record, file = hout)

    hout.close()


def filter_main(args):

    annotate_anno(args.variant_file, args.output_file, args.gnomad_file, args.simple_repeat_file)

