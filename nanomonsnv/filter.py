#! /usr/bin/env python

import sys, argparse
import os, time
import pysam
import subprocess

def annotate_anno(variant_file, output_file, gnomad_genome, simple_repeat):

    target_rnames_gnomad = []
    if gnomad_genome is not None:
        gnomad_genome_db = pysam.Tabixfile(gnomad_genome, encoding="utf-8")

        tabix_cmd = ['tabix','-l',gnomad_genome]
        with subprocess.Popen(tabix_cmd, stdout = subprocess.PIPE) as proc:
            for chromosome_name in proc.stdout:
                chromosome_name = chromosome_name.decode().rstrip('\n')
                target_rnames_gnomad.append(chromosome_name)

    target_rnames_simple_repeat = []
    if simple_repeat is not None:
        simple_repeat_db = pysam.Tabixfile(simple_repeat, encoding="utf-8")

        tabix_cmd = ['tabix','-l',simple_repeat]
        with subprocess.Popen(tabix_cmd, stdout = subprocess.PIPE) as proc:
            for chromosome_name in proc.stdout:
                chromosome_name = chromosome_name.decode().rstrip('\n')
                target_rnames_simple_repeat.append(chromosome_name)
                

    header_end_flag = False
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
                if chrom in target_rnames_gnomad:
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
                                break
                nanosnv_record = nanosnv_record + "\t" + str(round(GNOMAD_GENOME, 4))
                

            if simple_repeat is not None:
            
                repeat_seq = []
                if "chr" + chrom in target_rnames_simple_repeat:
                    for record_line in simple_repeat_db.fetch("chr"+chrom, int(F[1]) - 1, int(F[1]) + 1):
                        record = record_line.split('\t')
                        if record[0] != "chr"+chrom: continue
                        if record[1] >= F[1]: continue
                        if record[2] < F[1]: continue
                        repeat_seq.append(record[15])

                nanosnv_record = nanosnv_record + "\t" + ",".join(repeat_seq)

            print(nanosnv_record, file = hout)

    hout.close()


def filter_main(args):

    annotate_anno(args.variant_file, args.output_file, args.gnomad_file, args.simple_repeat_file)

