from __future__ import print_function 
import sys, argparse
import os, time
import pysam
import subprocess
import gzip


def mnv_to_snv(vcf_database):
    # create UCSC to GRC chr name corresponding table

    header_end_flag = False
    # hout = open(output_file, 'w')
    with gzip.open(vcf_database, 'rt') as hin:
    # with open(vcf_database, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            F = line.split('\t')
            
            if line.startswith('#'):
                print(line)
            else:
                header_end_flag = True
            
            if not header_end_flag: continue
            
            TYPE = ""
            chrom = F[0]
            pos = F[1]
            ref = F[3]
            alt = F[4]
            infos = F[7].split(';')
            for info in infos:
                if info.startswith("TYPE="):
                    TYPE = info.replace("TYPE=", '')
                    break

            if TYPE == "MNV":
                # print("MNV")
                # print("pos:"+pos)
                l_ref = list(ref)
                for i, x in enumerate(l_ref):
                    print("\t".join([chrom, str(int(pos)+i), F[2], x, alt[i]]) +"\t"+ "\t".join(F[5:]))
                    
            # else
            elif TYPE == "SNV":
                # print("SNV")
                print(line)
                        

_vcf_database = sys.argv[1]
mnv_to_snv(_vcf_database)

