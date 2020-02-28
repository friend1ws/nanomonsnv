from __future__ import print_function 
import sys, argparse
import os, time
import pysam
import subprocess

def annotate_anno(variant_file, output_file, ccl_tabix):

    target_rnames_ccl = []
    tabix_cmd = ['tabix','-l',ccl_tabix]
    print(tabix_cmd)
    proc = subprocess.Popen(tabix_cmd, stdout = subprocess.PIPE)
    for chromosome_name in proc.stdout:
        chromosome_name = chromosome_name.decode().rstrip('\n')
        target_rnames_ccl.append(chromosome_name)
    proc.stdout.close()
    proc.wait()
    
    header_end_flag = False
    ccl_tabix_db = pysam.Tabixfile(ccl_tabix, encoding="utf-8")
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
            # chrom = "chr"+F[0]
            chrom = F[0]
            pos = F[1]
            ref = F[2]
            alt = F[3]
            if ccl_tabix is not None:
            
                called_by = False
                if chrom in target_rnames_ccl:
                    for record_line in ccl_tabix_db.fetch(chrom, int(pos) - 1, int(pos) + 1):
                        record = record_line.split('\t')
                        
                        TYPE = ""
                        continue_flag = True
                        infos = record[7].split(';')
                        for info in infos:
                            if info.startswith("TYPE="):
                                TYPE = info.replace("TYPE=", '')
                                if TYPE in ["SNV","MNV"]:
                                    continue_flag = False
                                    break
                        
                        if continue_flag: continue
                        if record[0] != chrom: continue
                        if record[1] != pos: continue
                        if record[4] != alt: continue

                        infos = record[7].split(';')
                        for info in infos:
                            if info.startswith("called_by="):
                                # called_by = info.replace("called_by=", '')
                                called_by = True
                                break

                nanosnv_record = nanosnv_record + "\t" + str(called_by)
                
            print(nanosnv_record, file = hout)  
            
            
variant_file = sys.argv[1]
output_file = sys.argv[2]
ccl_tabix = sys.argv[3]

annotate_anno(variant_file, output_file, ccl_tabix)

