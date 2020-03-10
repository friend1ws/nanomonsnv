from __future__ import print_function 
import sys, argparse
import os, time
import math
import subprocess
import pysam
import shutil

# TP / (TP + FN) = recall
def tpr(lvT, relevant_element):
    return len(lvT) / relevant_element

# TP / (TP + FP) = precision
def ppv(lvT, lvF):
    return len(lvT) / (len(lvT)+len(lvF))
    
def get_presigion_recall(variant_file, group, hout, anno_idx, relevant_element=42993):

    l_valid_T = []
    l_valid_F = []
    with open(variant_file, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            F = line.split('\t')
            if F[anno_idx] == "True":
                l_valid_T.append(1)
            else:
                l_valid_F.append(1)
            
            print (str(round(ppv(l_valid_T,l_valid_F),4)) +"\t"+ str(round(tpr(l_valid_T, relevant_element),4)) +"\t"+ group, file = hout)


def filter_result(in_file, out_file, filter_type):
    
    hout = open(out_file, 'w')
    with open(in_file, 'r') as hin:
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
            
            print('\t'.join(F), file = hout)
    hout.close()

def annotate_anno(variant_file, output_file, ccl_tabix, chrom_idx, pos_idx, ref_idx, alt_idx):

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

            if line.startswith('#'): continue
        
            line = line.rstrip('\n')
            F = line.split('\t')
            nanosnv_record = line
            
            chrom = F[chrom_idx]
            pos = F[pos_idx]
            ref = F[ref_idx]
            alt = F[alt_idx]
            
            if ref == "-" or alt == "-": continue
            if chrom not in target_rnames_ccl: continue
                
            called_by = "False"
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

                called_by = "True"
                break

            nanosnv_record = nanosnv_record + "\t" + called_by
            
            print(nanosnv_record, file = hout)  

    hout.close()            


if __name__ == "__main__":

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    ccl_tabix = sys.argv[3]
    illumina_file = sys.argv[4]
    illumina_file2 = sys.argv[5]
    
    cmd = ["sort", "-k", "11", "-r", "-n", input_file]
    with open(output_file+".sorted.txt", 'w') as hout:
        subprocess.check_call(cmd, stdout=hout)
    
    cmd = ["sort", "-k", "18", "-r", "-n", illumina_file]
    with open(output_file+".sorted.illumina.txt", 'w') as hout:
        subprocess.check_call(cmd, stdout=hout)
    
    if illumina_file2 != "None":
        cmd = ["sort", "-k", "18", "-r", "-n", illumina_file2]
        with open(output_file+".sorted.illumina2.txt", 'w') as hout:
            subprocess.check_call(cmd, stdout=hout)
    
    annotate_anno(output_file+".sorted.txt", output_file+".anno.txt", ccl_tabix, 0,1,2,3)
    annotate_anno(output_file+".sorted.illumina.txt", output_file+".anno.illumina.txt", ccl_tabix, 0,1,3,4)
    if illumina_file2 != "None":
        annotate_anno(output_file+".sorted.illumina2.txt", output_file+".anno.illumina2.txt", ccl_tabix, 0,1,3,4)
        
    shutil.copy(output_file+".anno.txt", output_file+".filt1.txt")
    filter_result(output_file+".anno.txt", output_file+".filt2.txt", 2)
    filter_result(output_file+".anno.txt", output_file+".filt3.txt", 3)
    
    with open(output_file, 'w') as hout:
        get_presigion_recall(output_file+".filt1.txt","Method1", hout, 38)
        get_presigion_recall(output_file+".filt2.txt","Method2", hout, 38)
        get_presigion_recall(output_file+".filt3.txt","Method3", hout, 38)
        get_presigion_recall(output_file+".anno.illumina.txt","Illumina", hout, 31)
        if illumina_file2 != "None":
            get_presigion_recall(output_file+".anno.illumina2.txt","Illumina(subsampled)", hout, 31)
    
    os.remove(output_file+".sorted.txt")
    os.remove(output_file+".sorted.illumina.txt")
    os.remove(output_file+".anno.txt")
    os.remove(output_file+".anno.illumina.txt")
    os.remove(output_file+".filt1.txt")
    os.remove(output_file+".filt2.txt")
    os.remove(output_file+".filt3.txt")
    if illumina_file2 != "None":
        os.remove(output_file+".sorted.illumina2.txt")
        os.remove(output_file+".anno.illumina2.txt")
