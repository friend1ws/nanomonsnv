from __future__ import print_function
import subprocess
import sys, re
import pysam

# def check_pileup_record(ref, bases, base_qualities, mapping_qualities):
def check_pileup_record(bases, base_qualities, qnames):

    l_qname = qnames.split(',')
    ovar2qname = {'A': [], 'C': [], 'G': [], 'T': [], 'N': [], 'a': [], 'c': [], 'g': [], 't': [], 'n': []}
    base_ind = 0
    depth = 0
    check_count = 0
    while bases != '':
        if bases[0] in ['>', '<', '*']: 
            base_ind = base_ind + 1
            bases = bases[1:]
        elif bases[0] in '^':
            bases = bases[2:]
        elif bases[0] in '$':
            bases = bases[1:]
        elif bases[0] in ['A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n']:
            var_original = bases[0]

            ovar2qname[var_original].append(l_qname[base_ind])
            
            depth = depth + 1     
            base_ind = base_ind + 1
            bases = bases[1:]

        if len(bases) > 0 and bases[0] in ['+', '-']:
            match = re.search(r'^[\+\-](\d+)', bases)
            indel_size = int(match.group(1))
            bases = bases[(len(str(indel_size)) + indel_size + 1):]

    if len(base_qualities) != base_ind:
        print("Error???")
        sys.exit(1)

    # return([var2num, ovar2bq, ovar2mq, depth])
    return ovar2qname

def proc_pileup_line(chrom, pos, ref, alt, out_line, pileup_line, samfile, output_file_handle):

    
    F = pileup_line.rstrip('\n').split('\t')

    # Pileup format
    # chrom, pos, referece base, the number of reads, read bases, base qualities, mapping qualities,qname
    # [0],   [1], [2],           [3],                 [4],        [5],            [6],
    # seq2 156 A 11  .$......+2AG.+2AG.+2AGGG <975;:<<<<< [[[[[[[[[[[ ID1,ID2,ID3,ID4,ID5
    #
    # chrom, pos, referece base, the number of reads, read bases, base qualities, qname
    # [0],   [1], [2],           [3],                 [4],        [5],            [6]

    # var2num_tumor, ovar2bq_tumor, over2mq_tumor, depth_tumor = check_pileup_record(F[2], F[4], F[5], F[6])
    ovar2qname = check_pileup_record(F[4], F[5], F[6])
    
    h_haplotag_plus = {}
    h_haplotag_minus = {}
    # print(chrom)
    # print(pos)
    iter = samfile.fetch(chrom, int(pos)-1, int(pos))
    for read in iter:
        tag = ""
        try:
            tag = read.get_tag("HP")
        except KeyError:
            tag = "---"
        except:
            raise ValueError
            
        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]
        
        if flags[4] == "0":
            h_haplotag_plus[read.qname] = str(tag)
        else:
            h_haplotag_minus[read.qname] = str(tag)
    
    h_haplotag_ref = {"1":0,"2":0,"---":0}
    h_haplotag_alt = {"1":0,"2":0,"---":0}
    
    for x in ovar2qname[ref]:
        h_haplotag_ref[h_haplotag_plus[x]] += 1
    for x in ovar2qname[ref.lower()]:
        h_haplotag_ref[h_haplotag_minus[x]] += 1
    for x in ovar2qname[alt]:
        h_haplotag_alt[h_haplotag_plus[x]] += 1   
    for x in ovar2qname[alt.lower()]:
        h_haplotag_alt[h_haplotag_minus[x]] += 1   
        
    # print(ovar2qname)
    # print(h_haplotag_plus)
    # print(h_haplotag_minus)
    # print(h_haplotag_ref)
    # print(h_haplotag_alt)

    print (out_line 
    +"\t"+ str(h_haplotag_ref["1"]) 
    +"\t"+ str(h_haplotag_ref["2"]) 
    +"\t"+ str(h_haplotag_ref["---"]) 
    +"\t"+ str(h_haplotag_alt["1"]) 
    +"\t"+ str(h_haplotag_alt["2"]) 
    +"\t"+ str(h_haplotag_alt["---"]), file = output_file_handle)
    
    
def pileup_main(input_file, tumor_bam, output_file):
    
    samtools_options = "-Q 0 -q 40 -B --output-QNAME"
    samtools_option_list = samtools_options.split(' ')
    pysam_bam = pysam.AlignmentFile(tumor_bam, "rb")
    hout = open(output_file, 'w')

    with open(input_file, "r") as hin:
        for line in hin:
            line = line.rstrip('\n')
            F = line.split('\t')
            
            chrom = F[0]
            pos = F[1]
            ref = F[2]
            alt = F[3]
            region = chrom+":"+pos+"-"+pos
    
            samtools_cmd = ["samtools", "mpileup", tumor_bam, "-r", region] + samtools_option_list
            print(' '.join(samtools_cmd))
    
    
            with subprocess.Popen(samtools_cmd, stdout = subprocess.PIPE) as proc:
                try:
                    for pileup_line in proc.stdout:
                        proc_pileup_line(chrom, pos, ref, alt, line, pileup_line.decode(), pysam_bam, hout)    
                except:
                    raise ValueError
        
    pysam_bam.close()
    hout.close()        

if __name__ == "__main__":
    tumor_bam = "/home/ubuntu/bam/COLO829_phased2.chr22.bam"
    input_file = "COLO829BL_vallidateout_filter_chr22.txt"
    output_file = "COLO829BL_vallidateout_filter_chr22_haplotag.txt"
    
    tumor_bam = sys.argv[1]
    input_file = sys.argv[2]
    output_file = sys.argv[3]
    
    pileup_main(input_file, tumor_bam, output_file)

