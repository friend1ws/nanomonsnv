#! /usr/bin/env python

import re
import pysam

def check_pileup_record(ref, bases, qualities):

    var2num = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
    ovar2qual = {'A': [], 'C': [], 'G': [], 'T': [], 'N': [], 'a': [], 'c': [], 'g': [], 't': [], 'n': []}
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
        elif bases[0] in ['.', ',', 'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n']:
            if bases[0] not in ['.', ',']: 
                var_original = bases[0]
            else:
                var_original = ref.upper() if bases[0] == '.' else ref.lower()

            ovar2qual[var_original].append(ord(qualities[base_ind]) - 33)
            var = var_original.upper()
            var2num[var] = var2num[var] + 1
            depth = depth + 1     
            base_ind = base_ind + 1
            bases = bases[1:]

        if len(bases) > 0 and bases[0] in ['+', '-']:
            match = re.search(r'^[\+\-](\d+)', bases)
            indel_size = int(match.group(1))
            bases = bases[(len(str(indel_size)) + indel_size + 1):]

    if len(qualities) != base_ind:
        print("Error???")
        sys.exit(1)

    return([var2num, ovar2qual, depth])



def get_seq(reference, chr, start, end):

    seq = ""
    for item in pysam.faidx(reference, chr + ":" + str(start) + "-" + str(end)):
        # if item[0] == ">": continue
        seq = seq + item.rstrip('\n')
    seq = seq.replace('>', '')
    seq = seq.replace(chr + ":" + str(start) + "-" + str(end), '')

    if re.search(r'[^ACGTUWSMKRYBDHVNacgtuwsmkrybdhvn]', seq) is not None:
        print("The return value in get_seq function includes non-nucleotide characters:", file = sys.stderr)
        print(seq, file = sys.stderr)
        sys.exit(1)

    return seq

