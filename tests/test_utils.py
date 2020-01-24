#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import nanomonsnv
from nanomonsnv import utils as utils


class Test_utils_main(unittest.TestCase):

    def test1(self):
        ref = "A"
        bases = "^gtTTtT..-3CTG*t,+1c,-2ctT,TT*.tt.tCTTTtt*TT*tTTt..T,T,+1c,.T,,,..,..t*,,,.*,,,.T.T$"
        base_qualities = ",50)6.&+(/%;*5-,(/01)0/53,(%*+($28(%,-4-'2#4+&,)/.',*-28$%'0.'#256"
        var2num, ovar2bq, depth = utils.check_pileup_record(ref, bases, base_qualities)

        print(var2num)
        print(ovar2bq)
        print(depth)

        ans_var2num = {'A': 30, 'C': 1, 'G': 0, 'T': 29, 'N': 0}
        ans_ovar2bq = {'A': [13, 5, 7, 16, 4, 11, 2, 8, 14, 6, 11, 4, 2, 20], 
                        'C': [15], 
                        'G': [], 
                        'T': [20, 15, 21, 26, 20, 12, 14, 20, 18, 9, 10, 17, 23, 12, 12, 19, 17, 21],
                        'N': [],
                        'a': [14, 4, 9, 19, 6, 17, 10, 5, 11, 13, 17, 23, 3, 15, 13, 6],
                        'c': [],
                        'g': [],
                        't': [11, 8, 7, 14, 15, 8, 11, 7, 3, 7, 9],
                        'n': []}
        ans_depth = 60

        self.assertEqual(ans_var2num, var2num)
        self.assertEqual(ans_ovar2bq, ovar2bq)
        self.assertEqual(ans_depth, depth)
