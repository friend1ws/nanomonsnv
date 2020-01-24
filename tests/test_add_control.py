#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import nanomonsnv


class Test_add_control_main(unittest.TestCase):

    def setUp(self):
        self.parser = nanomonsnv.parser.create_parser()

    def test1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        variant_file =cur_dir + "/../data/COLO829_V600_test_detectout.txt"
        ctrl_bam1 = cur_dir + "/../data/BL2009_v600_test.bam"
        ctrl_bam2 = cur_dir + "/../data/HCC1954BL_v600_test.bam"
        reference = "/home/ubuntu/environment/database/GRCh37.fa"
        output_file = tmp_dir+"/COLO829_V600_test1_addctrlout.txt"
        args = self.parser.parse_args(["add_control", variant_file, output_file, reference, ctrl_bam1, ctrl_bam2])
        args.func(args)

        answer_file = cur_dir + "/../data/COLO829_V600_test_addctrlout.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)


