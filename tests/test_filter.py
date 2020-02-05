#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import nanomonsnv


class Test_filter_main(unittest.TestCase):

    def setUp(self):
        self.parser = nanomonsnv.parser.create_parser()

    def test1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        variant_file =cur_dir + "/../data/COLO829_V600_test_validateout.txt"
        output_file = tmp_dir+"/COLO829_V600_test1_filterout.txt"
        gnomad_file =cur_dir + "/../data/V600_AF.vcf.gz"
        simple_repeat_file =cur_dir + "/../data/V600_SR.bed.gz"
        args = self.parser.parse_args(["filter", variant_file, output_file, gnomad_file, simple_repeat_file])
        args.func(args)

        answer_file = cur_dir + "/../data/COLO829_V600_test_filterout.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)

