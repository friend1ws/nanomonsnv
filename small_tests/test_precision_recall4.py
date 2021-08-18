import unittest
import precision_recall4 as pr
import os, tempfile, shutil, filecmp

class TestPrecisionRecall4(unittest.TestCase):
    
    def setUp(self):
        pass
    
    def test_annotate_anno1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_mutation = cur_dir + "/data/variant_file1.txt"
        output_file = tmp_dir+"/variant_file1.txt"
        ccl_tabix = cur_dir + "/../../../comp/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.snv_mnv.vcf.gz"
        pr.annotate_anno(in_mutation, output_file, ccl_tabix, 0,1,2,3)

        answer_file = cur_dir + "/data/variant_file1.answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test_annotate_anno2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_mutation = cur_dir + "/data/variant_file2.txt"
        output_file = tmp_dir+"/variant_file2.txt"
        ccl_tabix = cur_dir + "/../../../comp/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.snv_mnv.vcf.gz"
        pr.annotate_anno(in_mutation, output_file, ccl_tabix, 0,1,2,3)

        answer_file = cur_dir + "/data/variant_file2.answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)

    def test_annotate_anno3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_mutation = cur_dir + "/data/variant_file3.txt"
        output_file = tmp_dir+"/variant_file3.txt"
        ccl_tabix = cur_dir + "/../../../comp/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.snv_mnv.vcf.gz"
        pr.annotate_anno(in_mutation, output_file, ccl_tabix, 0,1,2,3)

        answer_file = cur_dir + "/data/variant_file3.answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test_annotate_anno4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_mutation = cur_dir + "/data/variant_file4.txt"
        output_file = tmp_dir+"/variant_file4.txt"
        ccl_tabix = cur_dir + "/../../../comp/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.snv_mnv.vcf.gz"
        pr.annotate_anno(in_mutation, output_file, ccl_tabix, 0,1,2,3)

        answer_file = cur_dir + "/data/variant_file4.answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test_annotate_anno5(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_mutation = cur_dir + "/data/variant_file5.txt"
        output_file = tmp_dir+"/variant_file5.txt"
        ccl_tabix = cur_dir + "/../../../comp/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.snv_mnv.vcf.gz"
        pr.annotate_anno(in_mutation, output_file, ccl_tabix, 0,1,2,3)

        answer_file = cur_dir + "/data/variant_file5.answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test_annotate_anno6(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_mutation = cur_dir + "/data/variant_file6.txt"
        output_file = tmp_dir+"/variant_file6.txt"
        ccl_tabix = cur_dir + "/../../../comp/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.snv_mnv.vcf.gz"
        pr.annotate_anno(in_mutation, output_file, ccl_tabix, 0,1,2,3)

        answer_file = cur_dir + "/data/variant_file6.answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        

    def test_filter_result7(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_mutation = cur_dir + "/data/variant_file7.txt"
        output_file = tmp_dir+"/variant_file7.txt"
        pr.filter_result(in_mutation, output_file, 2)

        answer_file = cur_dir + "/data/variant_file7.answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        # shutil.rmtree(tmp_dir)
        
    def test_filter_result8(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_mutation = cur_dir + "/data/variant_file8.txt"
        output_file = tmp_dir+"/variant_file8.txt"
        pr.filter_result(in_mutation, output_file, 3)

        answer_file = cur_dir + "/data/variant_file8.answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        # shutil.rmtree(tmp_dir)

    def test_filter_result9(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_mutation = cur_dir + "/data/variant_file9.txt"
        output_file = tmp_dir+"/variant_file9.txt"
        pr.filter_result(in_mutation, output_file, 4)

        answer_file = cur_dir + "/data/variant_file9.answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        # shutil.rmtree(tmp_dir)
        
    def test_presigion_recall_result10(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_mutation = cur_dir + "/data/variant_file10.txt"
        output_file = cur_dir+"/variant_file10.txt"
        hout = open(output_file, 'w')
        pr.get_presigion_recall(in_mutation, "Test", hout, 0, 20)
        hout.close()

        answer_file = cur_dir + "/data/variant_file10.answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        # shutil.rmtree(tmp_dir)
        
if __name__ == "__main__":
    unittest.main()