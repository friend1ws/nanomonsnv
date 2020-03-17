import unittest
import pileup as pu
import os, tempfile, shutil, filecmp

class TestPileup(unittest.TestCase):
    
    def setUp(self):
        pass
    
    def test_check_pileup_record1(self):
        bases = "ACGT"
        base_qualities = "1234"
        qnames = 'e06069c3-26e3-48c7-95bc-4344dcc2b16d,' \
        +'b1228fe0-4692-41f6-a173-601fb443a4af,' \
        +'cab6765a-7e46-4ceb-941b-714c02a126a6,' \
        +'964fcfaf-9694-4a8e-a585-73b60c6cd39d'

        ovar2qname = pu.check_pileup_record(bases, base_qualities, qnames)

        l_result = {"A":['e06069c3-26e3-48c7-95bc-4344dcc2b16d'],
        "C":['b1228fe0-4692-41f6-a173-601fb443a4af'],
        "G":['cab6765a-7e46-4ceb-941b-714c02a126a6'],
        "T":['964fcfaf-9694-4a8e-a585-73b60c6cd39d'],
        "N":[],"a":[],"c":[],"g":[],"t":[],"n":[]
        }

    def test_check_pileup_record2(self):
        bases = "AACC"
        base_qualities = "1234"
        qnames = 'e06069c3-26e3-48c7-95bc-4344dcc2b16d,'\
        +'b1228fe0-4692-41f6-a173-601fb443a4af,' \
        +'cab6765a-7e46-4ceb-941b-714c02a126a6,' \
        +'964fcfaf-9694-4a8e-a585-73b60c6cd39d'

        ovar2qname = pu.check_pileup_record(bases, base_qualities, qnames)

        l_result = {"A":['e06069c3-26e3-48c7-95bc-4344dcc2b16d','b1228fe0-4692-41f6-a173-601fb443a4af'],
        "C":['cab6765a-7e46-4ceb-941b-714c02a126a6','964fcfaf-9694-4a8e-a585-73b60c6cd39d'],
        "G":[],
        "T":[],
        "N":[],"a":[],"c":[],"g":[],"t":[],"n":[]
        }
        
        assert all( (k,v) in l_result.items() for (k,v) in ovar2qname.items() )

    def test_check_pileup_record3(self):
        bases = "A*A*C*C"
        base_qualities = "1x2x3x4"
        qnames = 'e06069c3-26e3-48c7-95bc-4344dcc2b16d,' \
        +'hogehoge,' \
        +'b1228fe0-4692-41f6-a173-601fb443a4af,' \
        +'hogehoge,' \
        +'cab6765a-7e46-4ceb-941b-714c02a126a6,' \
        +'hogehoge,' \
        +'964fcfaf-9694-4a8e-a585-73b60c6cd39d'

        ovar2qname = pu.check_pileup_record(bases, base_qualities, qnames)

        l_result = {"A":['e06069c3-26e3-48c7-95bc-4344dcc2b16d','b1228fe0-4692-41f6-a173-601fb443a4af'],
        "C":['cab6765a-7e46-4ceb-941b-714c02a126a6','964fcfaf-9694-4a8e-a585-73b60c6cd39d'],
        "G":[],
        "T":[],
        "N":[],"a":[],"c":[],"g":[],"t":[],"n":[]
        }
        
        assert all( (k,v) in l_result.items() for (k,v) in ovar2qname.items() )
        

    def test_pileup_main4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        # input_bam = "/home/ubuntu/bam/COLO829_phased2.chr22.bam"
        input_bam = cur_dir + "/data/phased_test4.bam"
        input_txt = cur_dir + "/data/phased_test4.txt"
        output_file = tmp_dir+"/phased_test4.txt"
        pu.pileup_main(input_txt, input_bam, output_file)

        answer_file = cur_dir + "/data/phased_test4.answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))


    def test_pileup_main5(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        # input_bam = "/home/ubuntu/bam/COLO829_phased2.chr22.bam"
        input_bam = cur_dir + "/data/phased_test5.bam"
        input_txt = cur_dir + "/data/phased_test5.txt"
        output_file = tmp_dir+"/phased_test5.txt"
        pu.pileup_main(input_txt, input_bam, output_file)

        answer_file = cur_dir + "/data/phased_test5.answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        
        
if __name__ == "__main__":
    unittest.main()