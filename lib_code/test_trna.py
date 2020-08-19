#  tRNAExplorer v.1.0 - tRNA profiles and regulation networks
#  Copyright (C) 2020  Dr. Quanyuan He
#  School of Medicine, Hunan Normal University
#  Email: hqyone@hotmail.com
#  Freely distributed under the GNU General Public License (GPLv3)
#

from unittest import TestCase
from trna import tRNA,getTRFType
tRNA_str='chr6	27605637	27605745	tRNA-Leu-CAA-2-1	0	-	GTCAGGATGGCCGAGTGGTCTAAGGCGCCAGACTCAAGCTTACTGCTTCCTGTGTTCGGGTCTTCTGGTCTCCGTATGGAGGCGTGGGTTCGAATCCCACTTCTGACA	GTCAGGATGGCCGAGTGGTCTAAGGCGCCAGACTCAAGTTCTGGTCTCCGTATGGAGGCGTGGGTTCGAATCCCACTTCTGACA	intron: 39-62 (99-122)	CAA	95	97	Leu	TF_tRNA-Leu-CAA-2-1	AAGGCGAGCCGTAAGCCGTGATTACTCATCATGTATAGTTTCAAGGATTTTTGTCACAGTGTCAGGATGGCCGAGTGGTCTAAGGCGCCAGACTCAAGCTTACTGCTTCCTGTGTTCGGGTCTTCTGGTCTCCGTATGGAGGCGTGGGTTCGAATCCCACTTCTGACACAGCCTGATCTTTTGTTTTTCCCTACCAGATTTGTACCGAAAGTTTCGCCAAGGAATGCC	60	61	168	108	>>>>>>>..>>>...........<<<.>>>>>...............................<<<<<.>>>>....<<<<..>>>>>.......<<<<<<<<<<<<.	GTCAGGATGGCCGAGTGGTctAAGGCGCCAGACTCAAGcttactgcttcctgtgttcgggtcTTCTGGTCTCCGTATGGAGGCGTGGGTTCGAATCCCACTTCTGACA	35	37	78.0		9	12	23	26	27	32	63	68	-1	-1	-1	-1	83	88	95	100	0	7	100	107'
mapping_infor='I,53,87,102.0,TGTAGGCGTGGGTTCGGATCCCACTTCTGACACCA,TGGAGGCGTGGGTTCGAATCCCACTTCTGACACCA'

class TesttRNA(TestCase):
    def setUp(self):
        self.trna = tRNA()
        self.trna.LoadStr(tRNA_str)

class TestGetTrfType(TesttRNA):
    # def test_get_tab_title(self):
    #     self.fail()
    #
    # def test_get_tab_str(self):
    #     self.fail()

    def test_load_str(self):
        self.assertEqual(self.trna.chrom, 'chr6')
        #self.fail()
    #
    # def test_get_mature_seq(self):
    #     self.fail()
    #
    # def test_get_pre_mature_no_intron_seq(self):
    #     self.fail()
    #
    # def test_get_intron_location_in_i(self):
    #     self.fail()
    #
    def test_gettrftype(self):
        keysites = self.trna.GetKeySitesInfor('C')
        type = getTRFType(keysites,53,87)
        print(type)
        self.assertEqual(type, '3_tRF')

    def test_get_key_sites_infor(self):
        keysites = self.trna.GetKeySitesInfor('C')
        print(keysites)
        self.assertEqual(keysites, [1,36,38,87])
    #
    # def test_is_qualifiedt_rna(self):
    #     self.fail()
    #
    # def test_calculate_alignment_loc_in_i(self):
    #     self.fail()
    #
    # def test_calculate_mut_loc_ito_m(self):
    #     self.fail()
    #
    # def test_find_pos_type(self):
    #     self.fail()
    #
    # def test_find_cleavage_sites(self):
    #     self.fail()
    #
    # def test_detect_peak(self):
    #     self.fail()
    #
    # def test_create_trf_profiles(self):
    #     self.fail()
    #
    # def test_brief_mapping_infor2mutation_dic(self):
    #     self.fail()
    #
    # def test_create_mutation_dic_str(self):
    #     self.fail()
