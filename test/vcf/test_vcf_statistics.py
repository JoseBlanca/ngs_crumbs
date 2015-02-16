# Copyright 2012 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of ngs_crumbs.
# ngs_crumbs is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# ngs_crumbs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with ngs_crumbs. If not, see <http://www.gnu.org/licenses/>.

import unittest
from os.path import join
from tempfile import NamedTemporaryFile
from subprocess import check_output
from StringIO import StringIO

try:
    import pysam
except ImportError:
    pass

from vcf import Reader

from crumbs.vcf.snv import VCFReader
from crumbs.vcf.statistics import (VcfStats, HOM_REF, VCFcomparisons,
                                   _AlleleCounts2D, HOM_ALT, HOM, HET,
                                   draw_read_pos_stats,
                                   calc_snv_read_pos_stats)

from crumbs.utils.bin_utils import VCF_BIN_DIR
from crumbs.utils.test_utils import TEST_DATA_DIR

VARSCAN_VCF_PATH = join(TEST_DATA_DIR, 'sample.vcf.gz')
REF_PATH = join(TEST_DATA_DIR, 'sample_ref.fasta')
GATK_VCF_PATH = join(TEST_DATA_DIR, 'gatk_sample.vcf.gz')
FREEBAYES_VCF_PATH = join(TEST_DATA_DIR, 'freebayes_sample.vcf.gz')
FREEBAYES_MULTI_VCF_PATH = join(TEST_DATA_DIR, 'freebayes_multisample.vcf.gz')
GENERIC_VCF = join(TEST_DATA_DIR, 'generic.vcf.gz')


class TestVcfStats(unittest.TestCase):
    def test_vcf_stats(self):
        vcf_stats = VcfStats(VARSCAN_VCF_PATH, min_calls_for_pop_stats=2)
        assert vcf_stats.gt_quals(HET)[21] == 2
        assert vcf_stats.gt_quals(HOM)[3] == 25
        assert vcf_stats.gt_quals(HET).count == 53
        assert (0.28 - vcf_stats.heterozigosity_for_sample('pepo')) < 0.01
        assert vcf_stats.het_by_snp[0] == 46
        fpath = join(TEST_DATA_DIR, 'freebayes6.vcf.gz')
        vcf_stats = VcfStats(fpath)
        covertures = vcf_stats.depths.keys()
        for i in [8, 2, 100, 7]:
            assert i in covertures

    def test_only_gt_vcf(self):
        vcf_stats = VcfStats(GENERIC_VCF, min_calls_for_pop_stats=2)
        sample = 'BH_T_122_C4EGEACXX_6_250311606_X4'
        res = vcf_stats.heterozigosity_for_sample(sample)
        self.assertAlmostEqual(res, 0.16666666)


class AlleleCount2DTest(unittest.TestCase):
    def test_allele_count2d(self):
        allelecount = _AlleleCounts2D()
        allelecount.add(2, 3, (0, 0), 25)
        allelecount.add(2, 3, (1, 0), 25)
        allelecount.add(2, 3, (1, 1), 25)
        allelecount.add(2, 3, (0, 0), 25)
        allelecount.add(2, 3, (0, 1), 50)
        allelecount.add(2, 3, (2, 2), 25)
        allelecount.add(2, 3, (1, 0), 75)
        allelecount.add(2, 4, (1, 0), 75)

        assert allelecount.get_gt_count(2, 3, HOM_ALT) == 2
        assert allelecount.get_gt_count(2, 3, HOM_REF) == 2
        assert allelecount.get_gt_count(2, 3, HET) == 3
        assert allelecount.get_avg_gt_qual(2, 3, HOM_ALT) == 25
        assert allelecount.get_avg_gt_qual(2, 3, HOM_REF) == 25
        assert allelecount.get_avg_gt_qual(2, 3, HET) == 50

        allelecount.get_gt_depths_for_coverage(5)


class VCFcomparisonsTest(unittest.TestCase):
    def test_calculate_statistics(self):
        # with freebayes
        reader = Reader(filename=FREEBAYES_VCF_PATH)
        vcf_to_compare = VCFcomparisons(FREEBAYES_VCF_PATH)
        stats = vcf_to_compare.calculate_statistics(reader)
        assert stats['common'] == 944
        assert stats['uncalled'] == 0
        assert stats['different'] == 0
        assert stats['common_snps_prc'] == 100

        # with varscan
        reader = Reader(filename=VARSCAN_VCF_PATH)
        vcf_to_compare = VCFcomparisons(VARSCAN_VCF_PATH, samples=['mu16'])
        stats = vcf_to_compare.calculate_statistics(reader, samples=['mu16'])
        assert stats['common'] == 107
        assert stats['uncalled'] == 69
        assert stats['different'] == 0
        assert stats['common_snps_prc'] == 100

    def xtest_compare_vcfs_samples(self):
        binary = join(VCF_BIN_DIR, 'compare_vcfs_samples')
        assert 'usage' in check_output([binary, '-h'])
        samples_fhand = NamedTemporaryFile()
        samples_fhand.write('mu16\n')
        samples_fhand.flush()

        cmd = [binary, VARSCAN_VCF_PATH, '-s', samples_fhand.name,
               '-r', samples_fhand.name, '-v', VARSCAN_VCF_PATH]
        stats = check_output(cmd)
        result = 'common_snps_prc : 100.0\ndifferent : 0\ncommon : 107\n'
        result += 'uncalled : 69\n'
        assert stats == result

    def test_binary(self):
        binary = join(VCF_BIN_DIR, 'calculat_vcf_stats')


class ReadPosCoord(unittest.TestCase):

    def test_snv_read_pos_distrib(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
reference1\t187\trs6054257\tAA\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.
reference1\t210\t.\tA\tAA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
reference1\t215\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4
reference1\t230\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2
reference2\t350\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3
reference2\t400\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t./.:35:4\t0/2:17:2\t1/1:40:3
'''

        snvs = VCFReader(StringIO(vcf)).parse_snvs()
        bam_fpath = join(TEST_DATA_DIR, 'seqs.bam')
        sam = pysam.AlignmentFile(bam_fpath)
        stats = calc_snv_read_pos_stats(sam, snvs)
        print stats
        assert 'group1+454' in stats['5_read_pos_counts'].keys()
        assert '5_read_pos_boxplot' in stats
        assert '3_read_pos_boxplot' in stats

        fhand = NamedTemporaryFile(suffix='.png')
        draw_read_pos_stats(stats, fhand)
        # raw_input(fhand.name)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'ReadPosCoord']
    unittest.main()
