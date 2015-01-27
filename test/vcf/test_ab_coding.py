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
from StringIO import StringIO
from tempfile import NamedTemporaryFile

from crumbs.vcf.ab_coding import ABCoder, ENOUGH_SUPPORT, NOT_ENOUGH_SUPPORT


# Method could be a function
# pylint: disable=R0201
# Too many public methods
# pylint: disable=R0904
# Missing docstring
# pylint: disable=C0111


class ABCodingTest(unittest.TestCase):
    VCF_HEADER = '''##fileformat=VCFv4.1
##fileDate=20090805
##source=mysnpprogram
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Read Depth">
##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Read Depth">
'''
    vcf = '''#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6
20\t11\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/.\t1/.\t0/0\t0/0\t1/1\t1/1
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t./.\t./.\t1/1\t0/1\t0/1\t0/1
20\t15\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t0/0\t./.\t1/1\t1/1
20\t16\t.\tG\tA\t29\tPASS\tNS=3\tGT\t1/1\t0/0\t1/1\t1/1\t0/0\t0/0
20\t17\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t2/2\t2/2\t1/1\t1/1
20\t18\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t1/1
'''

    def test_ab_coding(self):
        fhand = StringIO(self.VCF_HEADER + self.vcf)

        coder = ABCoder(fhand, parents_a=['S1'], parents_b=['S2'],
                        parent_index_threshold=0.9)
        assert coder.offspring == ['S3', 'S4', 'S5', 'S6']
        try:
            list(coder.recode_genotypes())
            self.fail('RuntimeError expected')
        except RuntimeError:
            pass

        fhand = StringIO(self.VCF_HEADER + self.vcf)

        coder = ABCoder(fhand, parents_a=['S1'], parents_b=['S2'],
                        parent_index_threshold=0.9)
        result = coder.recode_genotypes(samples=coder.offspring)
        string = ''
        for snp, geno in result:
            string += str(snp.POS) + ' '
            string += ','.join(''.join(geno) for geno in geno.values())
            string += '\n'
        assert string == '''11 AA,AA,BB,BB
16 AA,AA,BB,BB
17 AA,AA,BB,BB
'''
        assert sum(coder.log.values()) == 6
        assert coder.log[NOT_ENOUGH_SUPPORT] == 2
        assert coder.log[ENOUGH_SUPPORT] == 3
        fhand = StringIO()
        coder.write_log(fhand)
        assert '6 SNPs ' in fhand.getvalue()

    def _ab_result_to_str(self, result):
        string = ''
        for snp, geno in result:
            string += str(snp.POS) + ' '
            genos = [gt_ if gt_ else ('.', '.') for gt_ in geno.values()]
            string += ','.join(''.join(geno) for geno in genos)
            string += '\n'
        return string

    def test_smooth(self):
        vcf = '''#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6
20\t11\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t1/1\t1/1\t0/0
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t0/0\t1/1\t1/1\t0/0
20\t15\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t0/0\t1/1\t1/1\t0/0
20\t16\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t0/0\t1/1\t1/1\t0/0
20\t17\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t0/0\t1/1\t1/1\t0/0
20\t18\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t0/0\t1/1\t0/0\t0/0
20\t19\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t0/0\t1/1\t1/1\t1/0
20\t20\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t0/0\t1/1\t1/1\t1/1
20\t21\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t0/0\t1/1\t1/1\t0/0
20\t22\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t0/0\t1/1\t1/1\t1/1
20\t23\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t0/0\t1/1\t1/1\t0/0
'''
        fhand = StringIO(self.VCF_HEADER + vcf)

        coder = ABCoder(fhand, parents_a=['S1'], parents_b=['S2'],
                        parent_index_threshold=0.9, smooth_threhsold=0.5,
                        window=7)
        result = coder.recode_genotypes(samples=coder.offspring)

        expected = '''11 AA,BB,BB,AA
14 AA,BB,BB,AA
15 AA,BB,BB,AA
16 AA,BB,BB,AA
17 AA,BB,BB,AA
18 AA,BB,BB,AA
19 AA,BB,BB,..
20 AA,BB,BB,..
21 AA,BB,BB,AA
22 AA,BB,BB,..
23 AA,BB,BB,AA
'''
        assert self._ab_result_to_str(result) == expected

        fhand = StringIO(self.VCF_HEADER + vcf)
        coder = ABCoder(fhand, parents_a=['S1'], parents_b=['S2'],
                        parent_index_threshold=0.9, smooth_threhsold=0.6,
                        recomb_threshold=2)
        result = coder.recode_genotypes(samples=coder.offspring)
        expected = '''11 AA,BB,BB,AB
14 AA,BB,BB,AB
15 AA,BB,BB,AB
16 AA,BB,BB,AB
17 AA,BB,BB,AB
18 AA,BB,BB,AB
19 AA,BB,BB,AB
20 AA,BB,BB,AB
21 AA,BB,BB,AB
22 AA,BB,BB,AB
23 AA,BB,BB,AB
'''
        assert self._ab_result_to_str(result) == expected

        fhand = NamedTemporaryFile(suffix='.png')
        coder.plot_smooth_hist(fhand)

if __name__ == "__main__":
#     import sys;sys.argv = ['', 'FilterTest.test_close_to_filter']
    unittest.main()
