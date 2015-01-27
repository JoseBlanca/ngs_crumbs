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
from tempfile import NamedTemporaryFile
from os.path import exists
from os import remove
from StringIO import StringIO

from crumbs.utils.file_utils import (compress_with_bgzip, uncompress_gzip,
                                     fhand_is_seekable,
                                     wrap_in_buffered_reader,
                                     index_vcf_with_tabix)
#                                          _build_template_fhand)
# Method could be a function
# pylint: disable=R0201
# Too many public methods
# pylint: disable=R0904
# Missing docstring
# pylint: disable=C0111


VCF = '''##fileformat=VCFv4.1
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number_of_Samples_With_Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total_Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele_Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral_Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership,build_129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2_membership">
##FILTER=<ID=q10,Description="Quality_below_10">
##FILTER=<ID=s50,Description="Less_than_50%_of_samples_have_data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype_Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read_Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype_Quality">
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
20 14370 rs6054257 G A 29 PASS NS=3;DP=14;AF=0.5;DB;H2 GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
20 17330 . T A 3 q10 NS=3;DP=11;AF=0.017 GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3 0/0:41:3
'''


class SeekableFileTest(unittest.TestCase):

    def test_is_seekable(self):
        'It tests wether the fhands are seekable or not'

        # StringIO
        fhand = StringIO('hola')
        assert fhand_is_seekable(fhand)

        # standard file
        fhand = NamedTemporaryFile()
        fhand.seek(0)
        assert fhand_is_seekable(fhand)

        # a wrapped BufferedReader
        fhand2 = wrap_in_buffered_reader(fhand)
        assert fhand_is_seekable(fhand2)

        # pylint: disable=R0903
        # pylint: disable=C0111
        class NonSeekable(object):
            'Just for testing'
            pass

        assert not fhand_is_seekable(NonSeekable())

        class NonSeekable2(object):
            def seek(self):
                pass

            def seekable(self):
                return False

        assert not fhand_is_seekable(NonSeekable2())


class CompressTest(unittest.TestCase):
    def test_bgzip_compression(self):
        orig = 'hola\ncaracola\n'
        orig_fhand = NamedTemporaryFile()
        orig_fhand.write(orig)
        orig_fhand.flush()

        compressed_fhand = NamedTemporaryFile(suffix='.gz')
        compress_with_bgzip(orig_fhand, compressed_fhand)

        compressed_fhand.seek(0)
        compressed = compressed_fhand.read()
        orig_fhand.seek(0)
        assert orig_fhand.read() == orig

        uncompressed_fhand = NamedTemporaryFile()
        uncompress_gzip(compressed_fhand, uncompressed_fhand)
        compressed_fhand.seek(0)
        assert compressed_fhand.read() == compressed

        uncompressed_fhand.seek(0)
        assert uncompressed_fhand.read() == orig

    def test_vcf_index(self):
        vcf = VCF.replace(' ', '\t')
        vcf_fhand = NamedTemporaryFile(suffix='.vcf')
        vcf_fhand.write(vcf)
        vcf_fhand.flush()
        compressed_fhand = NamedTemporaryFile(suffix='.vcf.gz')
        compress_with_bgzip(vcf_fhand, compressed_fhand)

        index_vcf_with_tabix(compressed_fhand.name)

        assert exists(compressed_fhand.name + '.tbi')
        remove(compressed_fhand.name + '.tbi')

if __name__ == "__main__":
#     import sys;sys.argv = ['', 'FilterTest.test_close_to_filter']
    unittest.main()
