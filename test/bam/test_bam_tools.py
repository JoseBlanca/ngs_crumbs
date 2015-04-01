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

import os.path
import unittest
import shutil
from subprocess import check_output, check_call
from tempfile import NamedTemporaryFile

from pysam import AlignmentFile

from crumbs.utils.test_utils import TEST_DATA_DIR
from crumbs.utils.bin_utils import BAM_BIN_DIR
from crumbs.bam.bam_tools import (filter_bam, calmd_bam, realign_bam,
                                  index_bam, merge_sams,
                                  _downgrade_edge_qualities)

# pylint: disable=C0111


class SortTest(unittest.TestCase):
    def test_sort_bam_bin(self):
        bin_ = os.path.join(BAM_BIN_DIR, 'sort_bam')
        assert 'usage' in check_output([bin_, '-h'])

        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        sorted_fhand = NamedTemporaryFile(suffix='.sorted.bam')
        check_call([bin_, bam_fpath, '-o', sorted_fhand.name])
        assert "@HD\tVN:1.4" in check_output(['samtools', 'view', '-h',
                                              sorted_fhand.name])
        assert os.path.exists(sorted_fhand.name + '.bai')
        os.remove(sorted_fhand.name + '.bai')
        sorted_fhand.close()

        # no index
        sorted_fhand = NamedTemporaryFile()
        check_call([bin_, bam_fpath, '-o', sorted_fhand.name, '--no-index'])
        assert not os.path.exists(sorted_fhand.name + '.bai')

        # sort the sam file
        fhand = NamedTemporaryFile()
        fhand.write(open(bam_fpath).read())
        fhand.flush()
        check_call([bin_, fhand.name])
        assert "@HD\tVN:1.4" in check_output(['samtools', 'view', '-h',
                                              bam_fpath])
        assert os.path.exists(fhand.name + '.bai')
        os.remove(fhand.name + '.bai')


class ToolsTest(unittest.TestCase):
    def test_index_bam(self):
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        index_bam(bam_fpath)
        index_bam(bam_fpath)

    def test_merge_sam(self):
        bam_fpath = os.path.join(TEST_DATA_DIR, 'sample.bam')
        fhand = NamedTemporaryFile(suffix='.bam')
        out_fpath = fhand.name
        fhand.close()
        try:
            merge_sams([bam_fpath, bam_fpath], out_fpath=out_fpath)
            samfile = AlignmentFile(out_fpath)
            assert len(list(samfile)) == 2
            assert os.stat(bam_fpath) != os.stat(out_fpath)
        finally:
            if os.path.exists(out_fpath):
                os.remove(out_fpath)


class FilterTest(unittest.TestCase):
    def test_filter_mapq(self):
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        out_fhand = NamedTemporaryFile()
        filter_bam(bam_fpath, out_fhand.name, min_mapq=100)
        assert len(open(out_fhand.name).read(20)) == 20


class RealignTest(unittest.TestCase):
    def test_realign_bamself(self):
        ref_fpath = os.path.join(TEST_DATA_DIR, 'CUUC00007_TC01.fasta')
        bam_fpath = os.path.join(TEST_DATA_DIR, 'sample.bam')
        out_bam = NamedTemporaryFile()
        realign_bam(bam_fpath, ref_fpath, out_bam.name)

    def test_realign_bin(self):
        bin_ = os.path.join(BAM_BIN_DIR, 'realign_bam')
        assert 'usage' in check_output([bin_, '-h'])

        bam_fpath = os.path.join(TEST_DATA_DIR, 'sample.bam')
        ref_fpath = os.path.join(TEST_DATA_DIR, 'CUUC00007_TC01.fasta')
        realigned_fhand = NamedTemporaryFile(suffix='.realigned.bam')
        check_call([bin_, bam_fpath, '-o', realigned_fhand.name, '-f',
                    ref_fpath])
        assert open(realigned_fhand.name).read()

        # in parallel
        realigned_fhand = NamedTemporaryFile(suffix='.realigned.bam')
        check_call([bin_, bam_fpath, '-o', realigned_fhand.name, '-f',
                    ref_fpath])  # , '-t', '2'])
        assert open(realigned_fhand.name).read()


class CalmdTest(unittest.TestCase):
    def test_calmd_bam(self):
        ref_fpath = os.path.join(TEST_DATA_DIR, 'CUUC00007_TC01.fasta')
        bam_fpath = os.path.join(TEST_DATA_DIR, 'sample.bam')
        orig_qual = AlignmentFile(bam_fpath).next().qual
        try:
            out_bam = NamedTemporaryFile()
            calmd_bam(bam_fpath, ref_fpath, out_bam.name)

            samfile = AlignmentFile(out_bam.name)
            calmd_qual = samfile.next().qual
            assert orig_qual != calmd_qual
            assert calmd_qual == 'HHHHHHBHGGH!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        finally:
            if os.path.exists(out_bam.name):
                out_bam.close()

    def test_calmd_no_out(self):
        ref_fpath = os.path.join(TEST_DATA_DIR, 'CUUC00007_TC01.fasta')
        bam_fpath = os.path.join(TEST_DATA_DIR, 'sample.bam')
        copied_fpath = os.path.join(TEST_DATA_DIR, 'sample_copy.bam')
        try:
            shutil.copy(bam_fpath, copied_fpath)
            orig_stats = os.stat(copied_fpath)
            calmd_bam(copied_fpath, ref_fpath)
            calmd_stats = os.stat(copied_fpath)
            assert calmd_stats != orig_stats
        finally:
            if os.path.exists(copied_fpath):
                os.remove(copied_fpath)

    def test_calmd_bin(self):
        bin_ = os.path.join(BAM_BIN_DIR, 'calmd_bam')
        assert 'usage' in check_output([bin_, '-h'])

        bam_fpath = os.path.join(TEST_DATA_DIR, 'sample.bam')
        ref_fpath = os.path.join(TEST_DATA_DIR, 'CUUC00007_TC01.fasta')
        calmd_fhand = NamedTemporaryFile(suffix='.calmd.bam')
        check_call([bin_, bam_fpath, '-o', calmd_fhand.name, '-f',
                    ref_fpath])
        assert open(calmd_fhand.name).read()


class DowngradeQuality(unittest.TestCase):

    def test_downngrade_read_edges(self):
        # With softclip
        bam_fpath = os.path.join(TEST_DATA_DIR, 'sample.bam')
        sam = AlignmentFile(bam_fpath)

        aligned_read = sam.next()
        _downgrade_edge_qualities(aligned_read, size=4, qual_to_substract=30)
        res = [9, 9, 9, 9, 9, 9, 3, 9, 8, 8, 9, 9, 9, 9, 9, 39,
               39, 39, 38, 38, 36, 33, 36, 38, 36, 38, 38, 38, 38, 39, 39, 38,
               38, 38, 9, 9, 9, 9]
        assert list(aligned_read.query_qualities) == res

        # without softclip
        sam = AlignmentFile(os.path.join(TEST_DATA_DIR, 'seqs.bam'))

        aligned_read = sam.next()
        _downgrade_edge_qualities(aligned_read, size=4, qual_to_substract=30)
        expected = [11, 13, 11, 11, 37, 43, 43, 46, 46, 57, 57, 48, 57, 57, 42,
                    41, 32, 35, 38, 38, 38, 38, 41, 41, 39, 37, 37, 44, 42, 48,
                    47, 57, 47, 47, 48, 47, 57, 57, 54, 48, 57, 48, 54, 50, 50,
                    50, 50, 50, 57, 59, 54, 54, 54, 57, 57, 59, 57, 52, 52, 52,
                    52, 57, 57, 57, 57, 52, 52, 52, 52, 29, 27, 27, 22]

        assert list(aligned_read.query_qualities) == expected

        # reverse
        # rev seqs (sam specification puts all the alignment query
        # forward(cigar, seq, qual, ...). Reverse is inly noted in the flag
        bam_fpath = os.path.join(TEST_DATA_DIR, 'sample_rev.bam')
        sam = AlignmentFile(bam_fpath)

        aligned_read = sam.next()
        aligned_read = sam.next()
        aligned_read = sam.next()
        _downgrade_edge_qualities(aligned_read, size=4,
                                  qual_to_substract=30)
        res = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0]
        assert list(aligned_read.query_qualities[:14]) == res

    def test_downgrade_read_edges_binary(self):
        binary = os.path.join(BAM_BIN_DIR, 'downgrade_bam_edge_qual')
        bam_fpath = os.path.join(TEST_DATA_DIR, 'sample_rev.bam')
        with NamedTemporaryFile() as out_fhand:
            cmd = [binary, '-o', out_fhand.name, bam_fpath]
            check_call(cmd)
            sam = AlignmentFile(out_fhand.name)
            res = [0, 0]
            read = sam.next()
            assert list(read.query_qualities[:2]) == res
            assert read.get_tag('dl') == '8)5B'
            assert read.get_tag('dr') == '8?>>'

if __name__ == "__main__":
    # import sys; sys.argv = ['', 'DowngradeQuality']
    unittest.main()
