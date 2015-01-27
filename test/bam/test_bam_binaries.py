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
import os

from subprocess import check_call, check_output
from tempfile import NamedTemporaryFile

from crumbs.utils.bin_utils import BAM_BIN_DIR
from crumbs.utils.test_utils import TEST_DATA_DIR


class BinTests(unittest.TestCase):
    def test_draw_coverage(self):
        bin_ = os.path.join(BAM_BIN_DIR, 'draw_coverage_hist')
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        fhand = NamedTemporaryFile(suffix='.png')
        out = check_output([bin_, bam_fpath, '-o', fhand.name])
        assert '147' in out

    def test_mapq_hist(self):
        bin_ = os.path.join(BAM_BIN_DIR, 'draw_mapq_hist')
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        fhand = NamedTemporaryFile(suffix='.png')
        null = open(os.devnull, 'w')
        check_call([bin_, bam_fpath, '-o', fhand.name], stdout=null)
        # raw_input(fhand.name)

        fhand = NamedTemporaryFile(suffix='.png')
        check_call([bin_, bam_fpath, '-o', fhand.name, '-t'], stdout=null)
        res = open(fhand.name).read()
        assert "[147 , 154[ (3):" in res

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
