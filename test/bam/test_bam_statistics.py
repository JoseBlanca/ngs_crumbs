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
from subprocess import check_output
from tempfile import NamedTemporaryFile

import pysam

from crumbs.utils.test_utils import TEST_DATA_DIR
from crumbs.utils.bin_utils import BAM_BIN_DIR
from crumbs.bam.statistics import (count_reads, ReferenceStats, ReadStats,
                                   CoverageCounter, _flag_to_binary,
                                   get_reference_counts,
                                   get_reference_counts_dict,
                                   get_genome_coverage, get_bam_readgroups,
                                   mapped_count_by_rg, BamCoverages)


# pylint: disable=R0201
# pylint: disable=R0904
# pylint: disable=C0111


class StatsTest(unittest.TestCase):
    def test_count(self):
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        bam = pysam.Samfile(bam_fpath)
        assert count_reads('reference1', [bam]) == 9
        assert count_reads('reference2', [bam]) == 9
        assert count_reads('reference1', [bam], start=1, end=10) == 0
        assert count_reads('reference1', [bam], start=0, end=500) == 9

    def test_reference_stats(self):
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        bam = pysam.Samfile(bam_fpath)
        refstats = ReferenceStats([bam], n_most_abundant_refs=1)
        rpkms = refstats.rpkms
        assert rpkms.min - 291715.28 < 0.1
        assert rpkms.max - 600240.1 < 0.1
        assert rpkms.average - 445977.7 < 0.1
        assert rpkms.median - 445977.7 < 0.1
        assert rpkms.variance - 23796889620.7 < 0.1
        assert rpkms.count == 2
        assert rpkms.sum - 891955.38 < 0.1
        assert refstats.most_abundant_refs[0]['reference'] == 'reference1'
        assert list(rpkms.calculate_distribution()['counts'])[0] == 1
        assert 'minimum:' in str(rpkms)
        assert 'Most represented' in str(refstats)

        refstats = ReferenceStats([bam, bam])
        n_max_expressed = len(set([i['reference'] for i in refstats.most_abundant_refs]))
        assert n_max_expressed == 2
        max_rpkm = refstats.most_abundant_refs[0]['rpkm']
        assert max_rpkm - refstats.rpkms.max < 0.1
        assert refstats.rpkms.max - 600240.1 < 0.1

    def test_ref_stats_bin(self):
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')

        bin_ = os.path.join(BAM_BIN_DIR, 'calculate_ref_stats')
        # help
        assert 'usage' in check_output([bin_, '-h'])

        assert 'RPKMs' in check_output([bin_, bam_fpath])

    def test_read_stats(self):
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        bam = pysam.Samfile(bam_fpath)
        stats = ReadStats([bam])
        assert stats.mapqs.count == 18
        assert stats.mapqs.min == 28
        assert stats.mapqs.max == 149
        assert stats.flag_counts['is_unmapped'] == 0

    def test_coverage_distrib(self):
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        bam = pysam.Samfile(bam_fpath)
        cov = CoverageCounter([bam])
        assert cov.count == 147
        assert cov.min == 6
        assert cov.max == 9

    def test_flag_to_binary(self):
        assert not _flag_to_binary(0)
        assert _flag_to_binary(1) == [0]
        assert _flag_to_binary(2) == [1]
        assert _flag_to_binary(1 | 2) == [0, 1]

    def test_ref_counts(self):
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        counts = list(get_reference_counts(bam_fpath))
        assert counts[2] == {'unmapped_reads': 0, 'reference': None,
                             'length': None, 'mapped_reads': 0}
        assert counts[1] == {'unmapped_reads': 0, 'reference': 'reference2',
                             'length': 1714, 'mapped_reads': 9}
        counts = get_reference_counts_dict([bam_fpath])
        assert None in counts.keys()
        assert 'reference2' in counts.keys()
        assert 'reference2' in counts.keys()

    def test_get_readgroup(self):
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        readgroups = get_bam_readgroups(pysam.Samfile(bam_fpath))
        assert readgroups == [{'LB': 'group1', 'ID': 'group1+454',
                               'PL': '454', 'SM': 'group1+454'},
                              {'LB': 'group2', 'ID': 'group2+454',
                               'PL': '454', 'SM': 'group2+454'}]

    def test_mapped_counts(self):
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        map_counts = mapped_count_by_rg([bam_fpath])
        assert map_counts['group1+454']['mapped'] == 9

        bam_fpath = os.path.join(TEST_DATA_DIR, 'sample_no_rg.bam')
        map_counts = mapped_count_by_rg([bam_fpath])
        assert map_counts['sample_no_rg']['mapped'] == 1
        assert 'bigger_mapqx' not in map_counts['sample_no_rg']
        map_counts = mapped_count_by_rg([bam_fpath], mapqx=20)
        assert map_counts['sample_no_rg']['bigger_mapqx'] == 1

        map_counts = mapped_count_by_rg([bam_fpath], mapqx=50)
        assert map_counts['sample_no_rg']['bigger_mapqx'] == 0

    def test_bin_mapped_counts(self):
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')

        bin_ = os.path.join(BAM_BIN_DIR, 'count_mapped_by_rg')
        cmd = [bin_, bam_fpath]
        output = check_output(cmd)
        assert "group2+454" in output


class GenomeCoverageTest(unittest.TestCase):

    def test_genome_cover(self):
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        scatter_group = get_genome_coverage([open(bam_fpath)])
        assert scatter_group.items() == [(0, 2400), (9, 144), (6, 3)]


class BamCoverageTest(unittest.TestCase):
    def test_bam_coverage(self):
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        cov = BamCoverages([bam_fpath], sampling_win_step=1)
        exp = {'group1+454': 9}
        assert cov._calculate_coverages_in_pos('reference1', 200) == exp
        res = cov.calculate_coverage_distrib_in_region(region=('reference1',
                                                               None, None))
        assert res['group1+454'] == {9: 73}
        cov = BamCoverages([bam_fpath], window=2, sampling_win_step=1)
        res = cov.calculate_coverage_distrib_in_region(region=('reference1',
                                                               None, None))
        assert res == {'group1+454': {9: 72, 5: 2}}

        cov = BamCoverages([bam_fpath], window=21, sampling_win_step=10)
        res = cov.calculate_coverage_distrib_in_region(region=('reference1',
                                                               None, None))
        assert res == {'group1+454': {9: 53, 3: 20, 6: 20}}

    def test_bin_draw_cov_hist(self):
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        binary = os.path.join(BAM_BIN_DIR, 'draw_coverage_hist')
        out_fhand = NamedTemporaryFile(suffix='.png')
        out_fhand2 = NamedTemporaryFile(suffix='.txt')
        cmd = [binary, bam_fpath, '-p', out_fhand.name, '-o', out_fhand2.name]
        check_output(cmd)
        assert 'group2+454' in open(out_fhand2.name).read()
        # raw_input(out_fhand.name)

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'BamCoverageTest']
    unittest.main()
