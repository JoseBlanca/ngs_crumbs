# This file is part of seq_crumbs.
# seq_crumbs is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# seq_crumbs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with seq_crumbs. If not, see <http://www.gnu.org/licenses/>.

import unittest
import os
from cStringIO import StringIO
from subprocess import check_output, call
from tempfile import NamedTemporaryFile

from Bio.bgzf import BgzfReader
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from crumbs.seq.pairs import (match_pairs, interleave_pairs,
                              deinterleave_pairs,
                              group_pairs, group_pairs_by_name,
                              _parse_pair_direction_and_name_from_title,
                              _parse_pair_direction_and_name)
from crumbs.iterutils import flat_zip_longest
from crumbs.utils.tags import FWD, SEQRECORD, SEQITEM
from crumbs.utils.bin_utils import SEQ_BIN_DIR
from crumbs.utils.test_utils import TEST_DATA_DIR
from crumbs.seq.seq import get_str_seq
from crumbs.seq.seqio import read_seqs, assing_kind_to_seqs
from crumbs.exceptions import (InterleaveError, PairDirectionError,
                               ItemsNotSortedError)
from crumbs.seq.seq import SeqWrapper, SeqItem
from crumbs.seq.utils.file_formats import set_format

# pylint: disable=R0201
# pylint: disable=R0904
# pylint: disable=C0111


class PairMatcherTest(unittest.TestCase):
    'It tests the mate pair checker'

    def test_pair_matcher(self):
        'It test the pair matcher function'
        # with equal seqs but the last ones
        file1 = os.path.join(TEST_DATA_DIR, 'pairend1.sfastq')
        file2 = os.path.join(TEST_DATA_DIR, 'pairend2.sfastq')
        fwd_seqs = read_seqs([open(file1)])
        rev_seqs = read_seqs([open(file2)])

        out_fhand = StringIO()
        orphan_out_fhand = StringIO()
        out_format = 'fastq'
        seqs = flat_zip_longest(fwd_seqs, rev_seqs)
        match_pairs(seqs, out_fhand, orphan_out_fhand, out_format)

        output = out_fhand.getvalue()
        assert '@seq1:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG' in output
        assert '@seq2:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in output
        orp = orphan_out_fhand.getvalue()
        assert '@seq8:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in orp

        # with the firsts seqs different
        file1 = os.path.join(TEST_DATA_DIR, 'pairend1.sfastq')
        file2 = os.path.join(TEST_DATA_DIR, 'pairend3.sfastq')
        fwd_seqs = read_seqs([open(file1)])
        rev_seqs = read_seqs([open(file2)])
        out_fhand = StringIO()
        orphan_out_fhand = StringIO()
        out_format = 'fastq'
        seqs = flat_zip_longest(fwd_seqs, rev_seqs)
        match_pairs(seqs, out_fhand, orphan_out_fhand, out_format)

        output = out_fhand.getvalue()
        assert '@seq4:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG' in output
        assert '@seq5:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in output
        orp = orphan_out_fhand.getvalue()
        assert '@seq1:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG' in orp
        assert '@seq3:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in orp
        assert '@seq6:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in orp

        file1 = os.path.join(TEST_DATA_DIR, 'pairend4.sfastq')
        file2 = os.path.join(TEST_DATA_DIR, 'pairend2.sfastq')
        fwd_seqs = read_seqs([open(file1)])
        rev_seqs = read_seqs([open(file2)])
        out_fhand = StringIO()
        orphan_out_fhand = StringIO()
        out_format = 'fastq'

        seqs = flat_zip_longest(fwd_seqs, rev_seqs)
        match_pairs(seqs, out_fhand, orphan_out_fhand, out_format)

        output = out_fhand.getvalue()
        assert '@seq8:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG' in output
        assert '@seq8:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in output
        orp = orphan_out_fhand.getvalue()
        assert '@seq1:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in orp
        assert '@seq2:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in orp

        # with reads with no direcction
        file1 = os.path.join(TEST_DATA_DIR, 'pairend7.sfastq')
        file2 = os.path.join(TEST_DATA_DIR, 'pairend2.sfastq')
        fwd_seqs = read_seqs([open(file1)])
        rev_seqs = read_seqs([open(file2)])
        out_fhand = StringIO()
        orphan_out_fhand = StringIO()
        out_format = 'fastq'

        seqs = flat_zip_longest(fwd_seqs, rev_seqs)
        match_pairs(seqs, out_fhand, orphan_out_fhand, out_format)
        output = out_fhand.getvalue()
        assert '@seq8:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG' in output
        assert '@seq8:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in output
        assert '@seq1:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in output
        assert '@seq1:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in output

        orp = orphan_out_fhand.getvalue()
        assert '@seq6:136:FC706VJ:2:2104:15343:197393.mpl_1' in orp
        assert '@seq7:136:FC706VJ:2:2104:15343:197393.hhhh' in orp
        assert '@seq2:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCAC' in orp

        # File is not sorted
        file1 = '''@s1.f
AACCAGTCAAC
+
CCCFFFFFGHH
@s2.f
AACCAGTCAAC
+
CCCFFFFFGHH
@s1.r
AACCAGTCAAC
+
CCCFFFFFGHH
'''
        file1 = StringIO(file1)
        set_format(file1, 'fastq')
        seqs = read_seqs([file1])
        out_fhand = StringIO()
        orphan_out_fhand = StringIO()
        out_format = 'fastq'

        try:
            match_pairs(seqs, out_fhand, orphan_out_fhand, out_format,
                        check_order_buffer_size=10)
            output = out_fhand.getvalue()
            self.fail('ItemsNotSortedError error expected')
        except ItemsNotSortedError:
            pass

    @staticmethod
    def test_all_orphan():
        'All reads end up in orphan'
        seqs = [SeqRecord(Seq('ACT'), id='seq1'),
                SeqRecord(Seq('ACT'), id='seq2')]
        seqs = list(assing_kind_to_seqs(SEQRECORD, seqs, None))
        out_fhand = StringIO()
        orphan_out_fhand = StringIO()
        match_pairs(seqs, out_fhand, orphan_out_fhand, out_format='fasta')
        assert orphan_out_fhand.getvalue() == '>seq1\nACT\n>seq2\nACT\n'

        out_fhand = StringIO()
        orphan_out_fhand = StringIO()
        match_pairs(seqs, out_fhand, orphan_out_fhand, ordered=False,
                    out_format='fasta')
        assert '>seq1\nACT\n' in orphan_out_fhand.getvalue()
        assert '>seq2\nACT\n' in orphan_out_fhand.getvalue()

    @staticmethod
    def test_mate_pair_unorderer_checker():
        'It test the mate pair function'
        # with equal seqs but the last ones
        file1 = os.path.join(TEST_DATA_DIR, 'pairend1.sfastq')
        file2 = os.path.join(TEST_DATA_DIR, 'pairend2.sfastq')
        fhand = NamedTemporaryFile()
        fhand.write(open(file1).read())
        fhand.write(open(file2).read())
        fhand.flush()
        seqs = read_seqs([fhand])

        out_fhand = StringIO()
        orphan_out_fhand = StringIO()
        out_format = 'fastq'
        match_pairs(seqs, out_fhand, orphan_out_fhand, out_format,
                    ordered=False)

        output = out_fhand.getvalue()
        assert '@seq1:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG' in output
        assert '@seq2:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in output
        orp = orphan_out_fhand.getvalue()
        assert '@seq8:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in orp

        # with the firsts seqs different
        file1 = os.path.join(TEST_DATA_DIR, 'pairend1.sfastq')
        file2 = os.path.join(TEST_DATA_DIR, 'pairend3.sfastq')
        fhand = NamedTemporaryFile()
        fhand.write(open(file1).read())
        fhand.write(open(file2).read())
        fhand.flush()
        seqs = read_seqs([fhand])

        out_fhand = StringIO()
        orphan_out_fhand = StringIO()
        out_format = 'fastq'
        match_pairs(seqs, out_fhand, orphan_out_fhand, out_format,
                    ordered=False)

        output = out_fhand.getvalue()
        assert '@seq4:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG' in output
        assert '@seq5:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in output
        orp = orphan_out_fhand.getvalue()
        assert '@seq1:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG' in orp
        assert '@seq3:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in orp
        assert '@seq6:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in orp

        file1 = os.path.join(TEST_DATA_DIR, 'pairend4.sfastq')
        file2 = os.path.join(TEST_DATA_DIR, 'pairend2.sfastq')
        fhand = NamedTemporaryFile()
        fhand.write(open(file1).read())
        fhand.write(open(file2).read())
        fhand.flush()
        seqs = read_seqs([fhand])

        out_fhand = StringIO()
        orphan_out_fhand = StringIO()
        out_format = 'fastq'

        match_pairs(seqs, out_fhand, orphan_out_fhand, out_format,
                    ordered=False)

        output = out_fhand.getvalue()
        assert '@seq8:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG' in output
        assert '@seq8:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in output
        orp = orphan_out_fhand.getvalue()
        assert '@seq1:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in orp
        assert '@seq2:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in orp

        # unordered file
        file1 = os.path.join(TEST_DATA_DIR, 'pairend1.sfastq')
        file2 = os.path.join(TEST_DATA_DIR, 'pairend2_unordered.sfastq')
        fhand = NamedTemporaryFile()
        fhand.write(open(file1).read())
        fhand.write(open(file2).read())
        fhand.flush()
        seqs = read_seqs([fhand])

        out_fhand = StringIO()
        orphan_out_fhand = StringIO()
        out_format = 'fastq'

        match_pairs(seqs, out_fhand, orphan_out_fhand, out_format,
                    ordered=False)
        output = out_fhand.getvalue()
        assert '@seq1:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG' in output
        assert '@seq2:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in output
        orp = orphan_out_fhand.getvalue()
        assert '@seq8:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in orp

        # with reads with no direcction
        file1 = os.path.join(TEST_DATA_DIR, 'pairend7.sfastq')
        file2 = os.path.join(TEST_DATA_DIR, 'pairend2.sfastq')
        fhand = NamedTemporaryFile()
        fhand.write(open(file1).read())
        fhand.write(open(file2).read())
        fhand.flush()
        seqs = read_seqs([fhand])

        out_fhand = StringIO()
        orphan_out_fhand = StringIO()
        out_format = 'fastq'

        match_pairs(seqs, out_fhand, orphan_out_fhand, out_format,
                    ordered=False)
        output = out_fhand.getvalue()
        assert '@seq8:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG' in output
        assert '@seq8:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in output
        assert '@seq1:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in output
        assert '@seq1:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in output

        orp = orphan_out_fhand.getvalue()
        assert '@seq6:136:FC706VJ:2:2104:15343:197393.mpl_1' in orp
        assert '@seq7:136:FC706VJ:2:2104:15343:197393.hhhh' in orp
        assert '@seq2:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCAC' in orp

    def test_pair_direction_and_name(self):
        'it test the pair_name parser'
        title = 'seq8:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG'
        name, dir_ = _parse_pair_direction_and_name_from_title(title)
        assert name == 'seq8:136:FC706VJ:2:2104:15343:197393'
        assert dir_ == FWD

        title = 'seq8:136:FC706VJ:2:2104:15343:197393/1'
        name, dir_ = _parse_pair_direction_and_name_from_title(title)
        assert name == 'seq8:136:FC706VJ:2:2104:15343:197393'
        assert dir_ == FWD

        title = 'seq8:136:FC706VJ:2:2104:15343:197393.f'
        name, dir_ = _parse_pair_direction_and_name_from_title(title)
        assert name == 'seq8:136:FC706VJ:2:2104:15343:197393'
        assert dir_ == FWD

        title = 'seq8:136:FC706VJ:2:2104:15343:197393.mp12'
        try:
            name, dir_ = _parse_pair_direction_and_name_from_title(title)
            self.fail()
        except PairDirectionError:
            pass

        title = r'seq8:136:FC706VJ:2:2104:15343:197393\1'
        name, dir_ = _parse_pair_direction_and_name_from_title(title)
        assert name == 'seq8:136:FC706VJ:2:2104:15343:197393'
        assert dir_ == FWD

        # With SeqRecord
        seq = SeqRecord(id=r'seq8:136:FC706VJ:2:2104:15343:197393\1',
                        seq=Seq('ACT'))
        name, dir_ = _parse_pair_direction_and_name(SeqWrapper(SEQRECORD, seq,
                                                               None))
        assert name == 'seq8:136:FC706VJ:2:2104:15343:197393'
        assert dir_ == FWD


class PairMatcherbinTest(unittest.TestCase):
    'It test the matepair binary'
    def test_pair_matcher_bin(self):
        'It test the pair matcher binary'
        pair_matcher_bin = os.path.join(SEQ_BIN_DIR, 'pair_matcher')
        assert 'usage' in check_output([pair_matcher_bin, '-h'])

        in_fpath = os.path.join(TEST_DATA_DIR, 'pairend5.sfastq')
        out_fhand = NamedTemporaryFile()
        orphan_fhand = NamedTemporaryFile()
        check_output([pair_matcher_bin, '-o', out_fhand.name,
                      '-p', orphan_fhand.name, in_fpath])

        result = open(out_fhand.name).read()
        assert '@seq1:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG' in result
        assert '@seq1:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in result

        orp = open(orphan_fhand.name).read()
        assert '@seq8:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in orp

        # compressed output
        in_fpath = os.path.join(TEST_DATA_DIR, 'pairend5.sfastq')
        out_fhand = NamedTemporaryFile()
        orphan_fhand = NamedTemporaryFile()
        check_output([pair_matcher_bin, '-o', out_fhand.name,
                      '-p', orphan_fhand.name, in_fpath, '-Z'])
        result = BgzfReader(out_fhand.name).read(2000)
        assert '@seq1:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG' in result
        assert '@seq1:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in result

        orp = BgzfReader(orphan_fhand.name).read(2000)
        assert '@seq8:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in orp

        # unordered file
        in_fpath = os.path.join(TEST_DATA_DIR, 'pairend6.sfastq')
        out_fhand = NamedTemporaryFile()
        orphan_fhand = NamedTemporaryFile()
        check_output([pair_matcher_bin, '-o', out_fhand.name,
                      '-p', orphan_fhand.name, in_fpath, '-u'])

        result = open(out_fhand.name).read()
        assert '@seq1:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG' in result
        assert '@seq1:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in result

        orp = open(orphan_fhand.name).read()
        assert '@seq8:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in orp


class InterleavePairsTest(unittest.TestCase):
    'It tests the interleaving and de-interleaving of pairs'
    def test_interleave(self):
        'It interleaves two iterators with paired reads'
        file1 = os.path.join(TEST_DATA_DIR, 'pairend1.sfastq')
        file2 = os.path.join(TEST_DATA_DIR, 'pairend2.sfastq')
        fwd_seqs = list(read_seqs([open(file1)], 'fastq'))
        rev_seqs = list(read_seqs([open(file2)], 'fastq'))

        try:
            list(interleave_pairs(fwd_seqs, rev_seqs))
            self.fail('InterleaveError expected')
        except InterleaveError:
            pass

        # we skip the tests
        seqs = list(interleave_pairs(fwd_seqs, rev_seqs, skip_checks=True))
        assert len(seqs) == 8

        file1 = os.path.join(TEST_DATA_DIR, 'pairend1.sfastq')
        file2 = os.path.join(TEST_DATA_DIR, 'pairend1b.sfastq')
        fwd_seqs = read_seqs([open(file1)], 'fastq')
        rev_seqs = read_seqs([open(file2)], 'fastq')

        seqs = list(interleave_pairs(fwd_seqs, rev_seqs))
        assert len(seqs) == 8

    def test_deinterleave(self):
        'It de-interleaves an iterator of alternating fwd and rev reads'

        fhand1 = os.path.join(TEST_DATA_DIR, 'pairend1.sfastq')
        fhand2 = os.path.join(TEST_DATA_DIR, 'pairend1b.sfastq')
        fwd_seqs = read_seqs([open(fhand1)], 'fastq')
        rev_seqs = read_seqs([open(fhand2)], 'fastq')

        seqs = interleave_pairs(fwd_seqs, rev_seqs)
        out_fhand1 = StringIO()
        out_fhand2 = StringIO()
        out_format = 'fastq'
        deinterleave_pairs(seqs, out_fhand1, out_fhand2, out_format)
        result1 = out_fhand1.getvalue()
        result2 = out_fhand2.getvalue()
        assert result1.strip() == open(fhand1).read().strip()
        assert result2.strip() == open(fhand2).read().strip()


class InterleaveBinTest(unittest.TestCase):
    'test of the interleave and deinterleave'

    def test_binaries(self):
        'It test the binaries'
        interleave_bin = os.path.join(SEQ_BIN_DIR, 'interleave_pairs')
        deinterleave_bin = os.path.join(SEQ_BIN_DIR, 'deinterleave_pairs')
        assert 'usage' in check_output([interleave_bin, '-h'])
        assert 'usage' in check_output([deinterleave_bin, '-h'])

        in_fpath1 = os.path.join(TEST_DATA_DIR, 'pairend1.sfastq')
        in_fpath2 = os.path.join(TEST_DATA_DIR, 'pairend1b.sfastq')
        out_fhand = NamedTemporaryFile()
        check_output([interleave_bin, '-o', out_fhand.name, in_fpath1,
                      in_fpath2])

        result = open(out_fhand.name).read()
        assert '@seq5:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG' in result
        assert '@seq5:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG' in result

        out_fhand1 = NamedTemporaryFile()
        out_fhand2 = NamedTemporaryFile()
        check_output([deinterleave_bin, '-o', out_fhand1.name, out_fhand2.name,
                      out_fhand.name])
        assert open(in_fpath1).read() == open(out_fhand1.name).read()
        assert open(in_fpath2).read() == open(out_fhand2.name).read()

        out_fhand1 = NamedTemporaryFile()
        out_fhand2 = NamedTemporaryFile()
        check_output([deinterleave_bin, '-o', out_fhand1.name, out_fhand2.name,
                      out_fhand.name, '-Z'])

        assert open(in_fpath1).read() == BgzfReader(out_fhand1.name).read(2000)
        assert open(in_fpath2).read() == BgzfReader(out_fhand2.name).read(2000)

        # skip checks
        in_fpath1 = os.path.join(TEST_DATA_DIR, 'pairend1.sfastq')
        in_fpath2 = os.path.join(TEST_DATA_DIR, 'pairend2.sfastq')
        out_fhand = NamedTemporaryFile()
        stderr = NamedTemporaryFile()
        ret_code = call([interleave_bin, '-o', out_fhand.name, in_fpath1,
                         in_fpath2], stderr=stderr)
        assert int(ret_code)
        assert 'read names from a pair do not matc' in open(stderr.name).read()
        check_output([interleave_bin, '-o', out_fhand.name, '-s', in_fpath1,
                      in_fpath2])
        result = open(out_fhand.name).read()
        assert 'seq4:136:FC706VJ:2:2104:15343:197393' in result
        assert 'seq3:136:FC706VJ:2:2104:15343:197393' in result

    def test_version(self):
        'It can return its version number'
        guess_bin = os.path.join(SEQ_BIN_DIR, 'interleave_pairs')
        stderr = NamedTemporaryFile()
        check_output([guess_bin, '--version'], stderr=stderr)
        assert 'from seq_crumbs version:' in open(stderr.name).read()

        guess_bin = os.path.join(SEQ_BIN_DIR, 'deinterleave_pairs')
        stderr = NamedTemporaryFile()
        check_output([guess_bin, '--version'], stderr=stderr)
        assert 'from seq_crumbs version:' in open(stderr.name).read()


def _build_some_paired_seqs():
    seq1 = SeqWrapper(SEQITEM, SeqItem('s1', ['>s1.f\n', 'A\n']), 'fasta')
    seq2 = SeqWrapper(SEQITEM, SeqItem('s1', ['>s1.r\n', 'C\n']), 'fasta')
    seq3 = SeqWrapper(SEQITEM, SeqItem('s2', ['>s2.f\n', 'T\n']), 'fasta')
    seq4 = SeqWrapper(SEQITEM, SeqItem('s2', ['>s2.r\n', 'G\n']), 'fasta')
    seqs = seq1, seq2, seq3, seq4
    return seqs


class PairGrouperTest(unittest.TestCase):

    def test_pair_grouper(self):
        seqs = _build_some_paired_seqs()
        paired_seqs = list(group_pairs(seqs))

        assert [get_str_seq(s) for s in paired_seqs[0]] == ['A', 'C']
        assert [get_str_seq(s) for s in paired_seqs[1]] == ['T', 'G']
        assert len(paired_seqs) == 2

        seqs = _build_some_paired_seqs()
        paired_seqs = list(group_pairs(seqs, n_seqs_in_pair=1,
                           check_name_matches=True))
        assert [get_str_seq(s) for pair in paired_seqs for s in pair] == ['A',
                                                                 'C', 'T', 'G']

    def test_name_check(self):
        seqs = _build_some_paired_seqs()
        try:
            list(group_pairs(seqs, n_seqs_in_pair=4))
            self.fail('InterleaveError expected')
        except InterleaveError:
            pass

        seqs = _build_some_paired_seqs()
        paired_seqs = list(group_pairs(seqs, n_seqs_in_pair=4,
                           check_name_matches=False))
        assert [get_str_seq(s) for s in paired_seqs[0]] == ['A', 'C', 'T', 'G']

    def test_n_seqs_check(self):
        seqs = _build_some_paired_seqs()
        seqs = seqs[:-1]
        try:
            list(group_pairs(seqs, n_seqs_in_pair=2))
            self.fail('InterleaveError expected')
        except InterleaveError:
            pass

        paired_seqs = list(group_pairs(seqs, n_seqs_in_pair=2,
                           check_all_same_n_seqs=False))
        assert [get_str_seq(s) for s in paired_seqs[0]] == ['A', 'C']
        assert [get_str_seq(s) for s in paired_seqs[1]] == ['T']

    @staticmethod
    def test_empty_iter():
        paired_seqs = list(group_pairs([]))
        assert not paired_seqs


class PairNameGrouperTest(unittest.TestCase):

    def test_pair_grouper(self):
        seqs = _build_some_paired_seqs()
        paired_seqs = list(group_pairs_by_name(seqs))
        assert [get_str_seq(s) for s in paired_seqs[0]] == ['A', 'C']
        assert [get_str_seq(s) for s in paired_seqs[1]] == ['T', 'G']
        assert len(paired_seqs) == 2

        seqs = seqs[0], seqs[2], seqs[1], seqs[3]
        paired_seqs = list(group_pairs_by_name(seqs))
        assert [get_str_seq(s) for s in paired_seqs[0]] == ['A']
        assert [get_str_seq(s) for s in paired_seqs[1]] == ['T']
        assert [get_str_seq(s) for s in paired_seqs[2]] == ['C']
        assert [get_str_seq(s) for s in paired_seqs[3]] == ['G']
        assert len(paired_seqs) == 4

        seqs = _build_some_paired_seqs()
        seqs = seqs[:-1]
        paired_seqs = list(group_pairs_by_name(seqs))
        assert [get_str_seq(s) for s in paired_seqs[0]] == ['A', 'C']
        assert [get_str_seq(s) for s in paired_seqs[1]] == ['T']

        seqs = _build_some_paired_seqs()
        seqs = seqs[:-1]
        try:
            paired_seqs = list(group_pairs_by_name(seqs,
                                                   all_pairs_same_n_seqs=True))
            self.fail('InterleaveError expected')
        except InterleaveError:
            pass

    def test_no_name(self):
        seqs = _build_some_paired_seqs()
        seq = SeqWrapper(SEQITEM, SeqItem('s', ['>s\n', 'N\n']), 'fasta')

        seqs = seqs[0], seqs[1], seqs[2], seq, seqs[3]
        paired_seqs = list(group_pairs_by_name(seqs))
        assert [get_str_seq(s) for s in paired_seqs[0]] == ['A', 'C']
        assert [get_str_seq(s) for s in paired_seqs[1]] == ['T']
        assert [get_str_seq(s) for s in paired_seqs[2]] == ['N']
        assert [get_str_seq(s) for s in paired_seqs[3]] == ['G']

        seqs = _build_some_paired_seqs()
        seqs = seqs[0], seq, seqs[1], seqs[2], seqs[3]
        paired_seqs = list(group_pairs_by_name(seqs))
        assert [get_str_seq(s) for s in paired_seqs[0]] == ['A']
        assert [get_str_seq(s) for s in paired_seqs[1]] == ['N']
        assert [get_str_seq(s) for s in paired_seqs[2]] == ['C']
        assert [get_str_seq(s) for s in paired_seqs[3]] == ['T', 'G']

        seqs = _build_some_paired_seqs()
        seqs = seq, seqs[0], seqs[1], seqs[2], seqs[3]
        paired_seqs = list(group_pairs_by_name(seqs))
        assert [get_str_seq(s) for s in paired_seqs[0]] == ['N']
        assert [get_str_seq(s) for s in paired_seqs[1]] == ['A', 'C']
        assert [get_str_seq(s) for s in paired_seqs[2]] == ['T', 'G']

    @staticmethod
    def test_empty_iter():
        paired_seqs = list(group_pairs_by_name([]))
        assert not paired_seqs

if __name__ == '__main__':
    #import sys;sys.argv = ['', 'PairGrouperTest']
    unittest.main()
