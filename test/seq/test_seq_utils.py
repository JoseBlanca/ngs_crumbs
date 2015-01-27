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
# pylint: disable=R0201
# pylint: disable=R0904
# pylint: disable=C0111

import unittest
from tempfile import NamedTemporaryFile
import os.path
from subprocess import check_output

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from crumbs.utils.bin_utils import SEQ_BIN_DIR
from crumbs.seq.seq import assing_kind_to_seqs, get_str_seq
from crumbs.seq.utils.seq_utils import (uppercase_length, ChangeCase,
                                        get_uppercase_segments)
from crumbs.utils.tags import SWAPCASE, UPPERCASE, LOWERCASE, SEQRECORD


class UppercaseLengthTest(unittest.TestCase):
    'It tests the uppercase character count'
    def test_uppercase_length(self):
        'It counts the number of uppercase letters in a string'
        assert uppercase_length('aCTaGGt') == 4
        assert uppercase_length('acagt') == 0


def _make_fhand(content=''):
    'It makes temporary fhands'
    fhand = NamedTemporaryFile()
    fhand.write(content)
    fhand.flush()
    return fhand


class MaskedSegmentsTest(unittest.TestCase):
    'It tests the lower case segments location functions'

    @staticmethod
    def test_masked_locations():
        'It test the masked locations function'

        assert list(get_uppercase_segments('aaATTTTTTaa')) == [(2, 8)]

        assert list(get_uppercase_segments('aaATTTaTTaa')) == [(2, 5), (7, 8)]

        assert list(get_uppercase_segments('AAATaaa')) == [(0, 3)]

        assert list(get_uppercase_segments('aaaaAAAA')) == [(4, 7)]

        seq = 'AATTaaTTaaTTT'
        assert list(get_uppercase_segments(seq)) == [(0, 3), (6, 7), (10, 12)]

        assert list(get_uppercase_segments('AATT')) == [(0, 3)]
        assert not list(get_uppercase_segments('aatt'))


class ChangeCaseTest(unittest.TestCase):
    'It tests the case change'
    def test_case_change(self):
        'It changes the case of the sequences'
        seqs = [SeqRecord(Seq('aCCg'), letter_annotations={'dummy': 'dddd'})]
        seqs = assing_kind_to_seqs(SEQRECORD, seqs, None)
        change_case = ChangeCase(action=UPPERCASE)
        strs = [get_str_seq(s) for s in change_case(seqs)]
        assert strs == ['ACCG']

        seqs = [SeqRecord(Seq('aCCg'))]
        seqs = assing_kind_to_seqs(SEQRECORD, seqs, None)
        change_case = ChangeCase(action=LOWERCASE)
        strs = [get_str_seq(s) for s in change_case(seqs)]
        assert strs == ['accg']

        seqs = [SeqRecord(Seq('aCCg'))]
        seqs = assing_kind_to_seqs(SEQRECORD, seqs, None)
        change_case = ChangeCase(action=SWAPCASE)
        strs = [get_str_seq(s) for s in change_case(seqs)]
        assert strs == ['AccG']

    def test_bin(self):
        'It tests the trim seqs binary'
        change_bin = os.path.join(SEQ_BIN_DIR, 'change_case')
        assert 'usage' in check_output([change_bin, '-h'])

        fastq = '@seq1\naTCgt\n+\n?????\n@seq2\natcGT\n+\n?????\n'
        fastq_fhand = _make_fhand(fastq)

        result = check_output([change_bin, '-a', 'upper', fastq_fhand.name])
        assert '@seq1\nATCGT\n+' in result


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'ChangeCaseTest.test_bin']
    unittest.main()
