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


# pylint: disable=C0111

from crumbs.seq.seq import get_str_seq
from crumbs.seq.pairs import group_pairs_by_name, group_pairs
from crumbs.seq.seqio import read_seqs, write_seqs
from crumbs.utils.tags import SEQITEM
from crumbs.iterutils import sorted_items, unique, unique_unordered


def _seqitem_pairs_equal(pair1, pair2):
    if len(pair1) != len(pair2):
        return False
    else:
        for read1, read2 in zip(pair1, pair2):
            if not get_str_seq(read1) == get_str_seq(read2):
                return False
        return True


def _read_pairs(in_fhands, paired_reads):
    seqs = read_seqs(in_fhands, prefered_seq_classes=[SEQITEM])
    if paired_reads:
        pairs = group_pairs_by_name(seqs)
    else:
        pairs = group_pairs(seqs, n_seqs_in_pair=1)
    return pairs


class _PairKeyGetter(object):
    def __init__(self, use_length=None):
        self._use_length = use_length

    def __call__(self, pair):
        key = []
        for read in pair:
            seq = get_str_seq(read)
            if self._use_length is not None:
                seq = seq[:self._use_length]
            key.append(seq)
        return tuple(key)


def filter_duplicates(in_fhands, out_fhand, paired_reads, use_length=None,
                      n_seqs_packet=None, tempdir=None):
    if not in_fhands:
        raise ValueError('At least one input fhand is required')
    pairs = _read_pairs(in_fhands, paired_reads)
    get_pair_key = _PairKeyGetter(use_length=use_length)
    if n_seqs_packet is None:
        unique_pairs = unique_unordered(pairs, key=get_pair_key)
    else:
        sorted_pairs = sorted_items(pairs, key=get_pair_key, tempdir=tempdir,
                                    max_items_in_memory=n_seqs_packet)
        unique_pairs = unique(sorted_pairs, key=get_pair_key)
    for pair in unique_pairs:
        write_seqs(pair, out_fhand)
