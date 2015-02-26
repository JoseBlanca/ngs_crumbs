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

from __future__ import division

import json
from os.path import join, abspath
from collections import Counter
from itertools import count
from array import array

from crumbs.vcf.prot_change import (get_amino_change, IsIndelError,
                                    BetweenSegments, OutsideAlignment)
from crumbs.vcf.snv import VCFReader
from crumbs.iterutils import RandomAccessIterator
from crumbs.utils.optional_modules import (parse_into_seqrecs, seq_index,
                                           CommOnly, RestrictionBatch,
                                           Analysis, Figure)
from crumbs.settings import get_setting
from crumbs.vcf.filters import _print_figure
from crumbs.statistics import calculate_dust_score
from crumbs.utils.tags import SEQRECORD
from crumbs.seq.seq import SeqWrapper
from crumbs.bam.statistics import calculate_window

DATA_DIR = abspath(join(__file__, '..', 'data'))
# Missing docstring
# pylint: disable=C0111

COMMON_ENZYMES = ['EcoRI', 'SmaI', 'BamHI', 'AluI', 'BglII', 'SalI', 'BglI',
                  'ClaI', 'TaqI', 'PstI', 'PvuII', 'HindIII', 'EcoRV',
                  'HaeIII', 'KpnI', 'ScaI', 'HinfI', 'DraI', 'ApaI', 'BstEII',
                  'ZraI', 'BanI', 'Asp718I']

MIN_READS = 6
MIN_READS_PER_ALLELE = 3


class BaseAnnotator(object):
    '''Base class for the Filter and Info annotator.

    It can be considered as a mapper for a SNP that modifies its filter and
    info fields.
    If an annotator object returns True to is_filter it will write its property
    name when the __call__ method considers it. If is_filter returns True the
    annotator won't write in the filter field but only in the info field.
    It will be also the __call__ method the one that should write in the info
    field. An annotator could write both in the info and the filter fields at
    the same time.
    The properties name and description are used to build the header
    description for the filters and the info property is used to write the
    header for the info field.
    '''

    return_modified_snv = False

    def __call__(self, record):
        raise NotImplementedError()

    @property
    def encode_conf(self):
        encoded_conf = json.dumps(self.conf).replace("\"", "'")
        return encoded_conf

    def _clean_filter(self, record):
        name = self.name
        if record.filters is not None and name in record.filters:
            record.filters.remove(name)

    def _clean_info(self, record):
        if self.info_id in record.infos:
            del(record.infos[self.info.id])

    @property
    def is_filter(self):
        raise NotImplementedError()

    @property
    def item_wise(self):
        raise NotImplementedError()

    @property
    def info(self):
        return {}

    @property
    def info_id(self):
        if self.info:
            return self.info['id']
        return None

    @staticmethod
    def _create_reader_from_snv(snv):
        orig_reader = snv.reader
        fpath = orig_reader.fhand.name
        min_calls_for_pop_stats = orig_reader.min_calls_for_pop_stats
        random_reader = VCFReader(open(fpath),
                               min_calls_for_pop_stats=min_calls_for_pop_stats)
        return random_reader


class CloseToSnv(BaseAnnotator):
    '''Filter snps with other close snvs.

    Allowed snv_types:  [snp, indel, unknown] '''

    def __init__(self, distance=60, max_maf_depth=None, snv_type=None):
        self.distance = distance
        self.random_reader = None
        self.max_maf_depth = max_maf_depth
        self.snv_type = snv_type
        self.conf = {'distance': distance, 'max_maf_depth': max_maf_depth,
                     'snv_type': snv_type}

    def __call__(self, snv):
        if self.random_reader is None:
            self.random_reader = self._create_reader_from_snv(snv)

        self._clean_filter(snv)
        chrom = snv.chrom
        pos = snv.pos
        start = pos - self.distance if pos - self.distance > 0 else 0
        end = pos + self.distance
        snv_type = self.snv_type
        max_maf_depth = self.max_maf_depth
        passed_snvs = 0
        for snv_in_window in self.random_reader.fetch_snvs(chrom, start, end):
            if snv_in_window.pos == pos:
                continue
            if max_maf_depth is None and snv_type is None:
                passed_snvs += 1
            elif max_maf_depth is None and snv_type is not None:
                if snv_type == snv_in_window.var_type:
                    passed_snvs += 1
            elif max_maf_depth is not None and snv_type is None:
                calculated_maf = snv_in_window.maf_depth
                if calculated_maf and calculated_maf <= max_maf_depth:
                    passed_snvs += 1
            else:
                calculated_maf = snv_in_window.maf_depth
                if (calculated_maf and calculated_maf <= max_maf_depth and
                    snv_type == snv_in_window.kind):
                    passed_snvs += 1
        if passed_snvs:
            snv.add_filter(self.name)

        if self.return_modified_snv:
            return snv

    @property
    def name(self):
        snv_type = '' if self.snv_type is None else self.snv_type[0]
        maf = self.max_maf_depth
        maf_str = '' if maf is None else '_{:.2f}'.format(maf)
        return 'cs{}{}{}'.format(snv_type, self.distance, maf_str)

    @property
    def description(self):
        snv_type = 'snv' if self.snv_type is None else self.snv_type
        maf = self.max_maf_depth
        maf_str = ', with maf:{:.2f}'.format(maf) if maf else ''
        desc = 'The snv is closer than {} nucleotides to another {}{} :: {}'
        return desc.format(self.distance, snv_type, maf_str, self.encode_conf)

    @property
    def is_filter(self):
        return True

    @property
    def item_wise(self):
        return True


# TODO: use fai if it is available
def _get_lengths(fhand):
    lengths = {}
    for seq in parse_into_seqrecs(fhand, 'fasta'):
        lengths[seq.id] = len(seq)
    return lengths


class HighVariableRegion(BaseAnnotator):
    'Filter depending on the variability of the region'

    def __init__(self, max_variability, window_in_bp, ref_fpath):
        self.max_variability = max_variability

        self._lengths = _get_lengths(open(ref_fpath))

        if not window_in_bp % 2:
            raise ValueError('Window in bp must be a odd number')
        self.window_in_bp = window_in_bp

        max_num_snps = (int(window_in_bp * max_variability) + 1) * 2
        if not max_num_snps % 2:
            max_num_snps += 1
        self._max_num_snps = max_num_snps

        self.conf = {'max_variability': max_variability,
                     'window_in_bp': window_in_bp}

    def __call__(self, snvs):
        # TODO: Randon acess iterator based on physical distance
        # RandomAcessRegionIterator(items, location_getter)
        # once we do this we can remove max_num_snps
        snvs = RandomAccessIterator(snvs, self._max_num_snps)
        half_win = (snvs._rnd_access_win - 1) // 2
        half_win_in_bp = (self.window_in_bp - 1) // 2
        for idx, snv in enumerate(snvs):
            self._clean_filter(snv)
            chrom = snv.chrom
            pos = snv.pos
            start = idx - half_win
            if start < 0:
                start = 0
            end = idx + half_win
            snvs_in_win = snvs[start:end]

            def snv_is_close(snv2):
                if snv2.chrom != chrom:
                    return False
                if abs(snv2.pos - pos) < half_win_in_bp:
                    return True
                else:
                    return False

            close_snvs = filter(snv_is_close, snvs_in_win)

            num_snvs = len(close_snvs)

            win_len = self.window_in_bp
            # The studied window could be smaller than expected if it is
            # located at the beginning of the chromosome
            dist_from_0 = pos
            if dist_from_0 < half_win_in_bp:
                win_len -= (half_win_in_bp - dist_from_0)

            # The studied window could be smaller than expected if it is
            # located at the end of the chromosome
            ref_len = self._lengths[chrom]
            end = ref_len - 1
            len_not_studied_at_end = pos + half_win_in_bp - end
            if len_not_studied_at_end > 0:
                win_len -= len_not_studied_at_end

            freq = num_snvs / win_len
            if freq >= self.max_variability:
                snv.add_filter(self.name)
            yield snv

    @property
    def name(self):
        return 'hv{}'.format(self.max_variability)

    @property
    def description(self):
        desc = 'The region has more than {} snvs per {} bases:: {}'
        return desc.format(int(self.max_variability * 100),
                           self.window_in_bp, self.encode_conf)

    @property
    def is_filter(self):
        return True

    @property
    def item_wise(self):
        return False


class CloseToLimit(BaseAnnotator):
    'Filter snps close to chromosome limits'

    def __init__(self, distance, ref_fpath):
        self.distance = distance
        self._lengths = _get_lengths(open(ref_fpath))
        self.conf = {'distance': distance}

    def __call__(self, snv):
        self._clean_filter(snv)
        chrom = snv.chrom
        pos = snv.pos
        seq_len = self._lengths[chrom]
        if pos < self.distance or seq_len - pos < self.distance:
            snv.add_filter(self.name)
        if self.return_modified_snv:
            return snv

    @property
    def name(self):
        return 'cl{}'.format(self.distance)

    @property
    def description(self):
        desc = 'The snv is closer than {} nucleotides to the reference edge'
        desc += ':: {}'
        return desc.format(self.distance, self.encode_conf)

    @property
    def is_filter(self):
        return True

    @property
    def item_wise(self):
        return True


class MafDepthLimit(BaseAnnotator):
    'Filter by maf'

    def __init__(self, max_maf, samples=None):
        self.max_maf = max_maf
        self.samples = samples
        self.conf = {'max_maf': max_maf, 'samples': samples}

    def __call__(self, snv):
        self._clean_filter(snv)
        maf = snv.maf_depth
        if max and maf > self.max_maf:
            snv.add_filter(self.name)
        if self.return_modified_snv:
            return snv

    @property
    def name(self):
        return 'maf{}'.format(self.max_maf)

    @property
    def description(self):
        samples = 'all' if self.samples is None else ','.join(self.samples)
        desc = 'The most frequent allele in samples: {}.'
        desc += 'frequency greater than {} :: {}'
        return desc.format(samples, self.max_maf, self.encode_conf)

    @property
    def is_filter(self):
        return True

    @property
    def item_wise(self):
        return True


class CapEnzyme(BaseAnnotator):
    'Filters by detectablity by restriction enzymes'

    def __init__(self, all_enzymes, ref_fpath):
        self.all_enzymes = all_enzymes
        self.ref_index = seq_index(ref_fpath, 'fasta')
        self.conf = {'all_enzymes': all_enzymes}
        self._last_chrom = None

    def __call__(self, snv):
        self._clean_filter(snv)
        alleles = snv.alleles
        enzymes = set()
        # we have to make all the posible conbinations
        chrom = snv.chrom
        last_chrom = self._last_chrom
        if last_chrom is not None and chrom == self._last_chrom[0]:
            ref = last_chrom[1]
        else:
            ref = self.ref_index[snv.chrom]
            self._last_chrom = chrom, ref
        used_combinations = []
        for i_index in range(len(alleles)):
            for j_index in range(len(alleles)):
                allelei = alleles[i_index]
                allelej = alleles[j_index]
                if (i_index == j_index or
                    (allelei, allelej) in used_combinations or
                    (allelej, allelei) in used_combinations):
                    continue
                used_combinations.append((allelei, allelej))
                i_j_enzymes = _cap_enzymes_between_alleles(allelei, allelej,
                                                           ref, snv.pos,
                                                           snv.end,
                                                  all_enzymes=self.all_enzymes)
                enzymes = enzymes.union(i_j_enzymes)
        if not enzymes:
            snv.add_filter(self.name)
        else:
            enzymes = [str(e) for e in list(enzymes)]
            snv.add_info(info=self.info['id'], value=','.join(enzymes))

        if self.return_modified_snv:
            return snv

    @property
    def name(self):
        return 'ce{}'.format('t' if self.all_enzymes else 'f')

    @property
    def description(self):
        enzs = 'all' if self.all_enzymes else 'cheap_ones'
        desc = "SNV is not a CAP detectable by the enzymes: {}:: {}"
        return desc.format(enzs, self.encode_conf)

    @property
    def info(self):
        return {'id': 'CAP', 'num': 1, 'type': 'String',
                'desc': 'Enzymes that can be use to differentiate the alleles'}

    @property
    def is_filter(self):
        return True

    @property
    def item_wise(self):
        return True


def _cap_enzymes_between_alleles(allele1, allele2, reference, start, end,
                                 all_enzymes=False):
    '''It looks in the enzymes that differenciate the given alleles.

    It returns a set.
    '''

    start += 1
    # we have to build the two sequences
    if all_enzymes:
        restriction_batch = CommOnly
    else:
        restriction_batch = RestrictionBatch(COMMON_ENZYMES)

    sseq = reference.seq
    post_seq_start = start - 100 if start - 100 > 0 else 0
    prev_seq = sseq[post_seq_start: start - 1]
    post_seq = sseq[end: end + 100]

    seq1 = prev_seq + allele1 + post_seq
    seq2 = prev_seq + allele2 + post_seq
    anal1 = Analysis(restriction_batch, seq1, linear=True)
    enzymes1 = set(anal1.with_sites().keys())
    anal1 = Analysis(restriction_batch, seq2, linear=True)
    enzymes2 = set(anal1.with_sites().keys())
    enzymes = set(enzymes1).symmetric_difference(set(enzymes2))
    return enzymes


class Kind(BaseAnnotator):
    '''Choose Snv by type

    type options = snp, indel, unknown
    '''

    def __init__(self, kind='snp'):
        self.kind = kind
        self.conf = {'kind': kind}

    def __call__(self, snv):
        self._clean_filter(snv)
        rec_type = snv.kind
        if rec_type != self.kind:
            snv.add_filter(self.name)

        if self.return_modified_snv:
            return snv

    @property
    def name(self):
        return 'vk{}'.format(self.kind[0].lower())

    @property
    def description(self):
        return 'It is not an {} :: {}'.format(self.kind, self.encode_conf)

    @property
    def is_filter(self):
        return True

    @property
    def item_wise(self):
        return True


class IsVariableAnnotator(BaseAnnotator):
    _ids = count(0)

    def __init__(self, samples=None, consider_reference=False, filter_id=None,
                 min_samples_for_non_var=1):
        if filter_id is None:
            filter_id = self._ids.next()
        self.min_samples_for_non_var = min_samples_for_non_var
        self.filter_id = filter_id
        self.samples = samples
        self.consider_reference = consider_reference
        self.conf = {'samples': samples,
                     'consider_reference': consider_reference}

    def __call__(self, snv):
        samples = self.samples
        if samples is None:
            calls = snv.calls
        else:
            calls = [snv.get_call(sample) for sample in samples]

        alleles = set()
        if self.consider_reference:
            alleles.add(0)

        n_samples_called = 0
        for call in calls:
            int_alleles = call.int_alleles
            if int_alleles:
                n_samples_called += 1
            alleles.update(int_alleles)

        if not n_samples_called:
            is_variable = None
        elif len(alleles) == 1:
            if len(samples) >= self.min_samples_for_non_var:
                is_variable = False
            else:
                is_variable = None
        else:
            is_variable = True

        snv.add_info(info=self.info_id, value=str(is_variable))

        if self.return_modified_snv:
            return snv

    @property
    def info(self):
        if self.samples is None:
            samples = 'all'
        else:
            samples = ','.join(self.samples)
        desc = 'True if the genotypes are variable. False if they are not. '
        desc += 'None if there is not enough data in the samples: {samples}'
        desc += ' :: {conf}'
        desc = desc.format(samples=samples, conf=self.encode_conf)

        return {'id': 'IV{}'.format(self.filter_id), 'num': 1,
                'type': 'String', 'desc': desc}

    @property
    def is_filter(self):
        return False

    @property
    def item_wise(self):
        return True


class IsVariableDepthAnnotator(BaseAnnotator):
    'Variable in readgroup using depths and not genotypes'
    _ids = count(0)

    def __init__(self, samples, filter_id=None, max_maf_depth=None,
                 min_reads=MIN_READS, in_union=True,
                 min_reads_per_allele=MIN_READS_PER_ALLELE,
                 in_all_groups=True, reference_free=True):
        if filter_id is None:
            filter_id = self._ids.next()
        self.filter_id = filter_id
        self.max_maf = max_maf_depth
        self.samples = samples
        self.min_reads = min_reads
        self.min_reads_per_allele = min_reads_per_allele
        self.in_union = in_union
        self.in_all_groups = in_all_groups
        self.reference_free = reference_free
        self.conf = {'max_maf': max_maf_depth, 'samples': samples,
                     'min_reads': min_reads, 'in_union': in_union,
                     'min_reads_per_allele': min_reads_per_allele,
                     'in_all_groups': in_all_groups,
                     'reference_free': reference_free}

    def __call__(self, snv):
        is_variable = variable_in_samples_depth(snv, self.samples,
                                                self.in_union,
                                                self.max_maf,
                                                self.in_all_groups,
                                                self.reference_free,
                                                self.min_reads,
                                                self.min_reads_per_allele)

        snv.add_info(info=self.info_id, value=str(is_variable))

        if self.return_modified_snv:
            return snv

    @property
    def info(self):
        samples = ','.join(self.samples)
        desc = 'True if the samples are variable. False if they are not. '
        desc += 'None if there is not enougth data in the samples: {samples}'
        desc += ' :: {conf}'
        desc = desc.format(samples=samples, conf=self.encode_conf)

        return {'id': 'IVD{}'.format(self.filter_id), 'num': 1,
                'type': 'String', 'desc': desc}

    @property
    def is_filter(self):
        return False

    @property
    def item_wise(self):
        return True


def variable_in_samples_depth(snv, samples, in_union=True, maf_depth=None,
                              in_all_groups=True, reference_free=True,
                              min_reads=6, min_reads_per_allele=3):
    allele_counts = {s: snv.get_call(s).allele_depths for s in samples}

    if in_union:
        allele_counts_union = Counter()
        for sample_counts in allele_counts.values():
            for allele, counts in sample_counts.items():
                allele_counts_union[allele] += counts
        allele_counts = {'all':  allele_counts_union}

    variable_in_samples = []
    for sample_counts in allele_counts.values():
        variable_in_sample = _is_variable_in_sample_depth(sample_counts,
                                                          maf=maf_depth,
                                                          min_reads=min_reads,
                                                reference_free=reference_free,
                                     min_reads_per_allele=min_reads_per_allele)

        variable_in_samples.append(variable_in_sample)

    if None in variable_in_samples or not variable_in_samples:
        return None

    if in_all_groups:
        return all(variable_in_samples)
    else:
        return any(variable_in_samples)


def _is_variable_in_sample_depth(allele_counts, reference_free, maf,
                                 min_reads, min_reads_per_allele):
    'It checks if the allele is variable'
    num_reads = sum(allele_counts.values())
    if min_reads > num_reads:
        return None
    num_alleles = 0
    # first remove  the bad alleles
    for allele in allele_counts:
        if allele_counts[allele] >= min_reads_per_allele:
            num_alleles += 1

    if not num_alleles:
        return None

    ref_allele = 0

    if (ref_allele in allele_counts and
            allele_counts[ref_allele] >= min_reads_per_allele):
        ref_in_alleles = True
    else:
        ref_in_alleles = False

    if num_alleles == 1:
        if reference_free:
            return False
        elif not reference_free and ref_in_alleles:
            return False

    if allele_counts:
        counts = allele_counts.values()
        maf_allele = max(counts) / float(sum(counts))
    else:
        maf_allele = None
    # maf_allele = _calculate_maf_for_counts(allele_counts)
    if maf and maf < maf_allele:
        return False

    return True


class HeterozigoteInSamples(BaseAnnotator):

    def __init__(self, samples=None, filter_id=None, min_percent_het_gt=0,
                 gq_threshold=0, min_num_called=0):
        self._samples = samples
        self._min_percent_het_gt = min_percent_het_gt
        self._gq_threshold = gq_threshold
        self._min_num_called = min_num_called
        if filter_id is None:
            filter_id = self._ids.next()
        self.filter_id = filter_id
        self.conf = {'samples': samples,
                     'min_percent_het_get': min_percent_het_gt,
                     'gq_threslhold': gq_threshold,
                     'min_num_called': min_num_called}

    def __call__(self, snv):
        call_is_het = []
        for call in snv.calls:
            sample_name = call.sample
            if self._samples and sample_name not in self._samples:
                continue
            if call.gt_qual >= self._gq_threshold:
                call_is_het.append(call.is_het)
        num_calls = len(call_is_het)
        num_hets = len(filter(bool, call_is_het))
        if not num_calls or num_calls < self._min_num_called:
            result = None
        else:
            percent = int((num_hets / num_calls) * 100)
            if percent >= self._min_percent_het_gt:
                result = True
            else:
                result = False

        snv.add_info(info=self.info_id, value=str(result))

        if self.return_modified_snv:
            return snv

    @property
    def info(self):
        if self._samples is None:
            samples = 'all'
        else:
            samples = ','.join(self._samples)
        description = 'True if at least {min_perc}% of the called samples are '
        description += 'het. False if not. None if not enough data in the '
        description += 'samples {samples} :: {conf}'
        description = description.format(min_perc=self._min_percent_het_gt,
                                         samples=samples,
                                         conf=self.encode_conf)

        return {'id': 'HIS{}'.format(self.filter_id), 'num': 1,
                'type': 'String', 'desc': description}

    @property
    def is_filter(self):
        return False

    @property
    def item_wise(self):
        return True


class AminoChangeAnnotator(BaseAnnotator):
    'The aminoacid changes with the alternative allele'

    def __init__(self, ref_fpath, orf_seq_fpath):
        self.orf_suffix = '_orf_seq'
        self.ref_index = seq_index(ref_fpath, 'fasta')
        self.orf_seq_index = seq_index(orf_seq_fpath, 'fasta')

        self.conf = {}

    def __call__(self, snv):
        self._clean_filter(snv)
        seq_name = snv.chrom
        seq_ref = self.ref_index[seq_name]
        aminos = None
        try:
            seq_estscan = self.orf_seq_index[seq_name + self.orf_suffix]
        except KeyError:
            return
        try:
            aminos = get_amino_change(seq_ref, seq_estscan, snv)
        except IsIndelError:
            snv.add_filter(self.name)
            snv.add_info(info=self.info_id, value='INDEL')
        except BetweenSegments:
            pass
        except OutsideAlignment:
            pass
        if aminos is None:
            return

        if set(aminos['ref_amino']) != set(aminos['alt_amino']):
            snv.add_filter(self.name)
            info_val = '{}->{}'.format(aminos['ref_amino'],
                                       ','.join(aminos['alt_amino']))
            snv.add_info(info=self.info_id, value=info_val)

        if self.return_modified_snv:
            return snv

    @property
    def name(self):
        return 'caa'

    @property
    def description(self):
        return "Alt alleles change the amino acid"

    @property
    def info(self):
        return {'id': 'AAC', 'num': 1, 'type': 'String',
                'desc': 'Amino acid change'}

    @property
    def is_filter(self):
        return False

    @property
    def item_wise(self):
        return True


def _parse_blossum_matrix(path):
    matrix = {}
    fhand = open(path)
    header = fhand.readline().strip().split('\t')[1:]
    header = [h.strip() for h in header]
    for line in fhand:
        items = line.strip().split('\t')
        aa1 = items[0].strip()
        for aa2, value in zip(header, items[1:]):
            matrix[(aa1, aa2)] = int(value)
    return matrix


class AminoSeverityChangeAnnotator(BaseAnnotator):
    'Check if the change is severe  or not'

    def __init__(self, ref_fpath, orf_seq_fpath):
        self.orf_suffix = '_orf_seq'
        self.ref_index = seq_index(ref_fpath, 'fasta')
        self.orf_seq_index = seq_index(orf_seq_fpath, 'fasta')
        blossum_path = join(DATA_DIR, 'blossum90.csv')
        self.blosum = _parse_blossum_matrix(blossum_path)

        self.conf = {}

    def _is_severe(self, ref_aa, alt_aa):
        return True if self.blosum[(ref_aa, alt_aa)] < 0 else False

    def __call__(self, snv):
        self._clean_filter(snv)
        seq_name = snv.chrom
        seq_ref = self.ref_index[seq_name]
        aminos = None
        try:
            seq_estscan = self.orf_seq_index[seq_name + self.orf_suffix]
        except KeyError:
            return
        try:
            aminos = get_amino_change(seq_ref, seq_estscan, snv)
        except IsIndelError:
            snv.add_filter(self.name)
        except BetweenSegments:
            pass
        except OutsideAlignment:
            pass

        if aminos is None:
            return

        ref_aa = aminos['ref_amino']
        alt_aas = aminos['alt_amino']
        if set(ref_aa) == set(alt_aas):
            return snv
        if any([self._is_severe(ref_aa, alt_aa) for alt_aa in alt_aas]):
            snv.add_filter(self.name)

        if self.return_modified_snv:
            return snv

    @property
    def name(self):
        return 'cas'

    @property
    def description(self):
        return "Alt alleles change the amino acid and the change is severe"

    @property
    def is_filter(self):
        return False

    @property
    def item_wise(self):
        return True


DEFATULT_DUST_THRESHOLD = get_setting('DEFATULT_DUST_THRESHOLD')
DEF_SNP_DUST_WINDOW = get_setting('DEF_SNP_DUST_WINDOW')


class LowComplexityRegionAnnotator(BaseAnnotator):
    def __init__(self, ref_fpath, dust_threshold=DEFATULT_DUST_THRESHOLD,
                 window=DEF_SNP_DUST_WINDOW):
        self.threshold = dust_threshold
        self.window = window
        self.ref_index = seq_index(ref_fpath, 'fasta')
        self._last_chrom = None
        self._scores = array('f')

    def __call__(self, snv):
        self._clean_filter(snv)
        # we have to make all the posible conbinations
        chrom = snv.chrom
        last_chrom = self._last_chrom
        if last_chrom is not None and chrom == self._last_chrom[0]:
            ref = last_chrom[1]
        else:
            ref = self.ref_index[snv.chrom]
            self._last_chrom = chrom, ref

        start, end = calculate_window(snv.pos, snv.end, self.window, len(ref))

        snv_win_seq = SeqWrapper(SEQRECORD, ref[start:end], None)
        score = calculate_dust_score(snv_win_seq)
        if score > self.threshold:
            snv.add_filter(self.name)

        self._scores.append(score)

        if self.return_modified_snv:
            return snv

    @property
    def name(self):
        return 'lcr'

    def draw_hist(self, fhand):
        fig = Figure()
        axes = fig.add_subplot(111)
        axes.hist(self._scores, fill=True, log=True, bins=20,
                  rwidth=1)
        axes.axvline(x=self.threshold)
        axes.set_xlabel('% complexity score')
        axes.set_ylabel('num. SNPs')
        _print_figure(axes, fig, fhand, plot_legend=False)

    @property
    def is_filter(self):
        return True

    @property
    def description(self):
        return "SNV is located in a low complexity region"

    @property
    def item_wise(self):
        return True


def mark_snv_as_filter_passed(snv):
    if snv.filters is None:
        snv.filters = []
    return snv
