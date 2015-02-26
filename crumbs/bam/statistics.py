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

import os.path
from subprocess import Popen, PIPE
from operator import itemgetter
from itertools import izip
from array import array
from collections import Counter
import random

from crumbs.utils.optional_modules import (histogram, zeros, median,
                                           sum as np_sum)

from crumbs.utils.optional_modules import AlignmentFile
from crumbs.statistics import (draw_histogram_ascii, IntCounter, LABELS,
                               BestItemsKeeper)

from crumbs.settings import get_setting
from crumbs.bam.flag import SAM_FLAG_BINARIES, SAM_FLAGS
from crumbs.utils.bin_utils import get_binary_path
from crumbs.collectionz import RecentlyAddedCache
from crumbs.iterutils import generate_windows


# pylint: disable=C0111


DEFAULT_N_BINS = get_setting('DEFAULT_N_BINS')
DEFAULT_N_MOST_ABUNDANT_REFERENCES = get_setting('DEFAULT_N_MOST_ABUNDANT_REFERENCES')


def count_reads(ref_name, bams, start=None, end=None):
    'It returns the count of aligned reads in the region'
    count = 0
    for bam in bams:
        count += bam.count(reference=ref_name, start=start, end=end)
    return count


class ArrayWrapper(object):
    'A thin wrapper around numpy to have the same interface as IntCounter'
    def __init__(self, array, bins=DEFAULT_N_BINS):
        self.array = array
        self.labels = LABELS.copy()
        self._bins = bins

    @property
    def min(self):
        return self.array.min()

    @property
    def max(self):
        return self.array.max()

    @property
    def average(self):
        return self.array.mean()

    @property
    def median(self):
        return median(self.array)

    @property
    def variance(self):
        return self.array.var()

    @property
    def count(self):
        return len(self.array)

    @property
    def sum(self):
        return np_sum(self.array)

    def calculate_distribution(self, bins=None, min_=None, max_=None):
        if min_ is None:
            min_ = self.min
        if max_ is None:
            max_ = self.max

        if bins is None:
            bins = self._bins

        counts, bins = histogram(self.array, bins=bins, range=(min_, max_))
        return {'bin_limits': bins, 'counts': counts}

    def update_labels(self, labels):
        'It prepares the labels for output files'
        self.labels.update(labels)

    def __str__(self):
        return self.write()

    def write(self, max_in_distrib=None):
        'It writes some basic stats of the values'
        if self.count != 0:
            labels = self.labels
            # now we write some basic stats
            format_num = lambda x: '{:,d}'.format(x) if isinstance(x, int) else '%.2f' % x
            text = '{}: {}\n'.format(labels['minimum'], format_num(self.min))
            text += '{}: {}\n'.format(labels['maximum'], format_num(self.max))
            text += '{}: {}\n'.format(labels['average'],
                                      format_num(self.average))

            if labels['variance'] is not None:
                text += '{}: {}\n'.format(labels['variance'],
                                          format_num(self.variance))
            if labels['sum'] is not None:
                text += '{}: {}\n'.format(labels['sum'],
                                          format_num(self.sum))
            if labels['items'] is not None:
                text += '{}: {}\n'.format(labels['items'], self.count)
            text += '\n'
            distrib = self.calculate_distribution(max_=max_in_distrib,
                                                  bins=self._bins)
            text += draw_histogram_ascii(distrib['bin_limits'], distrib['counts'])
            return text
        return ''


class ReferenceStats(object):
    def __init__(self, bams,
                 n_most_abundant_refs=DEFAULT_N_MOST_ABUNDANT_REFERENCES,
                 bins=DEFAULT_N_BINS):
        self._bams = bams
        self._bins = bins
        self._rpkms = None
        self._tot_reads = 0
        self._lengths = None
        self._n_most_expressed_reads = n_most_abundant_refs
        self._most_abundant_refs = None
        self._count_reads()

    def _count_reads(self):
        nreferences = self._bams[0].nreferences
        rpks = zeros(nreferences)
        references = []
        length_counts = IntCounter()

        first_bam = True
        n_reads = 0
        for bam in self._bams:
            if bam.nreferences != nreferences:
                msg = 'BAM files should have the same references'
                raise ValueError(msg)
            for index, count in enumerate(get_reference_counts(bam.filename)):
                n_reads += count['unmapped_reads'] + count['mapped_reads']
                if count['reference'] is None:
                    # some non-mapped reads have reference = None
                    continue
                kb_len = count['length'] / 1000
                rpk = count['mapped_reads'] / kb_len
                rpks[index] += rpk
                if first_bam:
                    # For the reference lengths we use the first BAM to make
                    references.append(count['reference'])
                    length_counts[count['length']] += 1
                else:
                    # the bams should be sorted with the references in the same
                    # order
                    if references[index] != count['reference']:
                        msg = 'The reference lengths do not match in the bams'
                        raise RuntimeError(msg)
            first_bam = False

        million_reads = n_reads / 1e6
        rpks /= million_reads  # rpkms
        self._rpkms = ArrayWrapper(rpks, bins=self._bins)

        abundant_refs = BestItemsKeeper(self._n_most_expressed_reads,
                                        izip(references, rpks),
                                        key=itemgetter(1))
        abundant_refs = [{'reference': i[0], 'rpkm': i[1]} for i in abundant_refs]
        self._most_abundant_refs = abundant_refs

        self._lengths = length_counts

    @property
    def lengths(self):
        return self._lengths

    @property
    def rpkms(self):
        return self._rpkms

    @property
    def most_abundant_refs(self):
        return self._most_abundant_refs

    def __str__(self):
        return self.write()

    def write(self, max_rpkm=None):
        result = 'RPKMs\n'
        result += '-----\n'
        result += self.rpkms.write(max_in_distrib=max_rpkm)
        result += '\n'
        result += 'Most represented references\n'
        result += '---------------------------\n'
        result += ''.join(['{reference:s}: {rpkm:.5f}\n'.format(**r) for r in self.most_abundant_refs])
        result += '\n'
        result += 'Lengths\n'
        result += '-----\n'
        result += str(self.lengths)
        return result


def _flag_to_binary(flag):
    'It returns the indexes of the bits sets to 1 in the given flag'
    return [index for index, num in enumerate(SAM_FLAG_BINARIES) if num & flag]


class ReadStats(object):
    def __init__(self, bams):
        # TODO flag, read_group
        self._bams = bams
        self._mapqs = IntCounter()
        self._flag_counts = {}
        self._count_mapqs()

    def _count_mapqs(self):
        mapqs = self._mapqs
        flag_counts = [0] * len(SAM_FLAG_BINARIES)
        for bam in self._bams:
            for read in bam:
                if not read.is_unmapped:
                    mapqs[read.mapq] += 1
                for flag_index in _flag_to_binary(read.flag):
                    flag_counts[flag_index] += 1

        for count, flag_bin in zip(flag_counts, SAM_FLAG_BINARIES):
            self._flag_counts[SAM_FLAGS[flag_bin]] = count

    @property
    def mapqs(self):
        return self._mapqs

    @property
    def flag_counts(self):
        return self._flag_counts


class CoverageCounter(IntCounter):
    def __init__(self, bams):
        self._bams = bams
        self._count_cov()

    def _count_cov(self):
        for bam in self._bams:
            for column in bam.pileup():
                self[len(column.pileups)] += 1


def get_reference_counts_dict(bam_fpaths):
    'It gets a list of bams and returns a dict indexed by reference'
    counts = {}
    for bam_fpath in bam_fpaths:
        for line in get_reference_counts(bam_fpath):
            ref_name = line['reference']
            length = line['length']
            mapped_reads = line['mapped_reads']
            unmapped_reads = line['unmapped_reads']
            if ref_name not in counts:
                counts[ref_name] = {'mapped_reads': 0, 'unmapped_reads': 0,
                                    'length': length}
            assert length == counts[ref_name]['length']
            counts[ref_name]['mapped_reads'] += mapped_reads
            counts[ref_name]['unmapped_reads'] += unmapped_reads
    return counts


def get_reference_counts(bam_fpath):
    'Using samtools idxstats it generates dictionaries with read counts'
    cmd = [get_binary_path('samtools'), 'idxstats', bam_fpath]
    idx_process = Popen(cmd, stdout=PIPE)
    # we're not using pysam.idxstats here because the stdout differed
    # depending on how the tests were run
    for line in idx_process.stdout:
        ref_name, ref_length, mapped_reads, unmapped_reads = line.split()
        if ref_name == '*':
            ref_name = None
            ref_length = None
        else:
            ref_length = int(ref_length)
        yield {'reference': ref_name, 'length': ref_length,
               'mapped_reads': int(mapped_reads),
               'unmapped_reads': int(unmapped_reads)}


MAPQS_TO_CALCULATE = (0, 20, 30, 40)


def get_rgs_from_samfiles(bams):
    rgs = {}
    for bam in bams:
        readgroups = get_bam_readgroups(bam)
        if not readgroups:
            continue
        for rg in readgroups:
            rgs[rg['ID']] = rg
    if not rgs:
        rgs[None] = {'LB': None, 'ID': None, 'PL': None, 'SM': None}
    return rgs


def calculate_window(start, end, window_len, seq_len):
    'Given a region it creates a window around with the desired length'
    snv_len = end - start
    if snv_len >= window_len:
        start = start
        end = end
    else:
        win_out_snv_len = window_len - snv_len
        to_add_left = win_out_snv_len // 2
        to_add_right = to_add_left
        if win_out_snv_len % 2:
            if random.choice([True, False]):
                to_add_left += 1
            else:
                to_add_right += 1
        start -= to_add_left
        end += to_add_right
    if start < 0:
        start = 0
        end += abs(start)
    if end > seq_len:
        start -= end - seq_len
        end = seq_len
    return start, end


class BamCoverages(object):
    def __init__(self, bam_fpaths, min_mapq=None, window=1,
                 sampling_win_step=1,
                 bam_pileup_stepper='all', bam_rg_field_for_vcf_sample='SM'):
        self.min_mapq = min_mapq
        self.window = window
        self.bam_pileup_stepper = bam_pileup_stepper
        self.bam_rg_field_for_vcf_sample = bam_rg_field_for_vcf_sample
        self._cov_cache = RecentlyAddedCache(window * 2)
        self.sampling_win_step = sampling_win_step
        self._bams = []
        self._rgs = {}
        self._ref_lens = {}
        self._prepare_bams(bam_fpaths)

    def _prepare_bams(self, bam_fpaths):
        bams = []
        rgs = {}
        for idx, bam_fpath in enumerate(bam_fpaths):
            bam = AlignmentFile(bam_fpath)
            rgs_ = get_bam_readgroups(bam)
            if rgs_ is None:
                rgs_ = [{'ID': None,
                         self.bam_rg_field_for_vcf_sample: str(None)}]
            bams.append({'bam': bam, 'rgs': rgs_})
            for read_group in rgs_:
                read_group['bam'] = idx
                rgs[read_group['ID']] = read_group
        self._bams = bams
        self._rgs = rgs
        bam = bams[0]['bam']
        # We have to assume that all bmas have the same references
        ref_lens = {ref: le_ for ref, le_ in zip(bam.references, bam.lengths)}
        self._ref_lens = ref_lens

    def _count_reads_in_column(self, column, min_mapq, bam):
        reads = Counter()
        for pileup_read in column.pileups:
            alig_read = pileup_read.alignment
            if min_mapq is not None:
                read_mapq = alig_read.mapq
                if read_mapq < min_mapq:
                    continue

            n_rgs = len(bam['rgs'])
            if not n_rgs:
                sample = str(None)
            elif n_rgs == 1:
                sample_field = self.bam_rg_field_for_vcf_sample
                sample = bam['rgs'][0][sample_field]
            else:
                rg_id = [tag[1] for tag in alig_read.tags if tag[0] == 'RG']
                if not rg_id:
                    sample = None
                else:
                    rg_id = rg_id[0]
                    sample_field = self.bam_rg_field_for_vcf_sample
                    sample = self._rgs[rg_id][sample_field]
            reads[sample] += 1
        return reads

    def _calculate_coverages_in_pos(self, chrom, pos):
        cache = self._cov_cache
        if (chrom, pos) in cache:
            return cache[(chrom, pos)]

        min_mapq = self.min_mapq
        start = pos
        end = pos + 1
        reads = Counter()
        for bam in self._bams:
            columns = bam['bam'].pileup(reference=chrom, start=start,
                                        end=end,
                                        stepper=self.bam_pileup_stepper,
                                        truncate=True)
            for column in columns:
                reads.update(self._count_reads_in_column(column, min_mapq,
                                                         bam))

        cache[(chrom, pos)] = reads
        return reads

    def _calculate_coverages_in_win(self, chrom, start, end):
        counts = Counter()
        n_pos = 0
        for pos in xrange(start, end, self.sampling_win_step):
            cnts_in_pos = self._calculate_coverages_in_pos(chrom, pos)
            counts.update(cnts_in_pos)
            n_pos += 1
        if n_pos:
            cnts = {sample: cnt / n_pos for sample, cnt in counts.items()}
        else:
            cnts = {}
        return cnts

    def calculate_coverage_in_pos(self, chrom, pos):
        start, end = calculate_window(pos, pos + 1, self.window,
                                      self._ref_lens[chrom])
        return self._calculate_coverages_in_win(chrom, start, end)

    def _calculate_complete_coverage_distrib(self, region):
        if region is None:
            chrom, start, end = None, None, None
        else:
            chrom, start, end = region

        min_mapq = self.min_mapq
        covs = {sample: IntCounter() for sample in self.samples}
        covs[None] = IntCounter()
        for bam in self._bams:
            columns = bam['bam'].pileup(reference=chrom, start=start, end=end,
                                        stepper=self.bam_pileup_stepper,
                                        truncate=True)
            for column in columns:
                col_counts = self._count_reads_in_column(column, min_mapq,
                                                         bam)
                for sample, sample_cov in col_counts.items():
                    covs[sample][sample_cov] += 1
        return covs

    def calculate_coverage_distrib_in_region(self, region=None):
        if region is None:
            if self.window == 1:
                regions = None
            else:
                regions = [(ref, 0, le_ - 1) for ref, le_ in self._ref_lens.items()]
        else:
            regions = [region]

        if self.window == 1:
            if regions is None:
                region = None
            else:
                region = regions[0]
            return self._calculate_complete_coverage_distrib(region)

        counts = {}
        for region in regions:
            chrom, start, end = region
            for start, end in generate_windows(self.window, start=0,
                                               end=self._ref_lens[chrom],
                                               step=1):
                counts_in_win = self._calculate_coverages_in_win(chrom, start,
                                                                 end)
                for sample, cnts_in_win in counts_in_win.items():
                    if sample not in counts:
                        counts[sample] = IntCounter()
                    counts[sample][int(round(cnts_in_win))] += 1

        return counts

    @property
    def samples(self):
        samples = []
        for rg_id in self._rgs:
            sample_field = self.bam_rg_field_for_vcf_sample
            samples.append(self._rgs[rg_id][sample_field])
        return samples


def get_genome_coverage(bam_fhands):
    coverage_hist = IntCounter()
    for bam_fhand in bam_fhands:
        bam_fpath = bam_fhand.name
        cmd = [get_binary_path('bedtools'), 'genomecov', '-ibam', bam_fpath]
        cover_process = Popen(cmd, stdout=PIPE)
        for line in cover_process.stdout:
            if line.startswith('genome'):
                cov, value = line.split('\t')[1: 3]
                coverage_hist[int(cov)] += int(value)
    return coverage_hist


def counter_to_scatter_group(coverage_hist):
    # convert histohgram to the format that scatter_draw understands
    scatter_group = {'x': array('l'), 'y': array('l')}
    for integer in range(0, coverage_hist.max + 1):
        scatter_group['x'].append(integer)
        scatter_group['y'].append(coverage_hist[integer])

    return scatter_group


def get_bam_readgroups(bam):
    header = bam.header
    if 'RG' not in header:
        return None
    readgroups = []
    for rg in header['RG']:
        readgroups.append(rg)
    return readgroups


def get_rg_from_alignedread(read):
    rgid = [value for key, value in read.tags if key == 'RG']
    return None if not rgid else rgid[0]


def mapped_count_by_rg(bam_fpaths, mapqx=None):
    do_mapqx = True if mapqx is not None else False
    counter_by_rg = {}
    for bam_fpath in bam_fpaths:
        bam = AlignmentFile(bam_fpath, 'rb')
        readgroups = get_bam_readgroups(bam)
        if readgroups is None:
            bam_basename = os.path.splitext(os.path.basename(bam_fpath))[0]
            readgroups = [bam_basename]
        else:
            readgroups = [rg['ID'] for rg in readgroups]
        for readgroup in readgroups:
            counter = IntCounter({'unmapped': 0, 'mapped': 0})
            if do_mapqx:
                counter['bigger_mapqx'] = 0
            counter_by_rg[readgroup] = counter

        for read in bam:
            rg = get_rg_from_alignedread(read)
            if rg is None:
                rg = bam_basename
            if do_mapqx and read.mapq >= mapqx:
                counter_by_rg[rg]['bigger_mapqx'] += 1
            if read.is_unmapped:
                counter_by_rg[rg]['unmapped'] += 1
            else:
                counter_by_rg[rg]['mapped'] += 1
    return counter_by_rg
