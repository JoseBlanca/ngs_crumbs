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
from numpy import histogram, zeros, median, sum as np_sum

import pysam
try:
    from pysam.csamtools import Samfile
except ImportError:
    from pysam import Samfile

from crumbs.statistics import (draw_histogram_ascii, IntCounter, LABELS,
                               BestItemsKeeper)

from crumbs.settings import get_setting
from crumbs.bam.flag import SAM_FLAG_BINARIES, SAM_FLAGS
from crumbs.utils.bin_utils import get_binary_path
from collections import Counter

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


class GenomeCoverages(object):
    def __init__(self, bam_fhands, mapqs=MAPQS_TO_CALCULATE):
        self._bams = [Samfile(bam_fhand.name) for bam_fhand in bam_fhands]
        self.mapqs_to_calculate = mapqs
        self.rgs = get_rgs_from_samfiles(self._bams)
        self._counters = {}
        for sample in self.samples:
            self._counters[sample] = {mapq: IntCounter() for mapq in mapqs}
        self._calculate()

    def __len__(self):
        return len([c for sample in self.samples for c in self._counters[sample].values()])

    def _calculate(self):
        for samfile in self._bams:
            for column in samfile.pileup(stepper='all', max_depth=100000):
                self._add(column)

    def _add(self, column):
        rgs = self.rgs
        column_coverages = {rg['SM']: Counter()for rg in rgs.values()}
        for pileup_read in column.pileups:
            alignment = pileup_read.alignment
            read_mapq = alignment.mapq
            rg = [tag[1] for tag in alignment.tags if tag[0] == 'RG']
            sample = rgs[rg[0]]['SM'] if rg else None
            for mapq_to_calc in self.mapqs_to_calculate:
                if read_mapq > mapq_to_calc:
                    column_coverages[sample][mapq_to_calc] += 1
        for sample, counts in column_coverages.items():
            for map_to_calc, count in counts.items():
                self._counters[sample][map_to_calc][count] += 1

    def get_per_sample_mapq_counters(self, sample):
        return self._counters.get(sample, None)

    @property
    def samples(self):
        return list(set([rg['SM'] for rg in self.rgs.values()]))


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
        bam = pysam.Samfile(bam_fpath, 'rb')
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
