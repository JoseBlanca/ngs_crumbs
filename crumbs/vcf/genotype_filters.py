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
from collections import Counter
from StringIO import StringIO

from crumbs.vcf.snv import VCFReader, VCFWriter
from crumbs.statistics import IntCounter
from crumbs.plot import HistogramPlotter

# Missing docstring
# pylint: disable=C0111
# Too few public methods
# pylint: disable=R0903

DEF_PROB_AA_THRESHOLD = 0.9999
HW = 'hw'
RIL_SELF = 'ril_self'


def run_genotype_filters(in_fhand, out_fhand, gt_filters, plots_dir=None,
                         reader_kwargs=None):
    if reader_kwargs is None:
        reader_kwargs = {}

    reader_kwargs['filename'] = 'pyvcf_bug_workaround'
    reader_kwargs['compressed'] = False
    reader = VCFReader(in_fhand, **reader_kwargs)

    templa_reader = VCFReader(StringIO(reader.header))
    writer = VCFWriter(out_fhand, template_reader=templa_reader)

    for snv in reader.parse_snvs():
        for mapper in gt_filters:
            snv = mapper(snv)
        try:
            writer.write_snv(snv)
        except IOError, error:
            # The pipe could be already closed
            if 'Broken pipe' in str(error):
                break
            else:
                raise


class LowQualityGenotypeFilter(object):

    def __init__(self, min_qual):
        self._min_qual = min_qual
        self._scores = IntCounter()

    def __call__(self, snv):
        for call in snv.calls:
            self._scores[call.gt_qual] += 1
        return snv.remove_gt_from_low_qual_calls(min_qual=self._min_qual)

    def draw_hist(self, fhand):
        plot = HistogramPlotter([self._scores], xlabel='Genotype quality',
                                ylabel='num genotypes',
                                titles=['Genotype quality distribution'])
        axe = plot.axes[0]
        axe.axvline(self._min_qual, color='r')
        plot.write_figure(fhand)


class HetGenotypeFilter(object):

    def __call__(self, snv):
        return snv.remove_gt_from_het_calls()


def prob_aa_given_n_a_reads_hw(num_a_reads, freq_a_in_pop):
    'It assumes HW'
    # TODO fix for Backcross
    proba = freq_a_in_pop
    res = proba
    res /= proba + 0.5 ** (num_a_reads - 1) * (1 - proba)
    return res

RIL_FREQ_AA_CACHE = {}


def _prob_aa_ril_self(n_generation):
    if n_generation in RIL_FREQ_AA_CACHE:
        return RIL_FREQ_AA_CACHE[n_generation]

    freq_aa = (2 ** n_generation - 1) / 2 ** (n_generation + 1)
    RIL_FREQ_AA_CACHE[n_generation] = freq_aa
    return freq_aa


def prob_aa_given_n_a_reads_ril_self(num_a_reads, n_generation):
    probaa = _prob_aa_ril_self(n_generation)
    res = probaa
    res /= probaa + (0.5 ** num_a_reads) * (1 - 2 * probaa)
    return res

PROB_AA_FUNCS = {RIL_SELF: prob_aa_given_n_a_reads_ril_self,
                 HW: prob_aa_given_n_a_reads_hw}


class LowEvidenceAlleleFilter(object):
    def __init__(self, prob_aa_threshold=DEF_PROB_AA_THRESHOLD,
                 genotypic_freqs_method=HW, genotypic_freqs_kwargs=None):
        self._min_prob = prob_aa_threshold
        self.genotypic_freqs_method = genotypic_freqs_method
        if genotypic_freqs_kwargs is None:
            genotypic_freqs_kwargs = {}
        self.genotypic_freqs_kwargs = genotypic_freqs_kwargs
        self.log = Counter()

    def __call__(self, snv):
        genotypic_freqs_method = self.genotypic_freqs_method

        if genotypic_freqs_method == HW:
            allele_freqs = snv.allele_freqs
            if not allele_freqs:
                num_samples = len(snv.record.samples)
                self.log['not_enough_individuals'] += num_samples
                self.log['tot'] += num_samples
                def set_all_gt_to_none(call):
                    return call.copy_setting_gt(gt=None,
                                                return_pyvcf_call=True)
                return snv.copy_mapping_calls(set_all_gt_to_none)

        calls = []
        kwargs = self.genotypic_freqs_kwargs
        log = self.log
        for call in snv.calls:
            if not call.called:
                filtered_call = call.call
                log['was_not_called'] += 1
            else:
                alleles = call.int_alleles
                if len(set(alleles)) > 1:
                    filtered_call = call.call
                    log['was_het'] += 1
                else:
                    allele = alleles[0]
                    allele_depths = call.allele_depths
                    if not allele_depths:
                        msg = 'Allele depths are required for the lowEvidence'
                        msg += 'Allele filter'
                        raise RuntimeError(msg)
                    depth = call.allele_depths[allele]
                    if genotypic_freqs_method in PROB_AA_FUNCS:
                        prob_aa_func = PROB_AA_FUNCS[genotypic_freqs_method]
                        if genotypic_freqs_method == HW:
                            kwargs['freq_a_in_pop'] = allele_freqs[allele]
                        kwargs['num_a_reads'] = depth
                        prob = prob_aa_func(**kwargs)
                    else:
                        msg = 'Method not implemented for genotypic freqs: '
                        msg += genotypic_freqs_method
                        raise NotImplementedError(msg)
                    if prob >= self._min_prob:
                        filtered_call = call.call
                        log['enough_evidence'] += 1
                    else:
                        geno = call.call.data.GT[:-1] + '.'
                        filtered_call = call.copy_setting_gt(gt=geno,
                                                        return_pyvcf_call=True)
                        log['not_enough_evidence'] += 1
                log['tot'] += 1
            calls.append(filtered_call)
        return snv.copy_mapping_calls(calls)
