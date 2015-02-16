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

from collections import namedtuple, Counter, OrderedDict
from itertools import imap
from array import array
from warnings import warn

import matplotlib
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1 import make_axes_locatable

from vcf import Reader as pyvcfReader

from crumbs.iterutils import RandomAccessIterator
from crumbs.vcf.ld import _count_biallelic_haplotypes
from crumbs.vcf.filters import _print_figure

# Missing docstring
# pylint: disable=C0111

DEF_AB_CODING_WIN = 31
MORE_THAN_2_CODINGS = 'More than 2 codings'
NOT_ENOUGH_SUPPORT = 'Not enough support'
ENOUGH_SUPPORT = 'Enough support'
NO_INFO = 'No information'
DEF_AB_CODER_THRESHOLD = 0.9

AlleleCoding = namedtuple('AlleleCoding', ['A', 'B'])


class GetCoding(object):
    suspicious_no_parent_alleles = 200

    def __init__(self, parents_a, parents_b):
        self.parents_a = parents_a
        self.parents_b = parents_b
        self._no_alleles = Counter()

    def _get_unique_allele(self, genotypes, samples):
        gts = []
        samples_with_geno = []
        for sample in samples:
            gt = genotypes.get(sample, None)
            if gt is not None:
                samples_with_geno.append(sample)
                gts.append(gt)
        alleles = set(allele for gt in gts if gt for allele in gt if allele != '.')

        for sample in samples_with_geno:
            self._no_alleles[sample] = 0

        if not alleles:
            for sample in samples:
                self._no_alleles[sample] += 1
            result = None
        elif len(alleles) > 1:
            result = None
        elif len(alleles) == 1:
            result = list(alleles)[0]

        for sample, count in self._no_alleles.items():
            if count > self.suspicious_no_parent_alleles:
                msg = 'Sample with more than %d contiguous snps with no '
                msg += 'genotypes: %s'
                msg %= (self.suspicious_no_parent_alleles, sample)
                self._no_alleles[sample] = 0
                warn(msg, RuntimeWarning)
        return result

    def __call__(self, snp):
        gts = {call.sample: call.gt_alleles for call in snp.samples if call.called}
        allele_a = self._get_unique_allele(gts, self.parents_a)
        allele_b = self._get_unique_allele(gts, self.parents_b)
        if allele_a is None and allele_b is None:
            coding = None
        elif (allele_a is not None and allele_b is not None):
            if allele_a != allele_b:
                coding = AlleleCoding(allele_a, allele_b)
            else:
                coding = None
        else:
            coding = None
        return coding


class ABCoder(object):
    def __init__(self, vcf_fhand, parents_a, parents_b,
                 parent_index_threshold=DEF_AB_CODER_THRESHOLD,
                 offspring=None, window=DEF_AB_CODING_WIN,
                 smooth_threhsold=None, recomb_threshold=None):
        # TODO a min number of genotypes to evaluate anything
        self._reader = pyvcfReader(vcf_fhand)
        self._offspring = offspring
        self.window = window
        self.parent_index_threshold = parent_index_threshold
        self.log = Counter()
        self.indexes = array('f')

        samples = self._reader.samples
        parents = set(parents_a + parents_b)
        missing_parents = parents.difference(samples)
        if missing_parents:
            msg = 'Some parents are not found in the vcf file: '
            msg += ','.join(missing_parents)
            raise ValueError(msg)
        self.parents_a = parents_a
        self.parents_b = parents_b

        self.smooth_threhsold = smooth_threhsold
        self.recomb_threshold = recomb_threshold
        self._recombs = array('B')
        self._smoothes = array('f')

    @property
    def offspring(self):
        if self._offspring is not None:
            return self._offspring
        off = set(self._reader.samples)
        offspring = off.difference(self.parents_a).difference(self.parents_b)
        offspring = list(sorted(offspring))
        self._offspring = offspring
        return offspring

    def _deduce_coding(self, snp_and_coding, snp1_calls, snp2_idxs):
        votes = Counter()
        offspring = self.offspring
        for snp2_idx in snp2_idxs:
            snp2, coding2 = snp_and_coding[snp2_idx]
            if coding2 is None:
                continue
            snp2_calls = [snp2.genotype(sample) for sample in offspring]
            haplos = _count_biallelic_haplotypes(snp1_calls, snp2_calls,
                                                 return_alleles=True)
            if haplos is None:
                continue
            else:
                haplo_cnt, alleles = haplos
            alleles_in_major_haplo = {alleles.b: alleles.a,
                                      alleles.B: alleles.A}
            if haplo_cnt is None:
                continue

            if (coding2.A not in alleles_in_major_haplo or
               coding2.B not in alleles_in_major_haplo):
                # The offspring alleles in snp2 do not match the alleles
                # in the parents
                continue
            # This is allele A in snp 1
            allele1A = alleles_in_major_haplo[coding2.A]
            # This is allele B in snp 1
            allele1B = alleles_in_major_haplo[coding2.B]
            voted_coding1 = AlleleCoding(allele1A, allele1B)

            recomb_rate = (haplo_cnt.aB + haplo_cnt.Ab) / sum(haplo_cnt)
            weight = 2 * (0.5 - recomb_rate) if recomb_rate < 0.5 else 0
            votes[voted_coding1] += weight
        if not votes or sum(votes.values()) == 0:
            deduced_coding1 = None
            self.log[NO_INFO] += 1
        elif len(votes) > 2:
            deduced_coding1 = None
            self.log[MORE_THAN_2_CODINGS] += 1
        else:
            deduced_coding1 = votes.most_common(1)[0][0]
            index = votes[deduced_coding1] / sum(votes.values())
            self.indexes.append(index)
            if index < self.parent_index_threshold:
                deduced_coding1 = None
                self.log[NOT_ENOUGH_SUPPORT] += 1
            else:
                self.log[ENOUGH_SUPPORT] += 1
        if deduced_coding1 is None:
            return None
        return {deduced_coding1.A: 'A', deduced_coding1.B: 'B'}

    def _map_to_ab(self, call, coding):
        if not call.called:
            return None
        try:
            mapped = [coding[allele] for allele in call.gt_alleles]
        except KeyError:
            msg = 'An allele is not found in the deduced AB coding: '
            msg += str(allele)
            raise RuntimeError(msg)
        return mapped

    def _recode_parent_genotypes(self, samples=None):
        get_coding = GetCoding(self.parents_a, self.parents_b)

        def mapper(snp):
            return snp, get_coding(snp)

        win = self.window
        snp_and_coding = RandomAccessIterator(imap(mapper, self._reader),
                                              rnd_access_win=win)
        offspring = self.offspring
        half_win = (win - 1) // 2
        for idx, (snp1, coding1) in enumerate(snp_and_coding):
            snp1_calls = [snp1.genotype(sample) for sample in offspring]

            start = idx - half_win
            if start < 0:
                start = 0
            end = idx + half_win + 1

            snp2_idxs = []
            for snp2_idx in range(start, end):
                try:
                    snp2_chrom = snp_and_coding[snp2_idx][0].CHROM
                except IndexError:
                    continue
                if snp2_chrom == snp1.CHROM:
                    snp2_idxs.append(snp2_idx)

            coding1 = self._deduce_coding(snp_and_coding, snp1_calls,
                                          snp2_idxs)
            if coding1 is None:
                # We haven't manage to deduce the AB coding for this snp
                continue
            coding1['.'] = '.'
            if samples is None:
                calls = snp1.samples
            else:
                calls = [snp1.genotype(sample) for sample in samples]
            recoded = OrderedDict((call.sample, self._map_to_ab(call, coding1)) for call in calls)
            yield snp1, recoded

    def _count_recombinations(self, sample_genotypes):
        last_gt = None
        recombs = 0
        for gt in sample_genotypes:
            if last_gt is None:
                last_gt = gt
                continue
            if gt is None:
                continue
            if last_gt != gt:
                recombs += 1
                last_gt = gt
        return recombs

    def _smooth(self, snp_idx_to_smooth, snp_gts, samples):
        snp1 = snp_gts[snp_idx_to_smooth][0]     
        snp1_calls = [snp1.genotype(sample) for sample in samples]
        chrom = snp1.CHROM

        # remove snps from other chromosomes
        snp_gts = [snp_gt for snp_gt in snp_gts if snp_gt[0].CHROM == chrom]

        snps, gt_for_snps_in_win = zip(*snp_gts)

        # we need the recomb rates
        # TODO, this is already calculated in ab coding, we have to cache
        weights = []
        for snp2 in snps:
            snp2_calls = [snp2.genotype(sample) for sample in samples]
            haplos = _count_biallelic_haplotypes(snp1_calls, snp2_calls,
                                                 return_alleles=True)
            if haplos is None:
                weight = 0
            else:
                haplo_cnt, alleles = haplos
                recomb_rate = (haplo_cnt.aB + haplo_cnt.Ab) / sum(haplo_cnt)
                weight = 2 * (0.5 - recomb_rate) if recomb_rate < 0.5 else 0
            weights.append(weight)

        # we have to transpose, we want the genotype for each indi not for
        # each snp
        indis_gts = {indi: [] for indi in samples}
        for indi in samples:
            for snp_gt in gt_for_snps_in_win:
                snp_indi_gt =  snp_gt.get(indi, None)
                if snp_indi_gt is not None:
                    snp_indi_gt = tuple(sorted(snp_indi_gt))
                indis_gts[indi].append(snp_indi_gt)

        # Now we can do the smoothing
        recomb_thres = self.recomb_threshold
        smooth_threhsold = self.smooth_threhsold
        smoothed_genos = []
        for indi in samples:
            indi_gt = indis_gts[indi]
            n_recombs = self._count_recombinations(indi_gt)
            counts = Counter()
            for weight, geno in zip(weights, indi_gt):
                counts[geno] += weight
            if None in counts:
                del counts[None]
            if counts:
                smoothed_geno, vote = counts.most_common(1)[0]
            else:
                smoothed_geno, vote = None, 0
            tot_index = sum(counts.values())
            if tot_index:
                index = vote / tot_index
            else:
                index = 0
            self._recombs.append(n_recombs)
            self._smoothes.append(index)
            if recomb_thres is None:
                if index > smooth_threhsold:
                    geno = smoothed_geno
                else:
                    geno = None
            else:
                if n_recombs > recomb_thres:
                    # We're assuming diploid here
                    geno = ('A', 'B')
                else:
                    if index > smooth_threhsold:
                        geno = smoothed_geno
                    else:
                        geno = None
            smoothed_genos.append(geno)
        return smoothed_genos

    def _smooth_genotypes(self, snp_ab_genotypes, samples):
        win = self.window
        snp_ab_genotypes = RandomAccessIterator(snp_ab_genotypes,
                                                rnd_access_win=win)

        half_win = (win - 1) // 2
        for idx, (snp, ab_genotype) in enumerate(snp_ab_genotypes):
            chrom = snp.CHROM

            start = idx - half_win
            if start < 0:
                start = 0
            end = idx + half_win + 1

            snp_gts_in_win = snp_ab_genotypes[start: end]
            smoothed_genos = self._smooth(idx - start, snp_gts_in_win, samples)
            smoothed_genos = OrderedDict(zip(samples, smoothed_genos))
            yield snp, smoothed_genos

    def recode_genotypes(self, samples=None):
        snps_ab_coding = self._recode_parent_genotypes(samples=samples)
        if self.smooth_threhsold is not None:
            snps_ab_coding = self._smooth_genotypes(snps_ab_coding,
                                                    samples=samples)
        return snps_ab_coding

    def plot_parent_coding_hist(self, fhand):
        fig = Figure()
        axes = fig.add_subplot(111)
        axes.hist(self.indexes, fill=True, log=True, bins=20, rwidth=1)
        axes.axvline(x=self.parent_index_threshold)
        axes.set_xlabel('support index for parental coding')
        axes.set_ylabel('num. SNPs')
        _print_figure(axes, fig, fhand, plot_legend=False)

    def plot_smooth_hist(self, fhand):
        bins = 20
        fig = Figure()
        axes2 = fig.add_subplot(111)
        x = self._smoothes
        y = self._recombs
        image = axes2.hist2d(x, y, bins=bins,
                             norm=matplotlib.colors.LogNorm())[-1]

        axes2.tick_params(
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off',
        right='off', left='off', labelleft='off')
        xmin2d, xmax2d = axes2.get_xlim()
        ymin2d, ymax2d = axes2.get_ylim()
        axes2.vlines(self.smooth_threhsold, ymin2d, ymax2d, color='r')
        if self.recomb_threshold is not None:
            axes2.hlines(self.recomb_threshold, xmin2d, xmax2d, color='r')

        divider = make_axes_locatable(axes2)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(image, cax=cax)

        axes = divider.append_axes('bottom', size=2, pad=0.1, sharex=axes2)

        #axes = fig.add_subplot(224)
        #print axes2.get_position().bounds
        axes.hist(x, fill=True, log=True, bins=bins, rwidth=1)
        axes.set_xlabel('Smooth index')
        axes.set_ylabel('Num. SNPs')
        axes.set_xlim(xmin2d, xmax2d)
        ymin, ymax = axes.get_ylim()
        axes.vlines(self.smooth_threhsold, ymin, ymax, color='r')

        axes = divider.append_axes('left', size=2, pad=0.1, sharey=axes2)
        #axes = fig.add_subplot(221)
        axes.hist(y, orientation='horizontal', fill=True, log=True, bins=bins,
                  rwidth=1)
        axes.set_ylabel('Num. recombs.')
        axes.set_xlabel('Num. SNPs')
        _print_figure(axes, fig, fhand, plot_legend=False)
        axes.set_ylim(ymin2d, ymax2d)
        xmin, xmax = axes.get_xlim()
        if self.recomb_threshold is not None:
            axes.hlines(self.recomb_threshold, xmin, xmax, color='r')

    def write_log(self, fhand):
        tot = sum(self.log.values())
        fhand.write('%d SNPs divided in:\n' % tot)
        enough = self.log[ENOUGH_SUPPORT]
        fhand.write('Enough support: %d (%.1f%%)\n' % (enough,
                                                       enough / tot * 100))
        not_enough = self.log[NOT_ENOUGH_SUPPORT]
        fhand.write('Not enough support: %d (%.1f%%)\n' % (not_enough,
                                                           not_enough / tot * 100))
        no_info = self.log[NO_INFO]
        fhand.write('No information: %d (%.1f%%)\n' % (no_info,
                                                       no_info / tot * 100))
