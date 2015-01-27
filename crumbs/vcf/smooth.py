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

from collections import Counter, OrderedDict
from array import array

from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1 import make_axes_locatable

from crumbs.vcf.ab_coding import DEF_AB_CODING_WIN
from crumbs.iterutils import RandomAccessIterator
from crumbs.vcf.ld import _count_biallelic_haplotypes
from crumbs.vcf.filters import _print_figure

# Missing docstring
# pylint: disable=C0111


class Smoother(object):
    def __init__(self, smooth_threhsold, recomb_threshold=None,
                 window=DEF_AB_CODING_WIN):
        # TODO a min number of genotypes to evaluate anything
        self.window = window
        self.smooth_threhsold = smooth_threhsold
        self.recomb_threshold = recomb_threshold
        self._recombs = array('B')
        self._smoothes = array('f')

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
        snps, gt_for_snps_in_win = zip(*snp_gts)

        # we need the recomb rates
        # TODO, this is already calculated in ab coding, we have to cache
        snp1 = snps[snp_idx_to_smooth]
        snp1_calls = [snp1.genotype(sample) for sample in samples]
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
                snp_indi_gt =  snp_gt.get(indi, (indi, None))
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
            smoothed_geno, vote = counts.most_common(1)[0]
            index = vote / sum(counts.values())
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

    def smooth_genotypes(self, snp_ab_genotypes, samples):
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

            # remove snps in other chromosomes
            snp_gts_in_win = [snp_gt for snp_gt in snp_ab_genotypes[start: end] if snp_gt[0].CHROM == chrom]

            smoothed_genos = self._smooth(idx - start, snp_gts_in_win, samples)
            smoothed_genos = OrderedDict(zip(samples, smoothed_genos))
            yield snp, smoothed_genos

    def plot_hist(self, fhand):
        bins = 20
        fig = Figure()
        axes2 = fig.add_subplot(111)
        x = self._smoothes
        y = self._recombs
        image = axes2.hist2d(x, y, bins=bins)[-1]
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
