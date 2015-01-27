from __future__ import division

import random
from collections import OrderedDict
from math import isinf, isnan
from os.path import join as pjoin
from os.path import exists
from os import mkdir
from collections import namedtuple, Counter
import array
from StringIO import StringIO
from operator import itemgetter
from itertools import chain
import warnings

import numpy

from scipy.optimize import curve_fit
from scipy.stats.distributions import t

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

from crumbs.iterutils import group_in_packets
from crumbs.iterutils import RandomAccessIterator

from crumbs.vcf.snv import VCFReader, VCFWriter, DEF_MIN_CALLS_FOR_POP_STATS
from crumbs.vcf.ld import calc_recomb_rate
from crumbs.vcf import snv


# Missing docstring
# pylint: disable=C0111
# Too few pulic methods
# pylint: disable=R0903

PASSED = 'passed'
FILTERED_OUT = 'filtered_out'

SNPS_PER_FILTER_PACKET = 50
DEF_SNV_WIN = 101
DEF_MIN_NUM_CHECK_SNPS_IN_WIN = 50
DEF_MAX_FAILED_FREQ = 0.5
DEF_MAX_DIST = 1000000
DEF_MIN_DIST = 125000
DEF_ALPHA = 0.05
DEF_MAX_ZERO_DIST_RECOMB = 0.01
DEF_NUM_SNPS_IN_WIN_FOR_WEIRD_RECOMB = 51
DEF_MIN_NUM_SNPS_WEIRD_RECOMB = 20
DEF_MAX_RECOMB_RATE_WEIRD_RECOMB = 0.25


def group_in_filter_packets(items, items_per_packet):
    for packet in group_in_packets(items, items_per_packet):
        yield {PASSED: packet, FILTERED_OUT: []}


def _write_log(log_fhand, tot_snps, passed_snps):
    tot_snps = int(tot_snps)
    log_fhand.write('SNVs processed: ' + str(tot_snps) + '\n')
    good_snps = int(sum(passed_snps.values()))
    msg = 'SNVs passsed: ' + str(good_snps) + '\n'
    log_fhand.write(msg)
    msg = 'SNVs filtered out: ' + str(int(tot_snps - good_snps)) + '\n'
    log_fhand.write(msg)
    log_fhand.write('Number of SNVs that passed each filter\n')
    for filter_, count in passed_snps.items():
        msg = filter_ + ': ' + str(int(count)) + '\n'
        log_fhand.write(msg)
    log_fhand.flush()


def filter_snvs(in_fhand, out_fhand, filters, filtered_fhand=None,
                log_fhand=None, reader_kwargs=None):
    '''IT filters an input vcf.

    The input fhand has to be uncompressed. The original file could be a
    gzipped file, but in that case it has to be opened with gzip.open before
    sending it to this function.
    '''
    if reader_kwargs is None:
        reader_kwargs = {}
    # The input fhand to this function cannot be compressed
    reader_kwargs.update({'compressed': False,
                         'filename': 'pyvcf_bug_workaround'})

    reader = VCFReader(in_fhand, **reader_kwargs)

    template_reader = VCFReader(StringIO(reader.header))
    writer = VCFWriter(out_fhand, template_reader=template_reader)
    if filtered_fhand:
        filtered_writer = VCFWriter(filtered_fhand,
                                    template_reader=template_reader)
    else:
        filtered_writer = None

    packets = group_in_filter_packets(reader.parse_snvs(),
                                      SNPS_PER_FILTER_PACKET)
    tot_snps = 00.01
    passed_snps = OrderedDict()
    broken_pipe = False
    for packet in packets:
        tot_snps += len(packet[PASSED]) + len(packet[FILTERED_OUT])
        for filter_ in filters:
            packet = filter_(packet)
            filter_name = filter_.__class__.__name__
            if filter_name not in passed_snps:
                passed_snps[filter_name] = 0
            passed_snps[filter_name] += len(packet[PASSED])

        for snv in packet[PASSED]:
            if not _safe_write_snv(writer, snv):
                broken_pipe = True
                break
        if filtered_writer:
            for snv in packet[FILTERED_OUT]:
                if not _safe_write_snv(filtered_writer, snv):
                    broken_pipe = True
                    break
        if broken_pipe:
            break

    if log_fhand:
        _write_log(log_fhand, tot_snps, passed_snps)

    writer.flush()


def _safe_write_snv(writer, snv):
    try:
        writer.write_snv(snv)
        return True
    except IOError, error:
        # The pipe could be already closed
        if 'Broken pipe' not in str(error):
            raise
    return False


class _BaseFilter(object):
    def __init__(self, samples_to_consider=None, reverse=False):
        self.reverse = reverse
        self.samples_to_consider = samples_to_consider

    def _setup_checks(self, filterpacket):
        pass

    def _do_check(self, snv):
        raise NotImplementedError()

    def __call__(self, filterpacket):
        self._setup_checks(filterpacket)
        reverse = self.reverse
        items_passed = []
        filtered_out = filterpacket[FILTERED_OUT][:]
        samples_to_consider = self.samples_to_consider
        for snv in filterpacket[PASSED]:
            if samples_to_consider is not None:
                snv_to_check = snv.filter_calls_by_sample(samples=samples_to_consider,
                                                          reverse=False)
            else:
                snv_to_check = snv

            passed = self._do_check(snv_to_check)
            if reverse:
                passed = not passed
            if passed:
                items_passed.append(snv)
            else:
                filtered_out.append(snv)

        return {PASSED: items_passed, FILTERED_OUT: filtered_out}


class MonomorphicFilter(_BaseFilter):
    "Filters monomorphic snvs"
    def __init__(self, reverse=False, samples_to_consider=None,
                 freq_threshold=1):
        parent_kwargs = {'reverse': reverse,
                         'samples_to_consider': samples_to_consider}
        super(MonomorphicFilter, self).__init__(**parent_kwargs)
        self._freq_threslhold = freq_threshold

    def _do_check(self, snv):
        return snv.is_polymorphic(self._freq_threslhold)


class CallRateFilter(_BaseFilter):
    'Filter by the min. number of genotypes called'

    def __init__(self, min_calls=None, min_call_rate=None, reverse=False,
                 samples_to_consider=None):
        parent_kwargs = {'reverse': reverse,
                         'samples_to_consider': samples_to_consider}
        super(CallRateFilter, self).__init__(**parent_kwargs)
        if min_calls is not None and min_call_rate is not None:
            msg = 'Both min_calls and min_call rate cannot be given'
            msg += 'at the same time'
            raise ValueError(msg)
        elif min_calls is None and min_call_rate is None:
            msg = 'min_calls or call_rate should be given'
            raise ValueError(msg)
        self.min_calls = min_calls
        self.min_call_rate = min_call_rate

    def _do_check(self, snv):
        if self.min_calls:
            if snv.num_called >= self.min_calls:
                return True
            else:
                return False
        else:
            if snv.call_rate >= self.min_call_rate:
                return True
            else:
                return False


class BiallelicFilter(_BaseFilter):
    'Filter the biallelic SNPs'

    def __init__(self, reverse=False, samples_to_consider=None):
        parent_kwargs = {'reverse': reverse,
                         'samples_to_consider': samples_to_consider}
        super(BiallelicFilter, self).__init__(**parent_kwargs)

    def _do_check(self, snv):
        if len(snv.alleles) == 2:
            return True
        else:
            return False


class IsSNPFilter(_BaseFilter):
    def _do_check(self, snv):
        return snv.is_snp


class SnvQualFilter(_BaseFilter):
    def __init__(self, min_qual, reverse=False, samples_to_consider=None):
        parent_kwargs = {'reverse': reverse,
                         'samples_to_consider': samples_to_consider}
        super(SnvQualFilter, self).__init__(**parent_kwargs)
        self.min_qual = min_qual

    def _do_check(self, snv):
        qual = snv.qual
        if qual is None:
            return False
        else:
            return qual >= self.min_qual


class ObsHetFilter(_BaseFilter):
    def __init__(self, min_het=None, max_het=None, remove_nd=True,
                 reverse=False, samples_to_consider=None):
        parent_kwargs = {'reverse': reverse,
                         'samples_to_consider': samples_to_consider}
        super(ObsHetFilter, self).__init__(**parent_kwargs)
        self.min_het = min_het
        self.max_het = max_het
        self.remove_nd = remove_nd

    def _do_check(self, snv):
        min_het = self.min_het
        max_het = self.max_het
        het = snv.obs_het
        if het is None and self.remove_nd:
            return False
        if min_het is not None and het < min_het:
            return False
        if max_het is not None and het > max_het:
            return False
        return True


class MafFilter(_BaseFilter):
    def __init__(self, min_maf=None, max_maf=None, remove_nd=True,
                 reverse=False, samples_to_consider=None):
        parent_kwargs = {'reverse': reverse,
                         'samples_to_consider': samples_to_consider}
        super(MafFilter, self).__init__(**parent_kwargs)
        if min_maf is None and max_maf is None:
            msg = 'At least one value should be given for the min or max het'
            raise ValueError(msg)
        self.min_maf = min_maf
        self.max_maf = max_maf
        self.remove_nd = remove_nd

    def _do_check(self, snv):
        min_maf = self.min_maf
        max_maf = self.max_maf
        maf = snv.maf
        if maf is None and self.remove_nd:
            return False
        if min_maf is not None and maf < min_maf:
            return False
        if max_maf is not None and maf > max_maf:
            return False
        return True

FISHER_CACHE = {}


def _fisher_extact_rxc(counts_obs, counts_exp):
    if (counts_obs, counts_exp) in FISHER_CACHE:
        return FISHER_CACHE[(counts_obs, counts_exp)]
    import rpy2.robjects as robjects
    env = robjects.r.baseenv()
    env['obs'] = robjects.IntVector(counts_obs)
    env['expected'] = robjects.IntVector(counts_exp)
    pvalue = robjects.r('fisher.test(cbind(obs, expected))$p.value')[0]

    FISHER_CACHE[(counts_obs, counts_exp)] = pvalue
    return pvalue


class WeirdSegregationFilter(object):
    def __init__(self, alpha=DEF_ALPHA, num_snvs_check=DEF_SNV_WIN,
                 max_failed_freq=DEF_MAX_FAILED_FREQ, samples=None,
                 win_width=DEF_MAX_DIST, win_mask_width=DEF_MIN_DIST,
                 min_num_snvs_check_in_win=DEF_MIN_NUM_CHECK_SNPS_IN_WIN,
                 debug_plot_dir=None):
        # We're assuming that most snps are ok and that a few have a weird
        # segregation
        self.alpha = alpha
        self.num_snvs_check = num_snvs_check
        self.max_failed_freq = max_failed_freq
        self.win_width = win_width
        if win_mask_width < 1:
            msg = 'You should mask at least a window of 3 bp to avoid the '
            msg += 'comparison with itself'
            raise ValueError(msg)
        self.win_mask_width = win_mask_width
        self.min_num_snvs_check_in_win = min_num_snvs_check_in_win
        self._failed_freqs = array.array('f')
        self.tot_snps = 0
        self.passed_snps = 0
        self.samples = samples
        if debug_plot_dir is not None and not exists(debug_plot_dir):
            mkdir(debug_plot_dir)
        self.plot_dir = debug_plot_dir

    def filter_vcf(self, vcf_fpath, min_samples=DEF_MIN_CALLS_FOR_POP_STATS):
        reader = VCFReader(open(vcf_fpath),
                           min_calls_for_pop_stats=min_samples)
        snvs = reader.parse_snvs()
        random_reader = VCFReader(open(vcf_fpath))

        for snv_1 in snvs:
            self.tot_snps += 1
            loc = snv_1.pos

            if self.plot_dir:
                chrom = str(snv_1.chrom)
                fname = chrom + '_' + str(loc) + '.png'
                chrom_dir = pjoin(self.plot_dir, chrom)
                if not exists(chrom_dir):
                    mkdir(chrom_dir)
                plot_fhand = open(pjoin(chrom_dir, fname), 'w')
                debug_plot_info = []
            else:
                plot_fhand = None

            win_1_start = loc - (self.win_width / 2)
            if win_1_start < 0:
                win_1_start = 0
            win_1_end = loc - (self.win_mask_width / 2)
            if win_1_end < 0:
                win_1_end = 0
            if win_1_end != 0:
                snvs_win_1 = random_reader.fetch_snvs(snv_1.chrom,
                                                      start=int(win_1_start),
                                                      end=int(win_1_end))
            else:
                snvs_win_1 = []

            win_2_start = loc + (self.win_mask_width / 2)
            win_2_end = loc + (self.win_width / 2)
            snvs_win_2 = random_reader.fetch_snvs(snv_1.chrom,
                                                  start=win_2_start,
                                                  end=win_2_end)
            snvs_in_win = list(chain(snvs_win_1, snvs_win_2))
            if len(snvs_in_win) > self.num_snvs_check:
                snvs_in_win = random.sample(snvs_in_win, self.num_snvs_check)
            if len(snvs_in_win) < self.min_num_snvs_check_in_win:
                # Not enough snps to check
                continue

            orig_snp = snv_1
            if self.samples is not None:
                snv_1 = snv_1.filter_calls_by_sample(self.samples)

            exp_cnts = snv_1.biallelic_genotype_counts

            if exp_cnts is None:
                continue

            test_values = []
            for snv_2 in snvs_in_win:
                if self.samples is not None:
                    snv_2 = snv_2.filter_calls_by_sample(self.samples)
                obs_cnts = snv_2.biallelic_genotype_counts
                if obs_cnts is None:
                    continue
                test_values.append(_fisher_extact_rxc(obs_cnts, exp_cnts))

                if plot_fhand:
                    debug_plot_info.append({'pos': snv_2.pos,
                                            'AA': obs_cnts[0],
                                            'Aa': obs_cnts[1],
                                            'aa': obs_cnts[2],
                                            'close_snp': True})
            alpha2 = self.alpha/len(test_values)
            results = []
            for idx, val in enumerate(test_values):
                result = False if val is None else val > alpha2
                results.append(result)

                if plot_fhand:
                    debug_plot_info[idx]['result'] = result

            if len(test_values) < self.min_num_snvs_check_in_win:
                # few snps can be tested for segregation
                continue

            tot_checked = len(test_values)
            if tot_checked > 0:
                failed_freq = results.count(False) / tot_checked
                passed = self.max_failed_freq > failed_freq
            else:
                failed_freq = None
                passed = False
            if failed_freq is not None:
                self._failed_freqs.append(failed_freq)

            if plot_fhand:
                debug_plot_info.append({'pos': snv_1.pos,
                                        'AA': exp_cnts[0],
                                        'Aa': exp_cnts[1],
                                        'aa': exp_cnts[2],
                                        'result': passed,
                                        'close_snp': False})
                self._plot_segregation_debug(debug_plot_info, plot_fhand)
            if passed:
                self.passed_snps += 1
                yield orig_snp

    @staticmethod
    def _plot_segregation_debug(plot_info, fhand):
        greens = ['#2d6e12', '#76b75b', '#164900']
        reds = ['#8f0007', '#fb1f2a', '#b90009']
        grays = ['0.5', '0.25', '0.75']
        fig = Figure(figsize=(10, 4))
        axes1 = fig.add_subplot(211)
        axes2 = fig.add_subplot(212)
        plot_info.sort(key=itemgetter('pos'))

        is_close_snp = [geno_info['close_snp'] for geno_info in plot_info]
        tested_snp_idx = is_close_snp.index(False)
        passed = [geno_info['result'] for geno_info in plot_info]
        total_cnts = [geno_info['AA'] + geno_info['Aa'] + geno_info['aa'] for geno_info in plot_info]

        num_snvs = len(passed)
        bottoms = [0] * num_snvs
        freq_bottoms = [0] * num_snvs
        lefts = range(num_snvs)
        for idx, geno in enumerate(('AA', 'Aa', 'aa')):
            cnts = [geno_info[geno] for geno_info in plot_info]
            backcolors = [greens[idx] if pss else reds[idx] for pss in passed]
            edgecolors = [grays[idx] for pss in passed]
            backcolors[tested_snp_idx], edgecolors[tested_snp_idx] = \
                edgecolors[tested_snp_idx], backcolors[tested_snp_idx]
            heights = cnts
            axes1.bar(lefts, heights, bottom=bottoms, color=backcolors,
                      edgecolor=edgecolors)

            freq_heights = [height/total_cnt for height, total_cnt in zip(heights, total_cnts)]
            axes2.bar(lefts, freq_heights, bottom=freq_bottoms,
                      color=backcolors, edgecolor=edgecolors)

            bottoms = [bot + heig for bot, heig in zip(bottoms, heights)]
            freq_bottoms = [bot + heig for bot, heig in zip(freq_bottoms, freq_heights)]

        axes1.tick_params(axis='y', which='both', left='off', right='off')
        axes2.tick_params(axis='y', which='both', left='off', right='off')

        canvas = FigureCanvas(fig)
        canvas.print_figure(fhand)
        fhand.flush()

    def plot_failed_freq_dist(self, fhand):
        fig = Figure()
        axes = fig.add_subplot(111)
        axes.hist(self._failed_freqs, fill=True, log=True, bins=20,
                  rwidth=1)
        axes.axvline(x=self.max_failed_freq)
        axes.set_xlabel('% of adjacent SNPs segregating differently')
        axes.set_ylabel('num. SNPs')
        _print_figure(axes, fig, fhand, plot_legend=False)

    def write_log(self, fhand):
        _write_log(fhand, self.tot_snps,
                   {self.__class__.__name__: self.passed_snps})


def _get_calls(snv, samples):
    if samples is None:
        calls = snv.calls
    else:
        calls = [call for call in snv.calls if call.sample in samples]
    return calls


RecombRate = namedtuple('RecombRate', ['index_in_vcf', 'pos', 'recomb_rate'])


def _calculate_segregation_rates(snvs, pop_type, snps_in_window,
                                 samples=None):
    half_win = (snps_in_window - 1) // 2
    prev_chrom = None
    for index1, snp1 in enumerate(snvs):
        calls1 = _get_calls(snp1, samples)
        start = index1 - half_win
        if start < 0:
            start = 0
        rates = []
        chrom = snp1.chrom
        if chrom != prev_chrom:
            recomb_cache = {}
            prev_chrom = chrom
        for index2 in range(start, index1 + half_win):
            try:
                snp2 = snvs[index2]
            except IndexError:
                continue
            if chrom != snp2.chrom:
                continue
            calls2 = _get_calls(snp2, samples)
            index = tuple(sorted([index1, index2]))
            if index1 == index2:
                recomb_rate = 0
            try:
                recomb_rate = recomb_cache[index]
            except KeyError:
                recomb_rate = _calc_recomb_rate(calls1, calls2, pop_type)
                if recomb_rate is None:
                    recomb_rate = float('nan')
                else:
                    recomb_rate = recomb_rate[0]
                recomb_cache[index] = recomb_rate
            rates.append(RecombRate(index2, snp2.pos, recomb_rate))
        yield snp1, chrom, snp1.pos, rates


class WeirdRecombFilter(object):
    def __init__(self, pop_type, max_zero_dist_recomb=DEF_MAX_ZERO_DIST_RECOMB,
                 alpha_recomb_0=DEF_ALPHA,
                 snps_in_window=DEF_NUM_SNPS_IN_WIN_FOR_WEIRD_RECOMB,
                 min_num_snps=DEF_MIN_NUM_SNPS_WEIRD_RECOMB,
                 max_recomb_curve_fit=DEF_MAX_RECOMB_RATE_WEIRD_RECOMB,
                 debug_plot_dir=None, samples=None):
        self.pop_type = pop_type
        self.max_zero_dist_recomb = max_zero_dist_recomb
        self.alpha_recomb_0 = alpha_recomb_0
        self.snps_in_window = snps_in_window
        self.min_num_snps = min_num_snps
        self.max_recomb_curve_fit = max_recomb_curve_fit
        self.debug_plot_dir = debug_plot_dir
        self.samples = samples
        self.tot_snps = 0
        self.passed_snps = 0
        self.not_fitted_counter = Counter()
        self.recomb_rates = {'ok': array.array('f'),
                             'ok_conf_is_None': array.array('f'),
                             'not_ok': array.array('f')}

    def filter_snvs(self, snvs):

        snps = RandomAccessIterator(snvs, rnd_access_win=self.snps_in_window)
        rates = _calculate_segregation_rates(snps, self.pop_type,
                                             self.snps_in_window,
                                             samples=self.samples)
        max_zero_dist = self.max_zero_dist_recomb
        for snp, chrom, pos, rates in rates:
            self.tot_snps += 1
            dists, recombs = zip(*[(rate.pos - pos, rate.recomb_rate) for rate in rates])
            if len(dists) < self.min_num_snps:
                continue
            if self.debug_plot_dir is None:
                plot_fhand = None
            else:
                chrom_dir = pjoin(self.debug_plot_dir, str(chrom))
                if not exists(chrom_dir):
                    mkdir(chrom_dir)
                fname = str(chrom) + '_' + str(pos) + '.png'
                plot_fhand = open(pjoin(chrom_dir, fname), 'w')
            res = _calc_ajusted_recomb(dists, recombs,
                                       max_recomb=self.max_recomb_curve_fit,
                                       max_zero_dist_recomb=max_zero_dist,
                                       alpha_recomb_0=self.alpha_recomb_0,
                                       plot_fhand=plot_fhand)
            self._store_log_info(*res)
            if res[1]:
                self.passed_snps += 1
                yield snp

    def _store_log_info(self, recomb_at_0, snp_ok, debug_info):
        if 'reason_no_fit' in debug_info:
            self.not_fitted_counter[debug_info['reason_no_fit']] += 1
            return
        index = snp_ok, debug_info['conf_interval'] is None
        if snp_ok:
            if debug_info['conf_interval'] is None:
                index = 'ok_conf_is_None'
            else:
                index = 'ok'
        else:
            index = 'not_ok'
        self.recomb_rates[index].append(recomb_at_0)

    def plot_recomb_at_0_dist_hist(self, fhand):
        fig = Figure()
        axes = fig.add_subplot(111)
        data = [self.recomb_rates['ok'], self.recomb_rates['ok_conf_is_None'],
                self.recomb_rates['not_ok']]
        labels = ['OK', 'OK conf. is none', 'Not OK']
        colors = [(0.3, 1, 0.3), (0.3, 1, 0.6), (1, 0.3, 0.3)]
        some_data = [bool(dataset) for dataset in data]
        labels = [label for draw, label in zip(some_data, labels) if draw]
        colors = [color for draw, color in zip(some_data, colors) if draw]
        data = [dataset for draw, dataset in zip(some_data, data) if draw]
        axes.hist(data, stacked=True, fill=True, log=True, bins=20,
                  label=labels, rwidth=1, color=colors)
        _print_figure(axes, fig, fhand)

    def write_log(self, fhand):
        _write_log(fhand, self.tot_snps,
                   {self.__class__.__name__: self.passed_snps})


def _kosambi(phys_dist, phys_gen_dist_conversion, recomb_at_origin):
    phys_gen_dist_conversion = abs(phys_gen_dist_conversion)
    # recomb rate should be in morgans per base
    d4 = numpy.absolute(phys_dist) * phys_gen_dist_conversion * 4
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try:
            ed4 = numpy.exp(d4)
        except Warning:
            raise RuntimeError('Numpy raised a warning calculating exp')
    return 0.5 * (ed4 - 1) / (ed4 + 1) + recomb_at_origin


def _fit_kosambi(dists, recombs, init_params):
    try:
        return curve_fit(_kosambi, dists, recombs, p0=init_params)
    except RuntimeError:
        return None, None
    except TypeError:
        # It happens when recombs is all nan
        return None, None


def _print_figure(axes, figure, plot_fhand, plot_legend=True):
    if figure is None:
        return
    if plot_legend:
        axes.legend()
    canvas = FigureCanvas(figure)
    canvas.print_figure(plot_fhand)
    plot_fhand.flush()


def _calc_ajusted_recomb(dists, recombs, max_recomb, max_zero_dist_recomb,
                         alpha_recomb_0, plot_fhand=None):
    # first rough interpolation
    # we remove the physical distances with high recombination rates because
    # they're not very informative. e.g. more than 40 cM will not discriminate
    # between false recombination due to hidden segregation in the parents and
    # true recombination

    if plot_fhand:
        fig = Figure()
        axes = fig.add_subplot(111)
        axes.set_axis_bgcolor((1, 0.6, 0.6))
        axes.scatter(dists, recombs, c='r', label='For 1st fit')
    else:
        axes = None
        fig = None

    dists = numpy.array(dists)
    recombs = numpy.array(recombs)
    recomb_rate = 1e-7
    popt, pcov = _fit_kosambi(dists, recombs, init_params=[recomb_rate, 0])
    if popt is None:
        _print_figure(axes, fig, plot_fhand)
        return None, False, {'kosambi_fit_ok': False,
                             'reason_no_fit': '1st fit failed'}

    est_dists = dists
    est_recombs = _kosambi(est_dists, popt[0], popt[1])

    if fig:
        axes.plot(est_dists, est_recombs, label='1st fit', c='r')

    # now we perform a second fit but only with those markers that are a
    # distance that results in a recombination fraction lower than max_recomb
    close_markers = est_recombs < max_recomb
    close_recombs = recombs[close_markers]
    close_dists = dists[close_markers]

    if plot_fhand:
        axes.scatter(close_dists, close_recombs, c='b', label='For 2nd fit')

    if len(close_dists) < 1:
        # This marker is so bad that their closest markers are at a large
        # distance
        _print_figure(axes, fig, plot_fhand)
        return None, False, {'kosambi_fit_ok': False,
                             'reason_no_fit': 'no close region left'}

    if len(close_dists) != len(dists):
        # If we've removed any points we fit again
        popt, pcov = _fit_kosambi(close_dists, close_recombs, init_params=popt)
    if popt is None:
        _print_figure(axes, fig, plot_fhand)
        return None, False, {'kosambi_fit_ok': False,
                             'reason_no_fit': '2nd fit failed'}

    est_close_recombs = _kosambi(close_dists, popt[0], popt[1])

    residuals = close_recombs - est_close_recombs
    if fig:
        axes.plot(close_dists, est_close_recombs, c='b', label='2nd_fit')

    # we exclude the markers with a residual outlier
    quartile_25, quartile_75 = numpy.percentile(residuals, [25, 75])
    iqr = quartile_75 - quartile_25
    outlayer_thrld = [quartile_25 - iqr * 1.5, quartile_75 + iqr * 1.5]
    ok_markers = [idx for idx, res in enumerate(residuals) if (not isnan(res) and (outlayer_thrld[0] < res < outlayer_thrld[1]))]
    ok_recombs = close_recombs[ok_markers]
    ok_dists = close_dists[ok_markers]

    if fig:
        axes.scatter(ok_dists, ok_recombs, c='g', label='For 3rd fit')

    if len(ok_dists) != len(close_dists):
        # If we've removed any points we fit again
        popt, pcov = _fit_kosambi(ok_dists, ok_recombs, init_params=popt)
    if popt is None:
        _print_figure(axes, fig, plot_fhand)
        return None, False, {'kosambi_fit_ok': False,
                             'reason_no_fit': '3rd fit failed'}
    var_recomb_at_dist_0 = pcov[1, 1]
    recomb_at_dist_0 = popt[1]
    ok_color = (0.3, 1, 0.6)
    if isinf(var_recomb_at_dist_0):
        conf_interval = None
        if abs(recomb_at_dist_0) < 0.01:
            # recomb is 0 for all points and the variance is inf
            snp_ok = True
        else:
            snp_ok = False
    else:
        if alpha_recomb_0 is None:
            conf_interval = None
            if abs(recomb_at_dist_0) <= max_zero_dist_recomb:
                snp_ok = True
                ok_color = (0.3, 1, 0.3)
            else:
                snp_ok = False
        else:
            num_data_points = len(ok_dists)
            num_params = len(popt)
            deg_of_freedom = max(0, num_data_points - num_params)
            tval = t.ppf(1.0 - alpha_recomb_0 / 2., deg_of_freedom)
            std_dev = var_recomb_at_dist_0 ** 0.5
            conf_interval = (recomb_at_dist_0 - std_dev * tval,
                             recomb_at_dist_0 + std_dev * tval)

            if abs(recomb_at_dist_0) <= max_zero_dist_recomb:
                snp_ok = True
                ok_color = (0.3, 1, 0.3)
            elif conf_interval[0] < 0 < conf_interval[1]:
                snp_ok = True
            else:
                snp_ok = False
            if plot_fhand:
                axes.vlines(0, conf_interval[0], conf_interval[1],
                            label='conf. interval')

    if plot_fhand:
        color = ok_color if snp_ok else (1, 0.3, 0.3)
        axes.set_axis_bgcolor(color)

    if popt is None:
        _print_figure(axes, fig, plot_fhand)
        return None, False, {'kosambi_fit_ok': False,
                             'reason_no_fit': '3rd fit failed'}

    est2_recombs = _kosambi(ok_dists, popt[0], popt[1])

    if fig:
        axes.plot(ok_dists, est2_recombs, c='g', label='3rd_fit')
        _print_figure(axes, fig, plot_fhand)
    return recomb_at_dist_0, snp_ok, {'kosambi_fit_ok': True,
                                      'conf_interval': conf_interval}
