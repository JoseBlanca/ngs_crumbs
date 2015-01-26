from __future__ import division

import sys
from collections import Counter, OrderedDict, namedtuple
import gzip
from operator import itemgetter
from io import BytesIO

from vcf import Reader as pyvcfReader
from vcf import Writer as pyvcfWriter
from vcf.model import make_calldata_tuple
# ouch, _Call is a private class, but we don't know how to modify a Call
from vcf.model import _Call as pyvcfCall
from vcf.model import _Record as pyvcfRecord
from crumbs.iterutils import generate_windows

from crumbs.seq.seqio import read_seqs
from crumbs.seq.seq import get_name, get_length
from crumbs.utils.file_utils import flush_fhand

# Missing docstring
# pylint: disable=C0111

VARSCAN = 'VarScan'
GATK = 'gatk'
FREEBAYES = 'freebayes'
GENERIC = 'generic'
HOM_REF = 0
HET = 1
HOM_ALT = 2
HOM = 3
MISSING_ALLELE_CHAR = '.'

DEF_MIN_CALLS_FOR_POP_STATS = 10
DEF_MIN_NUM_SNPS_IN_WIN = 5

# TODO check if SNV can be converted in a proxy using the recipes
# http://code.activestate.com/recipes/496741-object-proxying/
# http://code.activestate.com/recipes/496742-shelfproxy/


class _SNPQueue(object):
    def __init__(self):
        self.queue = []

    def empty(self):
        self.queue = []

    def pop(self, location):
        queue = self.queue
        index_to_remove = None
        for index, snp in enumerate(queue):
            if snp.pos < location:
                index_to_remove = index
            else:
                break
        if index_to_remove is not None:
            del queue[:index_to_remove + 1]

    def extend(self, snps):
        self.queue.extend(snps)


class _SNPSlidingWindow(object):
    def __init__(self, snp_reader, win_size, win_step, min_num_snps,
                 ref_fhand=None):
        self._snp_queue = _SNPQueue()
        self._reader = snp_reader
        self.win_size = win_size
        self.win_step = win_step
        self.min_num_snps = min_num_snps
        self._ref_fhand = ref_fhand

    def _snp_in_window(self, snp, win):
        if win[0] <= snp.pos < win[1]:
            return True
        else:
            return False

    def _get_chrom_lengths(self):
        chrom_lens = OrderedDict()
        if self._ref_fhand is None:
            vcf_fhand = gzip.open(self._reader.fhand.name)
            for line in vcf_fhand:
                line = line.strip()
                if line.startswith('#'):
                    continue
                items = line.split()
                chrom = items[0]
                loc = int(items[1])
                if chrom not in chrom_lens:
                    chrom_lens[chrom] = loc
                else:
                    if loc > chrom_lens[chrom]:
                        chrom_lens[chrom] = loc

        else:
            for read in read_seqs([self._ref_fhand]):
                chrom_lens[get_name(read)] = get_length(read)
        return chrom_lens

    def windows(self):
        chrom_lengths = self._get_chrom_lengths()
        snp_queue = self._snp_queue
        for chrom, chrom_length in chrom_lengths.items():
            wins = generate_windows(start=0,
                                    size=self.win_size, step=self.win_step,
                                    end=chrom_length + 1)
            snp_queue.empty()
            for win in wins:
                snp_queue.pop(win.start)
                if snp_queue.queue:
                    new_strech_start = snp_queue.queue[-1].pos + 1
                else:
                    new_strech_start = win.start
                new_snps = self._reader.fetch_snvs(chrom, new_strech_start,
                                                   win.end)
                snp_queue.extend(new_snps)
                if len(snp_queue.queue) >= self.min_num_snps:
                    yield {'chrom': chrom, 'start': win.start, 'end': win.end,
                           'snps': snp_queue.queue[:]}


class VCFReader(object):
    def __init__(self, fhand, compressed=None, filename=None,
                 min_calls_for_pop_stats=DEF_MIN_CALLS_FOR_POP_STATS):
        self.fhand = fhand
        self.pyvcf_reader = pyvcfReader(fsock=fhand, compressed=compressed,
                                        filename=filename)
        self.min_calls_for_pop_stats = min_calls_for_pop_stats
        self._snpcaller = None
        self._samples = self.pyvcf_reader.samples
        self._header_lines = self.pyvcf_reader._header_lines
        self._column_headers = self.pyvcf_reader._column_headers

    def parse_snvs(self):
        min_calls_for_pop_stats = self.min_calls_for_pop_stats
        last_snp = None
        try:
            counter =0
            for snp in self.pyvcf_reader:
                counter +=1
                snp = SNV(snp, reader=self,
                          min_calls_for_pop_stats=min_calls_for_pop_stats)
                last_snp = snp
                yield snp
        except Exception:
            from traceback import print_exception
            exc_type, exc_value, exc_traceback = sys.exc_info()

            print_exception(exc_type, exc_value, exc_traceback,
                            limit=20, file=sys.stderr)

            if last_snp is not None:
                chrom = str(last_snp.chrom)
                pos = last_snp.pos
                msg = 'Last parsed SNP was: {} {}\n'.format(chrom, pos + 1)
                sys.stderr.write(msg)
                raise

    def fetch_snvs(self, chrom, start, end=None):
        min_calls_for_pop_stats = self.min_calls_for_pop_stats
        try:
            snvs = self.pyvcf_reader.fetch(chrom, start + 1, end=end)
        except KeyError:
            snvs = []
        if snvs is None:
            snvs = []

        for snp in snvs:
            snp = SNV(snp, reader=self,
                      min_calls_for_pop_stats=min_calls_for_pop_stats)
            yield snp

    def sliding_windows(self, size, step=None, ref_fhand=None,
                        min_num_snps=DEF_MIN_NUM_SNPS_IN_WIN):
        random_snp_reader = VCFReader(open(self.fhand.name))
        sliding_window = _SNPSlidingWindow(snp_reader=random_snp_reader,
                                           win_size=size, win_step=step,
                                           min_num_snps=min_num_snps,
                                           ref_fhand=ref_fhand)
        for window in sliding_window.windows():
            yield window

    @property
    def snpcaller(self):
        if self._snpcaller is not None:
            return self._snpcaller

        metadata = self.pyvcf_reader.metadata
        if 'source' in metadata:
            if 'VarScan2' in metadata['source']:
                snpcaller = VARSCAN
            elif 'freebayes' in metadata['source'][0].lower():
                snpcaller = FREEBAYES
            else:
                snpcaller = GENERIC
        elif 'UnifiedGenotyper' in metadata:
            snpcaller = GATK
        else:
            snpcaller = GENERIC
        self._snpcaller = snpcaller
        return snpcaller

    @property
    def samples(self):
        return self._samples

    @property
    def filters(self):
        return self.pyvcf_reader.filters

    @property
    def infos(self):
        return self.pyvcf_reader.infos

    def _build_header(self, samples=None):
        if samples is None:
            samples = self.samples
        header = '\n'.join(self._header_lines)
        header += '\n#' + '\t'.join(self._column_headers)
        header += '\t' + '\t'.join(samples)
        return header

    @property
    def header(self):
        return self._build_header()

    def create_template_reader(self, samples):
        vcf = self._build_header(samples)
        fake_vcf = BytesIO(vcf)
        return VCFReader(fake_vcf)


class VCFWriter(pyvcfWriter):

    def __init__(self, stream, template_reader, lineterminator="\n"):
        template = template_reader.pyvcf_reader
        super(VCFWriter, self).__init__(stream, template,
                                        lineterminator=lineterminator)

    def write_snv(self, snv):
        super(VCFWriter, self).write_record(snv.record)

    def write_snvs(self, snvs):
        for snv in snvs:
            try:
                self.write_snv(snv)
            except IOError, error:
                # The pipe could be already closed
                if 'Broken pipe' in str(error):
                    break
                else:
                    raise

    def flush(self):
        flush_fhand(self.stream)


BiallelicGts = namedtuple('BiallelicGts', ['AA', 'Aa', 'aa'])


def get_ids(snp):
    if hasattr(snp, 'calls'):
        ids = snp.record.ID
    else:
        ids = snp.ID
    if ids is None:
        return []
    return ids.split(';')


def get_or_create_id(snp, prefix=''):
    ids = get_ids(snp)
    if ids:
        return ids[0]
    else:
        if hasattr(snp, 'calls'):
            chrom = snp.chrom
            pos = snp.pos
        else:
            chrom = snp.CHROM
            pos = snp.POS - 1
        return prefix + chrom + '_' + str(pos + 1)


class SNV(object):
    '''A proxy around the pyvcf _Record with some additional functionality'''
    def __init__(self, record, reader, ploidy=2,
                 min_calls_for_pop_stats=DEF_MIN_CALLS_FOR_POP_STATS):
        # min_calls_for_pop_stats is defined here because otherwise it would
        # behave like a global variable for all SNPs read from a common
        # reader
        self.record = record
        self.reader = reader
        self.ploidy = ploidy
        self.min_calls_for_pop_stats = min_calls_for_pop_stats
        self._mac_analyzed = False
        self._maf = None
        self._mac = None
        self._allele_counts = None
        self._maf_dp_analyzed = False
        self._allele_depths = None
        self._maf_depth = None
        self._depth = None

    @property
    def ids(self):
        return get_ids(self)

    def get_or_create_id(self, prefix=''):
        return get_or_create_id(self, prefix)

    @property
    def alleles(self):
        alleles = []
        for allele in self.record.alleles:
            try:
                allele = allele.sequence
            except AttributeError:
                pass
            alleles.append(allele)
        return alleles

    @property
    def obs_het(self):
        snp = self.record
        n_called = snp.num_called
        if n_called >= self.min_calls_for_pop_stats:
            return snp.num_het / n_called
        else:
            return None

    @property
    def exp_het(self):
        snp = self.record
        if snp.num_called >= self.min_calls_for_pop_stats:
            return snp.heterozygosity
        else:
            return None

    @property
    def calls(self):
        return [Call(sample, snv=self) for sample in self.record.samples]

    def get_call(self, sample):
        call = self.record.genotype(sample)
        return Call(call, self)

    def _calculate_maf_and_mac(self):
        if self._mac_analyzed:
            return
        self._mac_analyzed = True
        ploidy = self.ploidy
        snp = self.record

        n_chroms_sampled = 0
        allele_counts = Counter()
        for call in snp.samples:
            if call.called:
                genotype = map(int, call.gt_alleles)
                assert len(genotype) == ploidy
                n_chroms_sampled += ploidy
                for allele in genotype:
                    allele_counts[allele] += 1
        if not n_chroms_sampled:
            return
        max_allele_count = max(allele_counts.values())
        self._allele_counts = allele_counts
        self._maf = max_allele_count / n_chroms_sampled
        self._mac = n_chroms_sampled - max_allele_count

    @property
    def allele_counts(self):
        if self.record.num_called >= self.min_calls_for_pop_stats:
            self._calculate_maf_and_mac()
        return self._allele_counts

    @property
    def genotype_counts(self):
        snp = self.record
        if snp.num_called < self.min_calls_for_pop_stats:
            return None
        return Counter([tuple(sorted(call.int_alleles)) for call in self.calls if call.called])

    @property
    def genotype_freqs(self):
        counts = self.genotype_counts
        if counts is None:
            return None
        tot_genos = sum(counts.values())
        return {geno: cnt/tot_genos for geno, cnt in counts.items()}

    @property
    def biallelic_genotype_counts(self):
        gt_cnts = self.genotype_counts
        if gt_cnts is None:
            return

        gt_cnts = {gt: cnt for gt, cnt in gt_cnts.items() if None not in gt}
        al_cnts = Counter()
        for genotype, cnt in gt_cnts.items():
            for allele in genotype:
                al_cnts[allele] += cnt
        sorted_als = list(sorted(al_cnts.items(), key=itemgetter(1),
                                 reverse=True))
        al_A = sorted_als[0][0]

        gt_cnts_aa = {}
        for genotype, cnt in gt_cnts.items():
            genotype = tuple(sorted(('A' if allele == al_A else 'a' for allele in genotype)))
            gt_cnts_aa[genotype] = cnt

        return BiallelicGts(gt_cnts_aa.get(('A', 'A'), 0),
                            gt_cnts_aa.get(('A', 'a'), 0),
                            gt_cnts_aa.get(('a', 'a'), 0))

    @property
    def biallelic_genotype_freqs(self):
        gt_cnts = self.biallelic_genotype_counts
        if gt_cnts is None:
            return
        tot = sum(gt_cnts)
        return BiallelicGts(*[cnt / tot for cnt in gt_cnts])

    @property
    def maf(self):
        'Frequency of the most abundant allele'
        snp = self.record
        if snp.num_called >= self.min_calls_for_pop_stats:
            self._calculate_maf_and_mac()
            return self._maf
        else:
            return None

    @property
    def mac(self):
        'Sum of the allele count of all alleles but the most abundant'
        self._calculate_maf_and_mac()
        return self._mac

    @property
    def inbreed_coef(self):
        obs_het = self.obs_het
        if obs_het is None:
            return None
        try:
            inbreed_f = 1 - (obs_het / self.exp_het)
            return inbreed_f
        except ZeroDivisionError:
            return None

    def _calculate_mafs_dp(self):
        if self._maf_dp_analyzed:
            return None
        self._maf_dp_analyzed = True

        allele_depths = Counter()
        for call in self.calls:
            if not call.has_alternative_counts:
                continue
            for allele, depth in call.allele_depths.items():
                allele_depths[allele] += depth
        depth = sum(allele_depths.values())
        if not depth:
            return None
        maf_dp = max(allele_depths.values()) / depth
        self._depth = depth
        self._maf_depth = maf_dp
        self._allele_depths = allele_depths

    @property
    def maf_depth(self):
        self._calculate_mafs_dp()
        return self._maf_depth

    @property
    def allele_depths(self):
        self._calculate_mafs_dp()
        return self._allele_depths

    @property
    def depth(self):
        self._calculate_mafs_dp()
        return self._depth

    @property
    def chrom(self):
        return self.record.CHROM

    @property
    def pos(self):
        return self.record.POS - 1

    @property
    def end(self):
        return self.record.end

    @property
    def qual(self):
        return self.record.QUAL

    @property
    def ref(self):
        return self.record.REF

    @property
    def num_called(self):
        return self.record.num_called

    @property
    def call_rate(self):
        return self.record.num_called / len(self.reader.samples)

    @property
    def is_snp(self):
        return self.record.is_snp

    def __str__(self):
        return str(self.record)

    def __unicode__(self):
        return self.record.__unicode__()

    def __repr__(self):
        return repr(self.record)

    @property
    def filters(self):
        return self.record.FILTER

    @filters.setter
    def filters(self, value):
        self.record.FILTER = value

    def add_filter(self, filter_name):
        self.record.add_filter(filter_name)

    @property
    def infos(self):
        return self.record.INFO

    def add_info(self, *args, **kwargs):
        self.record.add_info(*args, **kwargs)

    @property
    def kind(self):
        # return snv type. [snp, indel, unknown]
        return self.record.var_type

    @property
    def is_indel(self):
        return self.record.is_indel

    @property
    def calldata_class(self):
        return make_calldata_tuple(self.record.FORMAT.split(':'))

    def copy_mapping_calls(self, call_mapper):
        calls = []
        sample_indexes = {}
        is_list = True if isinstance(call_mapper, list) else False
        for index, call in enumerate(self.calls):
            sample = call.sample
            if is_list:
                call = call_mapper[index]
            else:
                call = call_mapper(call)
            calls.append(call)
            sample_indexes[sample] = index
        record = self.record
        record = pyvcfRecord(record.CHROM, record.POS, record.ID,
                             record.REF, record.ALT, record.QUAL,
                             record.FILTER, record.INFO, record.FORMAT,
                             sample_indexes, calls)
        snv = SNV(record, self.reader,
                  min_calls_for_pop_stats=self.min_calls_for_pop_stats)
        return snv

    def copy(self):
        return self.copy_mapping_calls(lambda x: x.copy(return_pyvcf_call=True))

    def remove_gt_from_het_calls(self):
        def call_mapper(call):
            if call.is_het:
                call = call.copy_setting_gt(gt=None, return_pyvcf_call=True)
            else:
                call = call.call
            return call
        return self.copy_mapping_calls(call_mapper)

    def remove_gt_from_low_qual_calls(self, min_qual):
        'It returns a new SNV with low qual call set to uncalled'
        def call_mapper(call):
            if min_qual is not None and call.gt_qual < min_qual:
                call = call.copy_setting_gt(gt=None, return_pyvcf_call=True)
            else:
                call = call.call
            return call
        return self.copy_mapping_calls(call_mapper)

    def filter_calls_by_sample(self, samples, reverse=False, keep_info=False):
        calls = [self.get_call(sample).call for sample in samples]
        if reverse:
            all_samples = set([c.sample for c in self.calls])
            chosen_samples = set([c.sample for c in calls])
            unchosen_samples = all_samples.difference(chosen_samples)
            calls = [self.get_call(sample).call for sample in unchosen_samples]

        sample_indexes = {call.sample: idx for idx, call in enumerate(calls)}

        int_alleles = set()
        for call in calls:
            if call.called:
                int_alleles.update([int(allele) for allele in call.gt_alleles])

        record = self.record
        alt_alleles = [a for i, a in enumerate(record.ALT) if i + 1 in int_alleles]
        if not alt_alleles:
            alt_alleles.append(None)
        info = record.INFO if keep_info else {}
        record = pyvcfRecord(record.CHROM, record.POS, record.ID,
                             record.REF, alt_alleles, record.QUAL,
                             record.FILTER, info, record.FORMAT,
                             sample_indexes, calls)
        snv = SNV(record, self.reader,
                  min_calls_for_pop_stats=self.min_calls_for_pop_stats)
        return snv

    @property
    def allele_freqs(self):
        allele_counts = self.allele_counts
        if allele_counts is None:
            return None
        tot_obs = sum(allele_counts.values())
        return {alle: cnt / tot_obs for alle, cnt in allele_counts.viewitems()}

    def is_polymorphic(self, freq_threshold=1):
        allele_freqs = self.allele_freqs
        if not allele_freqs:
            return None
        if len(allele_freqs) == 1:
            return False
        else:
            if freq_threshold == 1:
                return True
            else:
                max_allele_freq = max(allele_freqs.values())
                if max_allele_freq > freq_threshold:
                    return False
                else:
                    return True


class Call(object):
    def __init__(self, call, snv):
        self.call = call
        self.snv = snv
        self._alt_sum_depths = None
        self._allele_depths = {}
        self._has_alternative_counts = None
        self._depths_analyzed = False

    def _get_allele_depths(self):
        if self._depths_analyzed:
            return
        self._depths_analyzed = True
        if not self.call.called:
            return
        snp_caller = self.snv.reader.snpcaller
        if snp_caller is None:
            msg = 'To parse the read depths you have to set the snpcaller'
            raise ValueError(msg)
        elif snp_caller == GATK:
            self._get_depths_gatk()
        elif snp_caller == VARSCAN:
            self._get_depths_varscan()
        elif snp_caller == FREEBAYES:
            self._get_depths_freebayes()
        elif snp_caller == GENERIC:
            self._get_depths_generic()
        else:
            msg = 'SNP caller not supported yet'
            raise NotImplementedError(msg)

    def _get_depths_generic(self):
        self._allele_depths = {}
        self._alt_sum_depths = None
        self._has_alternative_counts = None

    def _get_depths_gatk(self):
        call = self.call
        als = [int(a) for a in call.gt_alleles]
        al_counts = {al_: al_count for al_count, al_ in zip(call.data.AD, als)}
        self._allele_depths = al_counts
        sum_alt = sum(alc for al, alc in al_counts.items() if al != 0)
        self._alt_sum_depths = sum_alt
        self._has_alternative_counts = True

    def _get_depths_varscan(self):
        call = self.call
        data = call.data
        ref_depth = data.RD
        self._allele_depths = {0: ref_depth}
        self._alt_sum_depths = data.AD
        self._has_alternative_counts = False

    def _get_depths_freebayes(self):
        call = self.call
        data = call.data
        ref_depth = data.RO

        alt_depths = data.AO
        if isinstance(alt_depths, int):
            alt_depths = [alt_depths]
        # the number of alternative alleles should be all alleles - 1
        assert len(call.site.alleles) - 1 == len(alt_depths)

        al_dps = {allele + 1: count for allele, count in enumerate(alt_depths)}
        self._alt_sum_depths = sum(al_dps.values())
        al_dps[0] = ref_depth
        self._allele_depths = al_dps
        self._has_alternative_counts = True

    @property
    def depth(self):
        # In freebayes the allele depth (DP) and the sum of allele
        # observations do not match
        if self.snv.reader.snpcaller == FREEBAYES:
            al_dps = self.allele_depths
            depth = sum(al_dps.values())
        else:
            try:
                depth = self.call.data.DP
            except AttributeError:
                depth = None
        return depth

    @property
    def gt_qual(self):
        try:
            return self.call.data.GQ
        except AttributeError:
            return None

    @property
    def ref_depth(self):
        self._get_allele_depths()
        return self._allele_depths.get(0, None)

    @property
    def alt_sum_depths(self):
        self._get_allele_depths()
        return self._alt_sum_depths

    @property
    def allele_depths(self):
        self._get_allele_depths()
        return self._allele_depths

    @property
    def has_alternative_counts(self):
        self._get_allele_depths()
        return self._has_alternative_counts

    @property
    def maf_depth(self):
        al_dps = self.allele_depths
        depth = sum(al_dps.values())
        if depth:
            major_allele_depth = max(al_dps.values())
            return major_allele_depth / depth
        else:
            return None

    @property
    def sample(self):
        return self.call.sample

    @property
    def called(self):
        return self.call.called

    @property
    def gt_type(self):
        return self.call.gt_type

    @property
    def is_het(self):
        return self.call.is_het

    @staticmethod
    def _to_int(allele):
        try:
            return int(allele)
        except ValueError:
            if allele == MISSING_ALLELE_CHAR:
                return None
            raise

    @property
    def int_alleles(self):
        call = self.call
        if call.called:
            return [self._to_int(al) for al in call.gt_alleles]
        else:
            return []

    def __str__(self):
        return str(self.call)

    def __unicode__(self):
        return self.call.__unicode__()

    def __repr__(self):
        return repr(self.call)

    def copy(self, return_pyvcf_call=False):
        return self.copy_calldata({}, return_pyvcf_call)

    def copy_setting_gt(self, gt, return_pyvcf_call=False):
        return self.copy_calldata({'GT': gt}, return_pyvcf_call)

    def copy_calldata(self, dict_with_data, return_pyvcf_call=False):
        snv = self.snv
        calldata_class = snv.calldata_class
        call = self.call
        sampdat = []
        # Access to a protected member. In this case namedtuple _fields
        # is not a protected member
        # pylint: disable=W0212
        for field in calldata_class._fields:
            if field in dict_with_data:
                value = dict_with_data[field]
            else:
                value = getattr(call.data, field)
            sampdat.append(value)

        pyvcf_call = pyvcfCall(self.snv.record, self.call.sample,
                               calldata_class(*sampdat))
        if return_pyvcf_call:
            call = pyvcf_call
        else:
            call = Call(pyvcf_call, snv)
        return call
