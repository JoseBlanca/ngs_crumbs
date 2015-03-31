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

from os.path import join
import unittest
from tempfile import NamedTemporaryFile
from subprocess import check_output, Popen, PIPE, check_call
from StringIO import StringIO
import os

from crumbs.vcf.snv import VCFReader

from crumbs.vcf.filters import (PASSED, FILTERED_OUT, group_in_filter_packets,
                                CallRateFilter, BiallelicFilter, IsSNPFilter,
                                SnvQualFilter, ObsHetFilter, MafFilter,
                                filter_snvs, MonomorphicFilter,
                                WeirdSegregationFilter,
                                WeirdRecombFilter)
from crumbs.utils.bin_utils import VCF_BIN_DIR
from crumbs.utils.test_utils import TEST_DATA_DIR

# Method could be a function
# pylint: disable=R0201
# Too many public methods
# pylint: disable=R0904
# Missing docstring
# pylint: disable=C0111

VCF_PATH = join(TEST_DATA_DIR, 'sample.vcf.gz')
VCF_INDEL_PATH = join(TEST_DATA_DIR, 'sample_indel.vcf.gz')
VARI_VCF_PATH = join(TEST_DATA_DIR, 'vari_filter.vcf')

VCF_HEADER = '''##fileformat=VCFv4.1
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
'''
VCF = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.
20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4
20\t1230237\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t./.:35:4\t0/2:17:2\t1/1:40:3
'''


def filter_vcf(vcf_fhand, filter_):
    snps = VCFReader(vcf_fhand, min_calls_for_pop_stats=1).parse_snvs()
    packet = list(group_in_filter_packets(snps, 10))[0]
    filtered_packet = filter_(packet)
    return filtered_packet


def eval_prop_in_packet(packet, prop):
    eval_prop = lambda snps: [getattr(snp, prop)for snp in snps]
    packet = {PASSED: eval_prop(packet[PASSED]),
              FILTERED_OUT: eval_prop(packet[FILTERED_OUT])}
    return packet


class FiltersTest(unittest.TestCase):

    def test_group_in_filter_packets(self):
        items = list(range(10))
        packets = list(group_in_filter_packets(items, 4))
        res = [{FILTERED_OUT: [], PASSED: (0, 1, 2, 3)},
               {FILTERED_OUT: [], PASSED: (4, 5, 6, 7)},
               {FILTERED_OUT: [], PASSED: (8, 9)}]
        assert res == packets

    def test_biallelic(self):
        packet = filter_vcf(open(VCF_PATH), filter_=BiallelicFilter())
        res = eval_prop_in_packet(packet, 'alleles')
        assert not(res[FILTERED_OUT])
        assert len(res[PASSED]) == 10

        packet = filter_vcf(open(VCF_INDEL_PATH), filter_=BiallelicFilter())
        res = eval_prop_in_packet(packet, 'alleles')
        assert len(res[FILTERED_OUT]) == 1
        assert len(res[PASSED]) == 6

        # test with only some samples
        kwargs = {'samples_to_consider': ('pepo', 'mu16')}
        packet = filter_vcf(open(VCF_PATH), filter_=BiallelicFilter(**kwargs))
        res = eval_prop_in_packet(packet, 'alleles')
        assert not(res[FILTERED_OUT])
        assert len(res[PASSED]) == 10

        kwargs = {'samples_to_consider': ('pepo', 'mu16')}
        packet = filter_vcf(open(VCF_INDEL_PATH),
                            filter_=BiallelicFilter(**kwargs))

        res = eval_prop_in_packet(packet, 'alleles')
        assert len(res[PASSED]) == 7

    def test_monomorphic(self):
        fhand = NamedTemporaryFile()
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0/0:48:1:51,51\t0/0:48:8:51,51\t0/0:43:5:.,.
20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0/0:49:3:58,50\t0/1:3:5:65,3\t0/0:41:3
20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0/0:49:3:58,50\t0/1:3:5:65,3\t0/0:41:3
        '''
        fhand.write(VCF_HEADER + vcf)
        fhand.flush()
        packet = filter_vcf(open(fhand.name),
                            filter_=MonomorphicFilter())
        assert len(packet[PASSED]) == 2
        assert len(packet[FILTERED_OUT]) == 1

    def test_is_snp(self):
        packet = filter_vcf(open(VCF_PATH), filter_=IsSNPFilter())
        res = eval_prop_in_packet(packet, 'is_snp')
        assert not res[FILTERED_OUT]
        assert res[PASSED] == [True] * 10

        packet = filter_vcf(open(VCF_INDEL_PATH), filter_=IsSNPFilter())
        res = eval_prop_in_packet(packet, 'is_snp')
        assert res[FILTERED_OUT] == [False] * 7
        assert not res[PASSED]


class MafFilterTest(unittest.TestCase):
    def test_maf(self):
        packet = filter_vcf(open(VCF_PATH), filter_=MafFilter(min_maf=0.6))
        res = eval_prop_in_packet(packet, 'maf')
        assert res[FILTERED_OUT] == [0.5, 0.5, 0.5, None, 0.5]
        assert res[PASSED] == [0.75, 0.75, 1.0, 1.0, 1.0]

        packet = filter_vcf(open(VCF_PATH), filter_=MafFilter(max_maf=0.6))
        res = eval_prop_in_packet(packet, 'maf')
        assert res[FILTERED_OUT] == [0.75, None, 0.75, 1.0, 1.0, 1.0]
        assert res[PASSED] == [0.5, 0.5, 0.5, 0.5]

        packet = filter_vcf(open(VCF_PATH),
                            filter_=MafFilter(min_maf=0.6, max_maf=0.8))
        res = eval_prop_in_packet(packet, 'maf')
        assert res[FILTERED_OUT] == [0.5, 0.5, 0.5, None, 0.5, 1.0, 1.0, 1.0]
        assert res[PASSED] == [0.75, 0.75]

        # some samples
        packet = filter_vcf(open(VCF_PATH),
                            filter_=MafFilter(min_maf=0.6, max_maf=0.8,
                                              samples_to_consider=('pepo',
                                                                   'mu16')))
        res = eval_prop_in_packet(packet, 'maf')
        assert res[FILTERED_OUT] == [0.5, 0.5, 0.5, 0.75, None, 0.5, 0.75, 1.0,
                                     1.0, 1.0]
        assert res[PASSED] == []

    def test_maf_bin(self):
        binary = join(VCF_BIN_DIR, 'filter_vcf_by_maf')

        assert 'positional' in check_output([binary, '-h'])

        in_fhand = NamedTemporaryFile()
        in_fhand.write(VCF_HEADER + VCF)
        in_fhand.flush()
        log_fhand = NamedTemporaryFile()
        hist_fhand = NamedTemporaryFile(suffix='.png')
        cmd = [binary, '-m', '0.7', '-c', '2', '-l', log_fhand.name,
               '-t', hist_fhand.name, in_fhand.name]
        process = Popen(cmd, stderr=PIPE, stdout=PIPE)
        stdout = process.communicate()[0]
        assert "passsed: 2" in open(log_fhand.name).read()
        assert 'fileDate' in stdout

        # You can pipe and not set the template
        vcf_fpath = os.path.join(TEST_DATA_DIR, 'scaff000025.vcf.gz')
        cmd1 = ['zcat', vcf_fpath]
        process1 = Popen(cmd1, stdout=PIPE)
        log_fhand = NamedTemporaryFile()
        binary = join(VCF_BIN_DIR, 'filter_vcf_by_maf')
        cmd = [binary, '-m', '0.7', '-c', '2', '-l', log_fhand.name]
        process2 = Popen(cmd, stderr=PIPE, stdout=PIPE, stdin=process1.stdout)
        stdout = process2.communicate()[0]
        assert len(list(VCFReader(StringIO(stdout)).parse_snvs())) == 69

        # You cannot pipe a compressed file
        vcf_fpath = os.path.join(TEST_DATA_DIR, 'scaff000025.vcf.gz')
        cmd1 = ['cat', vcf_fpath]
        process1 = Popen(cmd1, stdout=PIPE)
        binary = join(VCF_BIN_DIR, 'filter_vcf_by_maf')
        cmd = [binary, '-m', '0.7', '-c', '2', '-l', log_fhand.name]
        process2 = Popen(cmd, stderr=PIPE, stdout=PIPE, stdin=process1.stdout)
        stdout, stderr = process2.communicate()
        assert not stdout
        assert 'tell member' in stderr


class CallRateFilterTest(unittest.TestCase):
    def test_missing_genotypes(self):
        packet = filter_vcf(open(VCF_PATH),
                            filter_=CallRateFilter(min_calls=2))
        res = eval_prop_in_packet(packet, 'num_called')
        expected = {FILTERED_OUT: [0, 1, 1, 1, 1], PASSED: [2, 2, 2, 2, 2]}
        assert res == expected

        packet = filter_vcf(open(VCF_PATH),
                            filter_=CallRateFilter(min_calls=1))
        res = eval_prop_in_packet(packet, 'num_called')
        expected = {FILTERED_OUT: [0], PASSED: [2, 2, 2, 2, 1, 2, 1, 1, 1]}
        assert res == expected

        packet = filter_vcf(open(VCF_PATH),
                            filter_=CallRateFilter(min_calls=1,
                                                   reverse=True))
        res = eval_prop_in_packet(packet, 'num_called')
        expected = {PASSED: [0], FILTERED_OUT: [2, 2, 2, 2, 1, 2, 1, 1, 1]}
        assert res == expected

        packet = filter_vcf(open(VCF_PATH),
                            filter_=CallRateFilter(min_calls=3))
        res = eval_prop_in_packet(packet, 'num_called')
        expected = {PASSED: [], FILTERED_OUT: [2, 2, 2, 2, 0, 1, 2, 1, 1, 1]}
        assert res == expected

        # some samples
        packet = filter_vcf(open(VCF_PATH),
                            filter_=CallRateFilter(min_calls=2,
                            samples_to_consider=('pepo', 'mu16')))
        res = eval_prop_in_packet(packet, 'num_called')
        expected = {PASSED: [], FILTERED_OUT: [2, 2, 2, 2, 0, 1, 2, 1, 1, 1]}
        assert res == expected

    def test_call_rate_bin(self):
        binary = join(VCF_BIN_DIR, 'filter_vcf_by_missing')

        assert 'positional' in check_output([binary, '-h'])

        in_fhand = NamedTemporaryFile()
        in_fhand.write(VCF_HEADER + VCF)
        in_fhand.flush()
        plot_fhand = NamedTemporaryFile(suffix='.png')
        cmd = [binary, '-m', '3', '-t', plot_fhand.name, in_fhand.name]
        process = Popen(cmd, stderr=PIPE, stdout=PIPE)
        stdout, stderr = process.communicate()
        assert "passsed: 5" in stderr
        assert 'fileDate' in stdout


class ObsHetFilterTest(unittest.TestCase):
    def test_obs_het(self):
        packet = filter_vcf(open(VCF_PATH),
                            filter_=ObsHetFilter(min_het=0.5))
        res = eval_prop_in_packet(packet, 'obs_het')
        assert res[FILTERED_OUT] == [0.0, 0.0, 0.0, None, 0.0, 0.0, 0.0]
        assert res[PASSED] == [0.5, 1.0, 0.5]

        packet = filter_vcf(open(VCF_PATH),
                            filter_=ObsHetFilter(max_het=0.5))
        res = eval_prop_in_packet(packet, 'obs_het')
        assert res[FILTERED_OUT] == [None, 1.0]
        assert res[PASSED] == [0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0]

        packet = filter_vcf(open(VCF_PATH), filter_=ObsHetFilter(min_het=0.1,
                                                                 max_het=0.9))
        res = eval_prop_in_packet(packet, 'obs_het')
        assert res[FILTERED_OUT] == [0.0, 0.0, 0.0, None, 1.0, 0.0, 0.0, 0.0]
        assert res[PASSED] == [0.5, 0.5]

        # some samples
        packet = filter_vcf(open(VCF_PATH),
                            filter_=ObsHetFilter(min_het=0.1, max_het=0.9,
                                                 samples_to_consider=('pepo',
                                                                      'mu16')))
        res = eval_prop_in_packet(packet, 'obs_het')
        assert res[FILTERED_OUT] == [0.0, 0.0, 0.0, 0.5, None, 1.0, 0.5, 0.0,
                                     0.0, 0.0]
        assert res[PASSED] == []

    def test_obs_het_bin(self):
        binary = join(VCF_BIN_DIR, 'filter_vcf_by_het')

        assert 'positional' in check_output([binary, '-h'])

        in_fhand = NamedTemporaryFile()
        in_fhand.write(VCF_HEADER + VCF)
        in_fhand.flush()
        plot_fhand = NamedTemporaryFile(suffix='.png')
        cmd = [binary, '-m', '0.5', '-c', '2', '-t', plot_fhand.name,
               in_fhand.name]
        process = Popen(cmd, stderr=PIPE, stdout=PIPE)
        stdout, stderr = process.communicate()
        assert "passsed: 3" in stderr
        assert 'fileDate' in stdout


class SnvQualTest(unittest.TestCase):
    def test_snv_qual(self):
        packet = filter_vcf(open(VCF_PATH), filter_=SnvQualFilter(20))
        res = eval_prop_in_packet(packet, 'qual')
        assert res[FILTERED_OUT] == [None] * 10
        assert not res[PASSED]

        fpath = join(TEST_DATA_DIR, 'freebayes_al_depth.vcf')
        filter_ = SnvQualFilter(1)
        packet = filter_vcf(open(fpath), filter_=filter_)
        res = eval_prop_in_packet(packet, 'qual')
        assert len(res[FILTERED_OUT]) == 1
        assert len(res[PASSED]) == 4

        plot_fhand = NamedTemporaryFile(suffix='.png')
        filter_.plot_hist(plot_fhand)


def get_snv_pos(vcf_fhand):
    pos = []
    for line in vcf_fhand:
        if line.startswith('#'):
            continue
        pos.append(int(line.split()[1]))
    return pos


class FilterVcfFunctTest(unittest.TestCase):
    def test_filter_fhand(self):
        in_fhand = StringIO(VCF_HEADER + VCF)
        out_fhand = StringIO()
        filter_snvs(in_fhand, out_fhand, filters=[])
        res = get_snv_pos(StringIO(out_fhand.getvalue()))
        in_pos = get_snv_pos(StringIO(in_fhand.getvalue()))
        assert in_pos == res

        in_fhand = StringIO(VCF_HEADER + VCF)
        out_fhand = StringIO()
        filtered_fhand = StringIO()
        log_fhand = StringIO()
        filter_snvs(in_fhand, out_fhand, filters=[BiallelicFilter()],
                    filtered_fhand=filtered_fhand, log_fhand=log_fhand)
        res = get_snv_pos(StringIO(out_fhand.getvalue()))
        filtered = get_snv_pos(StringIO(filtered_fhand.getvalue()))
        in_pos = get_snv_pos(StringIO(in_fhand.getvalue()))
        assert res == [14370, 17330, 1230237]
        assert filtered == [1110696, 1234567, 1234567]
        assert 'SNVs passsed: 3' in log_fhand.getvalue()


class BinaryFilterTest(unittest.TestCase):
    def test_biallelic_binary(self):
        binary = join(VCF_BIN_DIR, 'filter_vcf_by_biallelic')

        assert 'positional' in check_output([binary, '-h'])

        in_fhand = NamedTemporaryFile()
        in_fhand.write(VCF_HEADER + VCF)
        in_fhand.flush()
        out_fhand = NamedTemporaryFile()
        filtered_fhand = NamedTemporaryFile()
        cmd = [binary, '-o', out_fhand.name, '-f', filtered_fhand.name,
               in_fhand.name]
        process = Popen(cmd, stderr=PIPE, stdout=PIPE)
        stderr = process.communicate()[-1]
        assert "passsed: 3" in stderr

        res = get_snv_pos(open(out_fhand.name))
        filtered = get_snv_pos(open(filtered_fhand.name))
        assert res == [14370, 17330, 1230237]
        assert filtered == [1110696, 1234567, 1234567]

        # with stdout
        cmd = [binary, '-f', filtered_fhand.name, in_fhand.name]
        process = Popen(cmd, stderr=PIPE, stdout=PIPE)
        stdout, stderr = process.communicate()
        assert "passsed: 3" in stderr
        res = get_snv_pos(StringIO(stdout))
        assert res == [14370, 17330, 1230237]

        # with stdin
        in_fhand = NamedTemporaryFile()
        in_fhand.write(VCF_HEADER + VCF)
        in_fhand.flush()
        stdin = open(in_fhand.name)
        cmd = [binary, '-f', filtered_fhand.name]
        process = Popen(cmd, stderr=PIPE, stdout=PIPE, stdin=stdin)
        stdout, stderr = process.communicate()
        assert "passsed: 3" in stderr
        in_fhand.close()

    def test_by_sample_bin(self):
        binary = join(VCF_BIN_DIR, 'filter_vcf_by_sample')

        assert 'positional' in check_output([binary, '-h'])

        in_fhand = NamedTemporaryFile()
        in_fhand.write(VCF_HEADER + VCF)
        in_fhand.flush()
        cmd = [binary, '-s', 'NA00005', in_fhand.name]
        process = Popen(cmd, stderr=PIPE, stdout=PIPE)
        stderr = process.communicate()[-1]
        assert 'NA00005' in stderr

        in_fhand = NamedTemporaryFile()
        in_fhand.write(VCF_HEADER + VCF)
        in_fhand.flush()
        cmd = [binary, '-s', 'NA00002', in_fhand.name]
        process = Popen(cmd, stderr=PIPE, stdout=PIPE)
        stdout = process.communicate()[0]
        assert 'FORMAT\tNA00002\n' in stdout


def _create_vcf_file(vcf_string):
    vcf_fhand = NamedTemporaryFile(suffix='.vcf', delete=False)
    vcf_fhand.write(vcf_string)
    vcf_fhand.flush()
    compress_cmd = ['bgzip', vcf_fhand.name]
    check_call(compress_cmd)
    vcf_fpath = vcf_fhand.name + '.gz'
    tabix_cmd = ['tabix', '-p', 'vcf', vcf_fpath]
    tabix_index_fpath = vcf_fpath + '.tbi'
    check_call(tabix_cmd)
    return vcf_fpath, tabix_index_fpath


class ConsistentSegregationTest(unittest.TestCase):

    def test_consistent_segregation(self):
        vcf_fpath = os.path.join(TEST_DATA_DIR, 'scaff000025.vcf.gz')
        snv_filter = WeirdSegregationFilter(min_num_snvs_check_in_win=2,
                                            num_snvs_check=200)
        flt_snps = snv_filter.filter_vcf(vcf_fpath)
        num_flt_snps = len(list(flt_snps))
        assert num_flt_snps == 281
        plot_fhand = NamedTemporaryFile(suffix='.png')
        snv_filter.plot_failed_freq_dist(plot_fhand)

        fhand = StringIO()
        snv_filter.write_log(fhand)
        assert 'SNVs passsed: 281' in fhand.getvalue()

#         plot_dir = TemporaryDir()
#         snv_filter = WeirdSegregationFilter(min_num_snvs_check_in_win=2,
#                                             num_snvs_check=200,
#                                             debug_plot_dir=plot_dir.name)
#         flt_snps = snv_filter.filter_vcf(vcf_fpath)
#         list(flt_snps)
#         plot_dir.close()

    def test_bin(self):
        binary = join(VCF_BIN_DIR, 'filter_vcf_by_weird_segregation')
        cmd = [binary, '-h']
        process = Popen(cmd, stderr=PIPE, stdout=PIPE)
        stdout = process.communicate()[0]
        assert 'usage' in stdout

        vcf_fpath = os.path.join(TEST_DATA_DIR, 'scaff000025.vcf.gz')
        binary = join(VCF_BIN_DIR, 'filter_vcf_by_weird_segregation')
        cmd = [binary, '-n', '2', '-m', '200', '-s', '1_14_1_gbs', '-s',
               '1_17_1_gbs', '-s', '1_18_4_gbs', '-s', '1_19_4_gbs', '-s',
               '1_26_1_gbs', '-s', '1_27_1_gbs', '-s', '1_2_2_gbs', '-s',
               '1_35_13_gbs', '-s', '1_3_2_gbs', '-s', '1_50_1_gbs', '-s',
               '1_59_1_gbs', '-s', '1_63_4_gbs', '-s', '1_6_2_gbs', '-s',
               '1_70_1_gbs', '-s', '1_74_1_gbs', '-s', '1_79_1_gbs', '-s',
               '1_7_2_gbs', '-s', '1_81_10_gbs', '-s', '1_86_1_gbs', '-s',
               '1_8_2_gbs', '-s', '1_91_2_gbs', '-s', '1_94_4_gbs', '-s',
               '2_107_1_gbs', '-s', '2_10_2_gbs', '-s', '2_116_1_gbs', '-s',
               '2_11_1_gbs', '-s', '2_125_2_gbs', '-s', '2_13_1_gbs',
               vcf_fpath]
        process2 = Popen(cmd, stderr=PIPE, stdout=PIPE)
        stdout, stderr = process2.communicate()
        assert len(list(VCFReader(StringIO(stdout)).parse_snvs())) == 273
        assert 'SNVs processed:' in stderr


class ConsistentRecombinationTest(unittest.TestCase):

    def test_cons_recomb(self):
        vcf_fpath = os.path.join(TEST_DATA_DIR, 'scaff000025.vcf.gz')
        snvs = VCFReader(open(vcf_fpath)).parse_snvs()
        snv_filter = WeirdRecombFilter(pop_type='ril_self')
        flt_snvs = snv_filter.filter_snvs(snvs)
        assert len(list(flt_snvs)) == 258
        assert snv_filter.not_fitted_counter['no close region left'] == 10
        fhand = NamedTemporaryFile(suffix='.png')
        flt_snvs = snv_filter.plot_recomb_at_0_dist_hist(fhand)
        assert len(snv_filter.recomb_rates['ok']) == 245
        assert len(snv_filter.recomb_rates['ok_conf_is_None']) == 13
        assert len(snv_filter.recomb_rates['not_ok']) == 14

        snvs = VCFReader(open(vcf_fpath)).parse_snvs()
        snv_filter = WeirdRecombFilter(pop_type='ril_self',
                                       max_zero_dist_recomb=0.07,
                                       alpha_recomb_0=None)
        flt_snvs = snv_filter.filter_snvs(snvs)
        assert len(list(flt_snvs)) == 266
        assert snv_filter.not_fitted_counter['no close region left'] == 10
        fhand = NamedTemporaryFile(suffix='.png')
        flt_snvs = snv_filter.plot_recomb_at_0_dist_hist(fhand)
        assert len(snv_filter.recomb_rates['ok']) == 0
        assert len(snv_filter.recomb_rates['ok_conf_is_None']) == 266
        assert len(snv_filter.recomb_rates['not_ok']) == 6

        fhand = StringIO()
        snv_filter.write_log(fhand)
        assert 'SNVs processed: 282' in fhand.getvalue()

    def test_bin(self):
        binary = join(VCF_BIN_DIR, 'filter_vcf_by_weird_recomb')
        cmd = [binary, '-h']
        process = Popen(cmd, stderr=PIPE, stdout=PIPE)
        stdout = process.communicate()[0]
        assert 'usage' in stdout

        # You cannot pipe a compressed file
        vcf_fpath = os.path.join(TEST_DATA_DIR, 'scaff000025.vcf.gz')
        cmd1 = ['cat', vcf_fpath]
        process1 = Popen(cmd1, stdout=PIPE)
        binary = join(VCF_BIN_DIR, 'filter_vcf_by_weird_recomb')
        cmd = [binary, '--pop_type', 'ril_self', '--window', '60']
        process2 = Popen(cmd, stderr=PIPE, stdout=PIPE, stdin=process1.stdout)
        stdout, stderr = process2.communicate()
        assert not stdout
        assert 'RuntimeError' in stderr

        # You can pipe and not set the template
        vcf_fpath = os.path.join(TEST_DATA_DIR, 'scaff000025.vcf.gz')
        cmd1 = ['zcat', vcf_fpath]
        process1 = Popen(cmd1, stdout=PIPE)
        binary = join(VCF_BIN_DIR, 'filter_vcf_by_weird_recomb')
        cmd = [binary, '--pop_type', 'ril_self', '--window', '60']
        process2 = Popen(cmd, stderr=PIPE, stdout=PIPE, stdin=process1.stdout)
        stdout, stderr = process2.communicate()
        assert len(list(VCFReader(StringIO(stdout)).parse_snvs())) == 252
        assert 'SNVs processed:' in stderr

        # with a standard file
        vcf_fpath = os.path.join(TEST_DATA_DIR, 'scaff000025.vcf.gz')
        binary = join(VCF_BIN_DIR, 'filter_vcf_by_weird_recomb')
        cmd = [binary, '--pop_type', 'ril_self', vcf_fpath]
        process = Popen(cmd, stderr=PIPE, stdout=PIPE)
        stdout, stderr = process.communicate()
        assert len(list(VCFReader(StringIO(stdout)).parse_snvs())) == 258


if __name__ == "__main__":
    # import sys; sys.argv = ['', 'MafFilterTest']
    unittest.main()
