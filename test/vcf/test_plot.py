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
import itertools
import random
import unittest
from collections import OrderedDict

from crumbs.utils.file_utils import TemporaryDir
from crumbs.vcf.plot import plot_in_genome, plot_haplotypes
from crumbs.utils.test_utils import TEST_DATA_DIR
from tempfile import NamedTemporaryFile


class GenomePlot(unittest.TestCase):
    def gen_windows(self, chrom='ch1'):
        step = 10
        for x in range(0, 1000, step):
            start = x
            end = x + step
            value = random.random()
            value2 = random.random()
            yield {'start': start, 'end': end, 'chrom': chrom,
                   'values': {'val1': value, 'val2': value2}}

    def test_plot_window(self):
        iterator = itertools.chain(self.gen_windows(),
                                   self.gen_windows('ch2'))
        tempdir = TemporaryDir()
        out_base = join(tempdir.name, 'out')
        labels = OrderedDict({'val1': {'title': 'val1 title',
                                       'ylabel': 'val1 ylabel'},
                              'val2': {'title': 'val2 title',
                                       'ylabel': 'val2 ylabel'}})

        plot_in_genome(iterator, out_base=out_base, labels=labels)
        # raw_input(tempdir.name)
        tempdir.close()


class HaplotypePlot(unittest.TestCase):
    def test_plot_haplo(self):
        vcf_fhand = open(join(TEST_DATA_DIR, 'freebayes_multisample.vcf.gz'))
        plot_fhand = NamedTemporaryFile(suffix='.png')
        plot_haplotypes(vcf_fhand, plot_fhand)
        # raw_input(plot_fhand.name)

if __name__ == "__main__":
    #import sys; sys.argv = ['', 'HaplotypePlot']
    unittest.main()
