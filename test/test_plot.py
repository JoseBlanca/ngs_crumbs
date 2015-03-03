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

import unittest
from tempfile import NamedTemporaryFile
from array import array

import numpy

from crumbs.plot import (build_histogram, draw_scatter, draw_density_plot,
                         draw_histogram_in_fhand, draw_histograms,
                         draw_int_boxplot, HistogramPlotter)
from crumbs.statistics import IntCounter, IntBoxplot


class PlotTests(unittest.TestCase):
    def test_build_histogram(self):
        values = [1, 2, 3, 1, 2, 3, 2, 3, 2, 3, 2, 1, 4]
        fhand = NamedTemporaryFile(suffix='.png')
        build_histogram(values, fhand, bins=100, title='hola')
        # raw_input(fhand.name)

    def test_draw_scatter(self):
        fhand = NamedTemporaryFile(suffix='.png')
        data = [{'y': array('B', [40, 60, 13]),
                 'x': array('B', [20, 45, 45]),
                 'value': array('B', [9, 90, 90])}]

        draw_scatter(data, fhand, xlim=0, ylim=0)
        # raw_input(fhand.name)

    def test_draw_scatter_line(self):
        fhand = NamedTemporaryFile(suffix='.png')
        data = [{'y': array('B', [40, 60, 13]),
                 'x': array('B', [20, 45, 45])}]
        draw_scatter(data, fhand, plot_lines=True)
        # raw_input(fhand.name)

    def test_draw_density_plot(self):
        fhand = NamedTemporaryFile(suffix='.png')
        data = {'y': array('f', numpy.random.normal(40, 5, 40000)),
                'x': array('f', numpy.random.normal(40, 5, 40000))}
        draw_density_plot(data['x'], data['y'], fhand, n_bins=200)

    def test_int_boxplot(self):
        box = IntBoxplot()
        box.append(1, 50)
        box.append(1, 40)
        box.append(1, 30)
        box.append(1, 40)
        box.append(2, 30)
        box.append(2, 10)
        box.append(2, 20)
        box.append(2, 40)

        fhand = NamedTemporaryFile(suffix='.png')
        draw_int_boxplot(box, fhand=fhand)
        # raw_input(fhand.name)


class HistrogramPlotterTest(unittest.TestCase):
    def test_histogram_plotter(self):
        values = [1, 2, 3, 1, 2, 3, 2, 3, 2, 3, 2, 1, 4]
        counter = IntCounter(values)
        counters = [counter]
        histo_ploter = HistogramPlotter(counters)
        assert len(histo_ploter.axes) == len(counters)
        fhand = NamedTemporaryFile(suffix='.png')
        histo_ploter.write_figure(fhand)
        # raw_input(fhand.name)

        # Add more intcounters
        counters.append(IntCounter(values))
        counters.append(IntCounter(values))
        histo_ploter = HistogramPlotter(counters,
                                        distrib_labels=['1', '2', '3'])
        assert len(histo_ploter.axes) == len(counters)
        # raw_input(fhand.name)

    def test_draw_histogram_in_fhand(self):
        values = [1, 2, 3, 1, 2, 3, 2, 3, 2, 3, 2, 1, 4]
        fhand = NamedTemporaryFile(suffix='.png')
        counter = IntCounter(values)
        draw_histogram_in_fhand(counter, fhand=fhand)
        # raw_input(fhand.name)

    def tests_draw_histograms(self):
        fhand = NamedTemporaryFile(suffix='.png')
        values = [1, 2, 3, 1, 2, 3, 2, 3, 2, 3, 2, 1, 4]
        counters = []
        counters.append(IntCounter(values))
        counters.append(IntCounter(values))
        counters.append(IntCounter(values))
        titles = ['t1', 't2', 't3']
        draw_histograms(counters, fhand, titles=titles, plots_per_chart=2)
        # raw_input(fhand.name)

        values = [1, 2, 3, 1, 2, 3, 2, 3, 2, 3, 2, 1, 4, 0, 5, 4, 4, 4, 4, 4]
        counters = []
        counters.append(IntCounter(values))
        counters.append(IntCounter(values))
        counters.append(IntCounter(values))
        titles = ['t1', 't2', 't3']
        draw_histograms(counters, fhand, titles=titles, plots_per_chart=2,
                        ylimits=(0, 14))
        # raw_input(fhand.name)

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'HistrogramPlotterTest.tests_draw_histograms']
    unittest.main()
