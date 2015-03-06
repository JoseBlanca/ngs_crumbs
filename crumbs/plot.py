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

from os.path import splitext

from crumbs.exceptions import OptionalRequirementError
from crumbs.utils.optional_modules import FigureCanvas, Figure, cm, Normalize


FIGURE_SIZE = (15.0, 11.0)  # inche
COLORMAPS = ['Blues', 'Reds', 'Greens', 'Accent', 'Winter', 'Bone', 'Binary']
COLORS = ['blue', 'red', 'green', 'gray']
MARKERGLHYPS = ['^', 's', 'o', '.', 'D', 'h', 's']
MARKER_SIZE2 = 100
MARKER_SIZE3 = 5.0

BAR = 'bar'
LINE = 'line'


def _guess_output_for_matplotlib(fhand):
    'Given an fhand it guesses if we need png or svg'
    output = None
    if fhand is not None:
        output = splitext(fhand.name)[-1].strip('.')
    if not output:
        output = 'png'
    return output


def get_fig_and_canvas(num_rows=1, num_cols=1, figsize=None):
    if figsize is None:
        height = 5.0 * num_rows
        width = 7.5 * num_cols
        if height > 320.0:
            height = 320.0
        figsize = (width, height)
    try:
        fig = Figure(figsize=figsize)
        canvas = FigureCanvas(fig)
    except NameError:
        msg = 'Matplotlib module is required to draw graphical histograms'
        raise OptionalRequirementError(msg)
    return fig, canvas


class HistogramPlotter(object):
    def __init__(self, counters, plots_per_chart=1, kind=LINE, num_cols=1,
                 xlabel=None, ylabel=None, ylog_scale=False, ylimits=None,
                 distrib_labels=None, titles=None, xmax=None, xmin=None,
                 linestyles=None, figsize=None):
        if plots_per_chart > 1 and kind == BAR:
            error_msg = 'if kind is BAR only one plot per chart is allowed'
            raise ValueError(error_msg)
        self.kind = kind
        self.counters = counters
        self.num_cols = num_cols
        self.plots_per_chart = plots_per_chart
        self.ylog_scale = ylog_scale
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.num_plots, self.num_rows = self._get_plot_dimensions()
        fig, canvas = get_fig_and_canvas(num_rows=self.num_rows,
                                         num_cols=self.num_cols,
                                         figsize=figsize)
        self.figure = fig
        self.canvas = canvas
        axes = self._draw_plot(distrib_labels=distrib_labels, ylimits=ylimits,
                               titles=titles, xmax=xmax, xmin=xmin,
                               linestyles=linestyles)
        self.axes = axes

    def _get_plot_dimensions(self):
        num_plots, mod = divmod(len(self.counters), self.plots_per_chart)
        if mod != 0:
            num_plots += 1

        num_rows, mod = divmod(num_plots, self.num_cols)
        if mod != 0:
            num_rows += 1
        return num_plots, num_rows

    def write_figure(self, fhand):
        plot_format = _guess_output_for_matplotlib(fhand)
        self.canvas.print_figure(fhand, format=plot_format)
        fhand.flush()

    def _draw_histogram_in_axe(self, counter, axe, xmax=None, xmin=None,
                               title=None, distrib_label=None, linestyle=None,
                               ylimits=None):

        try:
            distrib = counter.calculate_distribution(max_=xmax, min_=xmin)
        except RuntimeError:
            axe.set_title(title + ' (NO DATA)')
            return axe
        except AttributeError as error:
            # if distributions is None
            err_msg = "'NoneType' object has no attribute "
            err_msg += "'calculate_distribution'"
            if err_msg in error:
                axe.set_title(title + ' (NO DATA)')
                return axe
            raise
        if distrib is None:
            axe.set_title(title + ' (NO DATA)')
            return axe

        counts = distrib['counts']
        bin_limits = distrib['bin_limits']
        if self.ylog_scale:
            axe.set_yscale('log')

        if self.xlabel:
            axe.set_xlabel(self.xlabel)
        if self.ylabel:
            axe.set_ylabel(self.ylabel)
        if title:
            axe.set_title(title)

        if self.kind == BAR:
            xvalues = range(len(counts))
            axe.bar(xvalues, counts)

            # the x axis label
            xticks_pos = [value + 0.5 for value in xvalues]

            left_val = None
            right_val = None
            xticks_labels = []
            for value in bin_limits:
                right_val = value
                if left_val:
                    xticks_label = (left_val + right_val) / 2.0
                    if xticks_label >= 10:
                        fmt = '%d'
                    elif xticks_label >= 0.1 and xticks_label < 10:
                        fmt = '%.1f'
                    elif xticks_label < 0.1:
                        fmt = '%.1e'
                    xticks_label = fmt % xticks_label
                    xticks_labels.append(xticks_label)
                left_val = right_val

            # we don't want to clutter the plot
            num_of_xlabels = 15
            step = int(len(counts) / float(num_of_xlabels))
            step = 1 if step == 0 else step
            xticks_pos = xticks_pos[::step]
            xticks_labels = xticks_labels[::step]
            axe.set_xticks(xticks_pos)
            axe.set_xticklabels(xticks_labels)
        elif self.kind == LINE:
            kwargs = {}
            if distrib_label is not None:
                kwargs['label'] = distrib_label
            if linestyle is not None:
                kwargs['linestyle'] = linestyle

            x_values = []
            for index, i in enumerate(bin_limits):
                try:
                    i2 = bin_limits[index + 1]
                except IndexError:
                    break
                x_values.append((i + i2) / 2.0)
            y_values = counts
            axe.plot(x_values, y_values, **kwargs)

        if ylimits is not None:
            axe.set_ylim(ylimits)

        return axe

    def _draw_plot(self, distrib_labels=None, titles=None, xmax=None,
                   xmin=None, linestyles=None, ylimits=None):
        counter_index = 0
        axes = []
        for plot_num in range(1, self.num_plots + 1):
            # print num_rows, num_cols, plot_num
            axe = self.figure.add_subplot(self.num_rows, self.num_cols,
                                          plot_num)
            for i in range(self.plots_per_chart):
                try:
                    counter = self.counters[counter_index]
                    if distrib_labels is None:
                        distrib_label = None
                    else:
                        distrib_label = distrib_labels[counter_index]
                    if linestyles is None:
                        linestyle = None
                    else:
                        linestyle = linestyles[counter_index]
                except IndexError:
                    break
                title = titles[counter_index] if titles else None
                self._draw_histogram_in_axe(counter, axe=axe, xmin=xmin,
                                            xmax=xmax, title=title,
                                            distrib_label=distrib_label,
                                            linestyle=linestyle,
                                            ylimits=ylimits)
                counter_index += 1

            if distrib_labels is not None:
                axe.legend()
            axes.append(axe)

        return axes


def draw_histogram_in_fhand(counter, fhand, title=None, xlabel=None, xmin=None,
                            xmax=None, ylabel=None, kind=BAR, ylimits=None,
                            ylog_scale=False):
    'It draws an histogram and if the fhand is given it saves it'
    plot_hist = HistogramPlotter([counter], xlabel=xlabel, ylabel=ylabel,
                                 xmax=xmax, xmin=xmin, titles=[title],
                                 ylimits=ylimits)
    plot_hist.write_figure(fhand)


def draw_histograms(counters, fhand, distrib_labels=None, num_cols=2,
                    plots_per_chart=3, xlabel=None, ylabel=None, titles=None,
                    kind=LINE, xmax=None, xmin=None, linestyles=None,
                    ylimits=None, ylog_scale=False):

    plot_hist = HistogramPlotter(counters, xlabel=xlabel, ylabel=ylabel,
                                 xmax=xmax, xmin=xmin, titles=titles,
                                 distrib_labels=distrib_labels, kind=kind,
                                 linestyles=linestyles, ylimits=ylimits,
                                 num_cols=num_cols, ylog_scale=ylog_scale,
                                 plots_per_chart=plots_per_chart)
    plot_hist.write_figure(fhand)


def get_canvas_and_axes(figure_size=FIGURE_SIZE, left=0.1, right=0.9, top=0.9,
                        bottom=0.1, plot_type=111):
    'It returns a matplotlib canvas and axes instance'
    try:
        fig = Figure(figsize=FIGURE_SIZE)
        canvas = FigureCanvas(fig)
    except NameError:
        msg = 'Matplotlib module is required to draw graphical histograms'
        raise OptionalRequirementError(msg)

    axes = fig.add_subplot(plot_type)
    fig.subplots_adjust(left=left, right=right, top=top, bottom=bottom)

    return canvas, axes


def build_histogram(values, fhand, bins=10, range_=None, stacked=False,
                    color=None, label=None, log=False, **kwargs):
    'It draws a histogram of a pandas Series into a file'
    canvas, axes = get_canvas_and_axes()
    plot_format = _guess_output_for_matplotlib(fhand)
    if color is None:
        if label is None:
            axes.hist(values, bins=bins, range=range_, stacked=stacked,
                      log=log)
        else:
            axes.hist(values, bins=bins, range=range_, stacked=stacked,
                      label=label, log=log)
    else:
        axes.hist(values, bins=bins, range=range_, stacked=stacked,
                  label=label, color=color, log=log)
    for key, value in kwargs.items():
        getattr(axes, 'set_{}'.format(key))(value)
    if label is not None:
        axes.legend()

    canvas.print_figure(fhand, format=plot_format)
    fhand.flush()


def draw_scatter(groups, fhand, plot_lines=False, **kwargs):
    # groups is a list of x,y and color_intensity
    canvas, axes = get_canvas_and_axes()
    plot_format = _guess_output_for_matplotlib(fhand)

    for key, value in kwargs.items():
        getattr(axes, 'set_{}'.format(key))(value)

    for index, group in enumerate(groups):
        x_vals = group['x']
        y_vals = group['y']
        sct_kwargs = {}
        sct_kwargs['marker'] = MARKERGLHYPS[index]
        # TODO. Review the API for the colors.
        # The following cases should be posible for color/value:
        #    - a tuple RGB
        #    - nothing. It has to chose one color by default
        #    - a string: 'blue', 'green', etc.
        #    - a list of numbers (intensities) and ColorMap
        # A possible option could be:
        # group['color] for anything that matplotlib could digest for a single
        # color: RGB tuple or string.
        # group['color_intensities'] for the values of the color. In this case
        # a color map should be given or a default one should be used.
        # In this is used 'color' and 'color_intensities' should be
        # incompatible
        # What about plot_lines? That should be incompatible with
        # color_intensities
        color = group.get('color', COLORS[index])
        if 'color_intensity' in group:
            # value is a list with color intensities
            sct_kwargs['c'] = group['color_intensity']
            sct_kwargs['norm'] = Normalize()
            sct_kwargs['cmap'] = cm.get_cmap(COLORMAPS[index])
            sct_kwargs['s'] = 50
            sct_kwargs['alpha'] = 0.5
            sct_kwargs['edgecolor'] = 'white'
        else:
            sct_kwargs['c'] = color

        if plot_lines:
            sct_kwargs['ms'] = MARKER_SIZE3
            sct_kwargs['mec'] = color
            sct_kwargs['mfc'] = color
            axes.plot(x_vals, y_vals, **sct_kwargs)
        else:
            sct_kwargs['s'] = MARKER_SIZE2
            axes.scatter(x_vals, y_vals, **sct_kwargs)

    canvas.print_figure(fhand, format=plot_format)
    fhand.flush()


def lines_to_plot(series, axes):
    for line in series:
        kwargs = {}
        if 'label' in line:
            kwargs['label'] = line['label']
        axes.plot(line['x'], line['y'], **kwargs)
        axes.set_title('title')
    axes.legend()


def draw_lines(series, fhand):
    ''''It draws a line plot with the given series.
    Series is a list of dictionaries, each dictionary at lest needs to contain
    x and y as lists'''
    canvas, axes = get_canvas_and_axes()
    plot_format = _guess_output_for_matplotlib(fhand)
    lines_to_plot(series, axes)

    canvas.print_figure(fhand, format=plot_format)
    fhand.flush()


def draw_density_plot(xs, ys, fhand=None, n_bins=40, canvas=None, axes=None,
                      range_=None, **kwargs):
    plot_format = _guess_output_for_matplotlib(fhand)
    if canvas is None and axes is None:
        canvas, axes = get_canvas_and_axes()
        using_given_axes = False
    elif canvas is not None and axes is not None and fhand is None:
        using_given_axes = True
    else:
        msg = 'If an axes is given the canvas is also required and no fhand'
        raise NotImplementedError(msg)

    for key, value in kwargs.items():
        getattr(axes, 'set_{}'.format(key))(value)

    # TODO check with norm=LogNorm() and/or normed=True
    if range_ is None:
        axes.hist2d(xs, ys, bins=n_bins)
    else:
        axes.hist2d(xs, ys, bins=n_bins, range=range_)

    if not using_given_axes:
        canvas.print_figure(fhand, format=plot_format)
        fhand.flush()


def draw_int_boxplot(boxplot, fhand=None, axes=None, title=None,
                     xlabel=None, ylabel=None):

    if axes is None:
        canvas, axes = get_canvas_and_axes()
        using_given_axes = False
    elif fhand is None:
        using_given_axes = True
    else:
        msg = 'If an axes is not given the fhand also required'
        raise NotImplementedError(msg)

    bar_width = 0.8

    x_vals = sorted(boxplot.counts.keys())

    # for the following line to work we would have to create all points in
    # memory, this is why we have reimplemented the boxplot
    # axes.boxplot([list(boxplot.counts[x_val].elements()) for x_val in x_vals]

    # we don't want more than 50 bars
    max_n_boxes = 40
    n_boxes = len(x_vals) - 1
    xpurge = n_boxes // max_n_boxes + 1

    xticks_lables = []
    xticks = []
    xpos = []
    for index, x_val in enumerate(x_vals):
        xpos.append(index + 0.5)
        if index % xpurge != 0:
            continue

        xticks_lables.append(str(x_val))
        # we add 0.5 to have the x ticks in  the
        # middle of the bar
        xticks.append(index + 0.5)

        intcounter = boxplot.counts[x_val]
        if intcounter.count < 5:
            # at least for values are required to calculate the qualties
            continue
        quart = intcounter.quartiles
        min_ = intcounter.min
        max_ = intcounter.max
        axes.bar(index + 0.5, quart[2] - quart[0], bottom=quart[0],
                 width=bar_width, facecolor='none', align='center')
        # median
        axes.plot([index + 0.5 - bar_width / 2, index + 0.5 + bar_width / 2],
                  [quart[1], quart[1]], color='red')
        # max
        axes.plot([index + 0.5 - bar_width / 4, index + 0.5 + bar_width / 4],
                  [max_, max_], color='black')
        axes.plot([index + 0.5, index + 0.5], [quart[0], min_], '--',
                  color='black')
        axes.plot([index + 0.5 - bar_width / 4, index + 0.5 + bar_width / 4],
                  [min_, min_], color='black')
        axes.plot([index + 0.5, index + 0.5], [quart[2], max_], '--',
                  color='black')

    max_x = len(x_vals) - 1
    axes.set_xlim(0, max_x + 1)

    axes.set_xticks(xticks)
    axes.set_xticklabels(xticks_lables)

    if title:
        axes.set_title(title)
    if xlabel:
        axes.set_xlabel(xlabel)
    if ylabel:
        axes.set_ylabel(ylabel)

    if not using_given_axes:
        canvas.print_figure(fhand)
        fhand.flush()
    return x_vals, xpos, xticks, xticks_lables

