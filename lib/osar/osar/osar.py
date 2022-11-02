#!/usr/bin/python
# -*-coding: utf-8 -*-
# Author: Joses Ho
# Email : joseshowh@gmail.com

from .munger import munger
from .misc import misc
from .track_plots import track_plots


def load(filename):
    '''Loads a saved OSAR object.'''
    import pickle as pick

    with open(filename, 'rb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        out = pick.load(f)

    return out

class osar(object):

    def __init__(self,
                 path_to_data,
                 driver,
                 countlog_folder='countlog',
                 border_folder='borders',
                 framerate_fps=25,
                 rolling_average_window_seconds=1,
                 traversal_time_window_seconds=1,
                 choicezone_width_mm=6.,
                 light_on_time_start=30,
                 light_off_time_start=60,
                 border_shift_mm=0
                 ):
        '''
        Create an `osar` object.

        Supply a `countlog_folder` with raw countlogs and a
        `border_folder` (default 'borders'). You must also enter the driver used
        in the experiment. All the countlogs in `countlog_folder` MUST share the
        same driver.

        Directory structure:
        path_to_data/
        |-- countlog_folder/
        |   |-- CountLog_date_genotype_light-intensity_light-color-1.csv
        |   |-- CountLog_date_genotype_light-intensity_light-color-2.csv
        |   .   any
        |   .   number
        |   .   of
        |   .   countlogs...
        |   |-- CountLog_date_genotype_light-intensity_light-color-n.csv
        |
        |-- border_folder/
        |   |-- Pattern 01.csv
        |   |-- Pattern 10.csv


        Keywords
        --------

        path_to_data: string
            Path to a folder with at least one countlog. All countlogs must share
            the same driver.

        driver: string
            The genetic driver used in the experiment.

        countlog_folder: string, default 'countlog'
            The name of the folder containing the countlogs.

        border_folder: string, default 'borders'
            The name of the folder containing the border files from CRITTA.

        framerate_fps: int, default 25
            The (estimated) framerate of the raw video data.

        rolling_average_window_seconds: integer, default 1
            The window size in seconds used to smooth the countlog data.

        traversal_time_window_seconds: integer, default 1
            The speeds before and after traversals will be computed for this
            time window. The default is 1 second before and after the
            traversal.

        choicezone_width_mm: float, default 6
            The total width of the choicezone. Therefore, the choicezone extends
            by half the value on _each side_ of the light-dark border.

        light_on_time_start: integer, default 30
            The start of the time window used to compute time preference metrics.
            This is in terms of the illumination epoch in seconds, where t = 0 is
            the start of the illumination.

        light_off_time_start: integer, default 60
            The end of the time window used to compute time preference metrics.
            This is in terms of the illumination epoch in seconds, where t = 0 is
            the start of the illumination.
            
        border_shift_mm: float, default 0
            How much to extend the lighted zones outward by.
        '''
        import warnings
        warnings.filterwarnings("ignore", message="numpy.dtype size changed")

        import numpy as np
        import pandas as pd

        from .munger import munger

        # Can comment out during rapid dev.
        from .contrasts import contrast_plotter
        from .timeseries import timeseries_plotter

        # Now, include with each created OSAR object, the version of OSAR used
        # to create it.
        self.osar_version = '0.23.5'

        # pixel length in mm. You shouldn't have to change this.
        self.PIXEL_LENGTH_MM = 0.104
        # Width of the entire choice zone in pix,
        # on both sides of the light-dark border.
        # So if ZONE_WIDTH_mm=10, there is a 5mm width
        # on each side of the light-dark border.
        # decimal point at end of number to make it a float.
        self.CHOICEZONE_WIDTH_MM = np.float(choicezone_width_mm)
        self.ZONE_WIDTH = (self.CHOICEZONE_WIDTH_MM/2) / self.PIXEL_LENGTH_MM
        ## Rolling average window size
        ## How many seconds to smooth.
        self.FRAME_RATE = framerate_fps
        self.ROLLING_AVERAGE_WIN_SIZE = np.float(rolling_average_window_seconds)
        self.WINDOW_SIZE = int(self.ROLLING_AVERAGE_WIN_SIZE * self.FRAME_RATE)
        # Before/after traversal window size.
        # This is the time window in seconds that we will extract before and after traversals
        # to compute speed before and after the crossing.
        self.TRAVERSAL_WINDOW_SIZE_SECONDS=np.float(traversal_time_window_seconds)
        # get window in terms of index size.
        self.TRAVERSAL_WINDOW_IDX=int(self.TRAVERSAL_WINDOW_SIZE_SECONDS * self.FRAME_RATE)
        # Window selection used for PI calculations.
        # This is in terms of the illumination epoch in seconds.
        # t = 0 is the start of the illumination.
        self.LIGHT_ON_TIME_START = light_on_time_start
        self.LIGHT_OFF_TIME_START = light_off_time_start
        self.driver = driver
        self.border_shift_mm = border_shift_mm


        print('Creating borders for each fly...')
        self.borders = munger.create_borders(path_to_data,border_folder)
        
        # Adjust for border_shift
        bshift_pixels = self.border_shift_mm / self.PIXEL_LENGTH_MM
        
        idxx = pd.IndexSlice
        
        p01_mz1and3 = self.borders.loc[idxx[:,['Pattern 01']], 
                                       idxx["midzone1", "midzone3"]]
        p01_mz2 = self.borders.loc[idxx[:,['Pattern 01']], 
                                       idxx["midzone2"]]

        self.borders.loc[idxx[:,['Pattern 01']], 
                         idxx["midzone1", "midzone3"]] = p01_mz1and3 - bshift_pixels           
        self.borders.loc[idxx[:,['Pattern 01']], 
                         idxx["midzone2"]] = p01_mz2 + bshift_pixels      



        p10_mz1and3 = self.borders.loc[idxx[:,['Pattern 10']], 
                                       idxx["midzone1", "midzone3"]]
        p10_mz2 = self.borders.loc[idxx[:,['Pattern 10']], 
                                       idxx["midzone2"]]

        self.borders.loc[idxx[:,['Pattern 10']], 
                         idxx["midzone1", "midzone3"]] = p10_mz1and3 + bshift_pixels           
        self.borders.loc[idxx[:,['Pattern 10']], 
                          idxx["midzone2"]] = p10_mz2 - bshift_pixels 
        
        
        print('Done.\n')
        self.flies, self.tracks, \
                self.results, self.weird_flies = munger.process_countlogs(
                                                  self.borders,
                                                  driver,
                                                  path_to_data,
                                                  countlog_folder,
                                                  self.TRAVERSAL_WINDOW_IDX,
                                                  self.LIGHT_ON_TIME_START,
                                                  self.LIGHT_OFF_TIME_START,
                                                  self.PIXEL_LENGTH_MM,
                                                  self.CHOICEZONE_WIDTH_MM,
                                                  self.WINDOW_SIZE,
                                                  self.ZONE_WIDTH)

        self.metrics = self.results.columns.tolist()[8:]

        # Comment these two lines below during development of either plot.
        # Passes an instance of `self` to plotter.
        self.contrasts = contrast_plotter.contrast_plotter(self)

        # Passes an instance of `self` to timeseries.
        self.timeseries = timeseries_plotter.timeseries_plotter(self)

    # Use the code below for rapid debugging without having to
    # re-create the entire OSAR object.
    def _contrast_dev(self):
        '''PROTOTYPING for Contrast plotting method for OSAR.'''
        from .contrasts import contrast_plotter
        return contrast_plotter.contrast_plotter(self)


    def _timeseries_dev(self):
        '''PROTOTYPING for timeseries plotting method (PI and speed over time)
        for OSAR.'''
        from .timeseries import timeseries_plotter
        # Passes an instance of `self` to timeseries.
        return timeseries_plotter.timeseries_plotter(self)


    def __repr__(self):
        rep_strs = ["Total number of flies: {}".format(len(self.results)),
                    "Light intensities used: {}".format(str(self.results.light_intensity.unique())),
                    "Choice zone width: {}mm".format(self.CHOICEZONE_WIDTH_MM),
                    "Choice time window: {}s".format(self.TRAVERSAL_WINDOW_SIZE_SECONDS),
                    "Smoothing window size: {}s".format(self.ROLLING_AVERAGE_WIN_SIZE),
                    "Created with OSAR v{}".format(self.osar_version)
                    ]
        return '\n'.join(rep_strs)



    def save(self, filename):
        '''Saves the current OSAR object as a Python pickle.'''
        import pickle as pick

        with open(filename, 'wb') as f:
            # To ensure compatibility with Py2, set protocol=2
            pick.dump(self, f, protocol=2)





    def flycounts(self):
        '''Returns the fly counts according to the various categories.'''
        import pandas as pd
        x = self.flies.groupby(['driver','opsin','status',
                        'light_intensity']).count().dropna().to_records()
        x = pd.DataFrame(x)
        # Rename one of the columns (any, really)..
        x.rename(columns={'original_count_log':'fly_count'}, inplace=True)
        out = x[['driver','opsin','status',
                 'light_intensity','fly_count']]
        return out


    def max_velocity_mm_per_s(self):
        '''Returns the maximum velocity observed.'''
        import numpy as np
        v = []
        for f in self.tracks.keys():
            track = self.tracks[f]
            if 'velocity' in track.columns:
                v.append(track.velocity.abs().max() * self.PIXEL_LENGTH_MM)
            else:
                v.append(track.cX_velocity.abs().max() * self.PIXEL_LENGTH_MM)
        return np.max(v)


    def heatmap(self, colormap = 'magma'):
        '''
        For each Light Intensity, 4 heatmaps will be created. There will be
        a speed and occupancy heatmap for sibling controls, and the same two
        heatmaps for offspring experimental flies.

        Keywords
        --------
        colormap: matplotlib colormap
            Choose between 'viridis', 'plasma', 'inferno', and 'magma'.
            See https://matplotlib.org/examples/color/colormaps_reference.html

        Returns a list of matplotlib Figures.
        '''
        from .heatmap_plotter import heatmap_plotter
        out = heatmap_plotter.generic_heatmap_plotter(self)

        return out


    def baseline_comparison(self, pattern='Pattern 01'):
        '''
        Returns a 3x1 array of AxesSubplot objects for each of the metrics:
            'pi_smoothed'
            'light_attraction_index'
            'log2_speed_ratio'
        '''
        import matplotlib.pyplot as plt
        import seaborn as sns

        res = self.results.copy()
        res['groups'] = res.status.str.cat(res.light_intensity, sep='; ')

        sns.set(style='ticks', context='poster')
        out = []
        for _metric in ['pi_smoothed',
                        'light_attraction_index',
                        'log2_speed_ratio']:

            f, axx = plt.subplots(nrows=3, figsize=(10, 12),
                                  sharex=True, sharey=True,
                                  gridspec_kw={'hspace':0.5})

            kw = dict(data=res, x='groups', hue='genotype')
            plot_y = ['{}_BASELINE'.format(_metric),
                      '{}_{}'.format(_metric, pattern),
                      '{}_{}_baseline_corrected'.format(_metric, pattern)]

            for k, a in enumerate(axx):
                sw = sns.swarmplot(ax=a, y=plot_y[k], **kw)
                sw.legend(loc='upper left', bbox_to_anchor=(1,1))
                for item in sw.get_xticklabels():
                    item.set_rotation(45)
                    item.set_horizontalalignment('right')
                a.set_title(plot_y[k], fontsize=15)
                a.set_ylabel('')
                a.set_xlabel('')
                # a.set_ylim(-1.1, 1.1)
                sns.despine(ax=a, trim=True)

            f.suptitle(self.driver, fontsize=20, y=0.92)
            out.append(axx)

        return out[0], out[1], out[2]



    def correlation_plot(self, x, y,
                         group_by,
                         include_intensities=None,
                         drop_nan_from=None,
                         decimal_places=3):
        '''
        Produces correlation plots for any 2 metrics produced by OSAR.

        Keywords
        --------
        x, y: str
            Metrics from OSAR to be plotted.

        group_by: str, "status",  "genotype", or "light_intensity"
            For each category in `group_by`, a correlation plot will be
            produced. Currently only accepts "status" or "genotype.

        include_intensities: list, default None
            Select intensities to include. If None, all intensities are plotted.

        drop_nan_from: boolean, list, default None
            List of columns (in addition to x and y) from which NaNs will be
            dropped.

        decimal_places: int, default 2
            Round off R2 and its 95% CI to the indicated number of decimal
            places.

        Returns a matplotlib Figure.
        '''
        from .plot_helpers import plot_helpers as pth
        import warnings
        import scipy as sp
        import numpy as np
        import pandas as pd
        import scikits.bootstrap as skb

        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle
        from matplotlib.lines import Line2D
        from matplotlib.text import Text

        import seaborn as sns

        allowed_group_bys = ['status', 'genotype', 'light_intensity']
        if group_by not in allowed_group_bys:
            raise ValueError('{} is not a recognized group_by. Please choose one of {}'\
                             .format(group_by, allowed_group_bys))

        results = self.results.copy()

        if include_intensities is not None:
            if isinstance(include_intensities, list):
                # check given intensities are actually in the results.
                not_in_results = [a for a in include_intensities
                                  if a not in
                                  results.light_intensity.cat.categories]
                if len(not_in_results) > 0:
                    raise IndexError('{} are not light intensities found in the osar object.'.format(not_in_results))

                results = results[results.light_intensity\
                                  .isin(include_intensities)]
            else:
                raise ValueError('{} is not a list.'.format(include_intensities))

        if drop_nan_from is not None:
            if isinstance(drop_nan_from, list):
                # check given intensities are actually in the results.
                not_in_results = [a for a in drop_nan_from
                                  if a not in
                                  results.columns]

                if len(not_in_results) > 0:
                    raise IndexError('{} are not found in the columns of `osar.results`.'.format(not_in_results))

                drop_nan_from.append(x)
                drop_nan_from.append(y)
            else:
                raise ValueError('{} is not a list.'.format(drop_nan_from))
        else:
            drop_nan_from = [x, y]

        allcols = allowed_group_bys + drop_nan_from
        results = results[allcols].dropna()

        facets = results[group_by].cat.categories.tolist()
        flies = self.flies
        driver = self.driver

        if len (flies.opsin.unique()) > 1:
            raise ValueError('The experiment has more than one opsin.')

        opsin = flies.opsin.unique()[0]
        if opsin == 'Chrimson':
            colmap = 'OrRd'
            edgecol = 'red'
        elif opsin == 'GtACR1':
            colmap = 'YlGn'
            edgecol = 'green'

        # light_intensities = flies.light_intensity.cat.categories.tolist()
        light_intensities = results.light_intensity.unique().tolist()

        # Custom palettes
        pal_colors = sns.color_palette(colmap, len(light_intensities))
        pal = dict(zip(light_intensities, pal_colors))

        # Custom legends
        legend_elements = []

        for l in light_intensities:
            legend_elements.append(Circle((0,0),
                                           label=l,
                                           facecolor=pal[l],
                                           radius=10))


        results.set_index([group_by, 'light_intensity'],
                          inplace=True)

        sns.set(style='ticks', font_scale=2)
        collated_reg = []
        facet_count = len(facets)
        fig, axx = plt.subplots(nrows=facet_count,
                                figsize=(6, 7*facet_count),
                                gridspec_kw={'hspace':0.5})
        xlims=(-1, 1)
        ylims=(-1, 1)

        for j, facet in enumerate(facets):
            facet_data = results.loc[facet]
            plot_ax = axx[j]

            # Set 5% padding on x and ylims.
            plot_ax.set_xlim(np.array(xlims)*1.05)
            plot_ax.set_ylim(np.array(ylims)*1.05)

            for q, light_intensity in enumerate(light_intensities[::-1]):
                facet_slice = facet_data.loc[light_intensity]

                plot_ax.scatter(x=facet_slice[x],
                                y=facet_slice[y],
                                marker='o',
                                edgecolor=edgecol,
                                color=pal[light_intensity],
                                alpha=0.8,
                                zorder=q,
                                s=80)

            # Draw the regression line.
            reg_data = facet_data[[x, y]].dropna()
            x_reg = reg_data[x]
            y_reg = reg_data[y]
            r = sns.regplot(x_reg, y_reg,
                            ci=95, ax=plot_ax,
                            marker='x',
                            scatter=False,
                            line_kws={'lw': 1,
                                      'color': 'k'})

            # Compute the r2 and slope.
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                corr = sp.stats.linregress(x_reg, y_reg)
                lowers, uppers = skb.ci((x_reg,y_reg),
                                        statfunction=pth.r2_and_slope,
                                        n_samples=1000)

            r2 = np.round(corr.rvalue**2, decimal_places)
            r2_lower = np.round(lowers[0], decimal_places)
            r2_upper = np.round(uppers[0], decimal_places)

            slope = np.round(corr.slope, decimal_places)
            slope_lower = np.round(lowers[1], decimal_places)
            slope_upper = np.round(uppers[1], decimal_places)
            pval = np.round(corr.pvalue, decimal_places)

            annot_text = r'$R^{2}=$' + \
            '{} [95CI {}, {}]\n'.format(r2, r2_lower, r2_upper) + \
            'slope = {} [95CI {}, {}]\n'.format(slope, slope_lower, slope_upper) + \
            r'$n=$' + str(len(reg_data))

            collated_reg.append(
                            pd.Series([driver, opsin, facet,
                                       len(reg_data), x, y,
                                       r2, r2_lower, r2_upper,
                                       slope, slope_lower, slope_upper,
                                       corr.pvalue]))

            plot_ax.annotate(annot_text,
                             (1.05, -0.5),  # xy coords
                             ha='left', va='center',
                             textcoords='data',
                             fontsize=20)

            r.legend(title='Light Intensity',
                     handles=legend_elements,
                     loc='upper left',
                     bbox_to_anchor=(1.07,1))

            # plot_ax.set_ylabel('Preference Index')
            # plot_ax.set_xlabel(xlabel_dict[x_metric])

            subplot_title = "{} {} {}".format(driver, opsin, facet)

            plot_ax.set_title(subplot_title, fontweight='bold', fontsize=20)
            plot_ax.tick_params(length=15)
            sns.despine(ax=plot_ax, trim=True, offset=2)

        collated_reg = pd.DataFrame(collated_reg)
        collated_reg.columns=['driver', 'opsin', group_by,
                              'Ns', 'x_metric', 'y_metric',
                              'r2', 'r2_lower','r2_upper',
                              'slope', 'slope_lower', 'slope_upper',
                              'pvalue']
        collated_reg.reset_index(inplace=True, drop=True)

        return fig, collated_reg
