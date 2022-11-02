#!/usr/bin/python
# -*-coding: utf-8 -*-
# Author: Joses Ho
# Email : joseshowh@gmail.com

"""
Produces a dabest object functions for espresso objects.
"""

class log2_speed_ratio_plotter:
    """
    Speed ratio Produces a dabest objectting class for osar object.

    Available methods
    -----------------
    light_speed_against_dark_speed
    """



    def __init__(self,plotter):
        self.__results = plotter._results
        self.__driver = plotter._driver



    def light_speed_against_dark_speed(self,
                                       epoch='first',
                                       baseline_corrected=True,
                                       compare_by=None,
                                       group_by=None,
                                       color_by=None,
                                       distance_cutoff=15,
                                       cutoff_in_window=False,
                                       clean=False):
        '''
        Produces a dabest object for the log2 ratio of speed in light versus speed in dark.
        The log2 ratio will be computed per fly, and mean differences computed
        for each group in `group_by`.

        A positive log2 speed ratio indicates that the fly's speed in the light
        was faster relative to its speed in the dark.

        Keywords
        --------
        epoch: string, {'first','last','both'}, default 'first'
            Which epoch to plot. If 'both', the classical PI (i.e. the average
            of both epochs) will be plotted.

        baseline_corrected: boolean, default True
            If True, the baseline corrected metric is plotted.

        compare_by: string, default None
            Accepts a categorical column in the `osar.results`  object. This
            column will be used as the factor for generating and visualizing
            contrasts.

        group_by: string, default None
            Accepts a categorical column in the `osar.results` object. Each
            group in this column will receive its own 'hub-and-spoke' plot.

        color_by: string, default None
            Accepts a categorical column in the espresso object. Each group in
            this column will be colored seperately.

        palette_type: string, default 'categorical'
            Accepts 'categorical' or 'sequential'.
            
        distance_cutoff: float, default 15
            Any flies that travelled less than this amount in millimeters across
            the either the entire assay or the illuminated epoch
            (see `cutoff_in_window`) will not be plotted and included in the
            contrast. The default of 15mm should ensure that a fly will have
            encountered a light-dark border at least once.

        cutoff_in_window: boolean, default True
            If True, only considers the distance travelled during the
            illumination epochs. If False, the entire distance travelled
            throughout the assay is considered.

        clean: boolean, default False
            If True, will only plot flies with valid values (not NaNs) for three
            metrics: pi_smoothed_baseline_corrected,
                     light_attraction_index_baseline_corrected,
                     log2_speed_ratio



        Returns
        -------
        A matplotlib Figure depicting the contrast, and a pandas DataFrame with
        the estimation statistics.
        '''
        from . import plot_helpers as __pth

        if baseline_corrected is True:
            bc_suffix = '_baseline_corrected'
        else:
            if epoch == 'both':
                bc_suffix = '_not_baseline_corrected'
            else:
                bc_suffix = ''

        pi_dict = {'first':'{}_Pattern 01{}'.format("pi_smoothed", bc_suffix),
                   'last':'{}_Pattern 10{}'.format("pi_smoothed", bc_suffix),
                   'both':'{}{}'.format("pi_smoothed", bc_suffix)}

        lai_dict = {'first':'{}_Pattern 01{}'.format('light_attraction_index',
                                                     bc_suffix),
                      'last':'{}_Pattern 10{}'.format('light_attraction_index',
                                                     bc_suffix),
                      'both':'{}{}'.format('light_attraction_index', bc_suffix)}

        speed_dict = {'first':'{}_Pattern 01{}'.format('log2_speed_ratio',
                                                     bc_suffix),
                      'last':'{}_Pattern 10{}'.format('log2_speed_ratio',
                                                     bc_suffix),
                      'both':'log2_speed_ratio'}

        if epoch not in pi_dict.keys():
            raise KeyError('{} is not a valid `epoch` option.',format(epoch) +
                           '\nPlease enter "first", "last", or "both".')

        if clean:
            rr = self.__results
            plot_df = rr[(~rr[pi_dict[epoch]].isna()) &
                        (~rr[lai_dict[epoch]].isna()) &
                        (~rr[speed_dict[epoch]].isna())]
        else:
            plot_df = self.__results

        dabest_object = __pth.generic_dabest_parser(
                         plot_df=plot_df,
                         yvar=speed_dict[epoch],
                         driver=self.__driver,
                         opsin=self.__results.opsin[1],
                         compare_by=compare_by,
                         group_by=group_by,
                         color_by=color_by,
                         distance_cutoff=distance_cutoff)
        return dabest_object




    # def cross_to_dark(self,
    #                   compare_by=None,
    #                   group_by=None,
    #                   color_by=None,
    #                   palette_type='categorical',
    #                   **contrastplot_kwargs):
    #     '''
    #     Produces a dabest object for the log2 ratio of speed in light versus speed in dark,
    #     for traversals towards the dark. The log2 ratio will be computed per fly,
    #     and mean differences computed for each group in `group_by`.
    #
    #     Keywords
    #     --------
    #     compare_by: string, default None
    #         Accepts a categorical column in the `osar.results`  object. This
    #         column will be used as the factor for generating and visualizing
    #         contrasts.
    #
    #     group_by: string, default None
    #         Accepts a categorical column in the `osar.results` object. Each
    #         group in this column will receive its own 'hub-and-spoke' plot.
    #
    #     color_by: string, default None
    #         Accepts a categorical column in the espresso object. Each group in
    #         this column will be colored seperately.
    #
    #     palette_type: string, default 'categorical'
    #         Accepts 'categorical' or 'sequential'.
    #
    #     You can also pass along other keywords accepted by `plot` from the
    #     `dabest` package.
    #
    #     Returns
    #     -------
    #     A matplotlib Figure depicting the contrast, and a pandas DataFrame with
    #     the estimation statistics.
    #     '''
    #     from . import plot_helpers as __pth
    #
    #     dabest_object = __pth.generic_dabest_parser(
    #                      plot_df=self.__results,
    #                      yvar='log2_speed_ratio_traversal_dark',
    #                      driver=self.__driver,
    #                      opsin=self.__results.opsin[1],
    #                      compare_by=compare_by,
    #                      group_by=group_by,
    #                      color_by=color_by,
    #                      palette_type=palette_type,
    #                      **contrastplot_kwargs)
    #     return dabest_object
    #
    #
    #
    #
    #
    # def cross_to_light(self,
    #                    compare_by=None,
    #                    group_by=None,
    #                    color_by=None,
    #                    palette_type='categorical',
    #                    **contrastplot_kwargs):
    #     '''
    #     Produces a dabest object for the log2 ratio of speed in light versus speed in dark,
    #     for traversals towards the light. The log2 ratio will be computed per
    #     fly, and mean differences computed for each group in `group_by`.
    #
    #     Keywords
    #     --------
    #     compare_by: string, default None
    #         Accepts a categorical column in the `osar.results`  object. This
    #         column will be used as the factor for generating and visualizing
    #         contrasts.
    #
    #     group_by: string, default None
    #         Accepts a categorical column in the `osar.results` object. Each
    #         group in this column will receive its own 'hub-and-spoke' plot.
    #
    #     color_by: string, default None
    #         Accepts a categorical column in the espresso object. Each group in
    #         this column will be colored seperately.
    #
    #     palette_type: string, default 'categorical'
    #         Accepts 'categorical' or 'sequential'.
    #
    #     You can also pass along other keywords accepted by `contrastplot` from the
    #     `bootstrap_contrast` package.
    #
    #     Returns
    #     -------
    #     A matplotlib Figure depicting the contrast, and a pandas DataFrame with
    #     the estimation statistics.
    #     '''
    #     from . import plot_helpers as __pth
    #
    #     dabest_object = __pth.generic_dabest_parser(
    #                      plot_df=self.__results,
    #                      yvar='log2_speed_ratio_traversal_light',
    #                      driver=self.__driver,
    #                      opsin=self.__results.opsin[1],
    #                      compare_by=compare_by,
    #                      group_by=group_by,
    #                      color_by=color_by,
    #                      palette_type=palette_type,
    #                      **contrastplot_kwargs)
    #     return dabest_object
