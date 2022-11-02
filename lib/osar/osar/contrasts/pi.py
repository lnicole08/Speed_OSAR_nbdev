#!/usr/bin/python
# -*-coding: utf-8 -*-
# Author: Joses Ho
# Email : joseshowh@gmail.com

"""
PI Produces a dabest object functions for osar objects.
"""


class pi_plotter(object):
    """
    PI Produces a dabest objectting class for espresso object.

    Available methods
    -----------------
    time_spent_in_light
    traversals_towards_light
    reversals_towards_light
    """

    #    #    #    #    #####
    #    ##   #    #      #
    #    # #  #    #      #
    #    #  # #    #      #
    #    #   ##    #      #
    #    #    #    #      #

    def __init__(self,plotter):
        self.__results = plotter._results
        self.__driver = plotter._driver



    def time_spent_in_light(self,
                            epoch='both',
                            baseline_corrected=True,
                            compare_by=None,
                            group_by=None,
                            color_by=None,
                            distance_cutoff=15,
                            cutoff_in_window=False,
                            clean=False):
        '''
        Produces a dabest object for the preference index for time spent in light. The
        PI for light will be computed per fly, and mean differences computed for
        each group in `group_by`.

        Keywords
        --------
        epoch: string, {'first','last','both'}, default 'first'
            Which epoch to plot. If 'both', the classical PI (i.e. the average
            of both epochs) will be plotted.

        baseline_corrected: boolean, default True
            If True, the baseline corrected value is plotted.

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
        A dabest object.
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
                                     yvar=pi_dict[epoch],
                                     driver=self.__driver,
                                     opsin=self.__results.opsin[1],
                                     compare_by=compare_by,
                                     group_by=group_by,
                                     color_by=color_by,
                                     distance_cutoff=distance_cutoff,
                                     cutoff_in_window=cutoff_in_window,
                                     clean=clean)
        return dabest_object



    def light_attraction_index(self,
                                epoch='both',
                                baseline_corrected=True,
                                compare_by=None,
                                group_by=None,
                                color_by=None,
                                palette_type='categorical',
                                distance_cutoff=15,
                                cutoff_in_window=False,
                                clean=False):
        '''
        Produces a dabest object for the light attraction index.  This will be computed
        per fly, and mean differences computed for each group in `group_by`.

        Keywords
        --------
        epoch: string, {'first','last','both'}, default 'first'
            Which epoch to plot. If 'both', the classical PI (i.e. the average
            of both epochs) will be plotted.

        baseline_corrected: boolean, default True
            If True, the baseline corrected value is plotted.

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
        A dabest object.
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
                                     yvar=lai_dict[epoch],
                                     driver=self.__driver,
                                     opsin=self.__results.opsin[1],
                                     compare_by=compare_by,
                                     group_by=group_by,
                                     color_by=color_by,
                                     distance_cutoff=distance_cutoff,
                                     cutoff_in_window=cutoff_in_window,
                                     clean=clean)
        return dabest_object



    def traversals_towards_light(self,
                                 compare_by=None,
                                 group_by=None,
                                 color_by=None,
                                 palette_type='categorical',
                                 distance_cutoff=15,
                                 cutoff_in_window=True):
        '''
        Produces a dabest object for the preference index for traversals to light. The
        PI will be computed per fly, and mean differences computed for each
        group in `group_by`.

        Keywords
        --------
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

        Returns
        -------
        A dabest object.
        '''
        from . import plot_helpers as __pth

        dabest_object = __pth.generic_dabest_parser(
                                    plot_df=self.__results,
                                    yvar='traversal_light_pi',
                                    driver=self.__driver,
                                    opsin=self.__results.opsin[1],
                                    compare_by=compare_by,
                                    group_by=group_by,
                                    color_by=color_by,
                                    distance_cutoff=distance_cutoff,
                                    cutoff_in_window=cutoff_in_window)
        return dabest_object




    def reversals_towards_light(self,
                                compare_by=None,
                                group_by=None,
                                color_by=None,
                                palette_type='categorical',
                                distance_cutoff=15,
                                cutoff_in_window=True):
        '''
        Produces a dabest object for the preference index for traversals to dark. The
        PI will be computed per fly, and mean differences computed for each
        group in `group_by`.

        Keywords
        --------
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

        You can also pass along any other keyword that is accepted by the
        `contrastplot` function from the `bootstrap_contrast` package.

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

        Returns
        -------
        A dabest object.
        '''
        from . import plot_helpers as __pth

        dabest_object = __pth.generic_dabest_parser(
                                    plot_df=self.__results,
                                    yvar='reversal_light_pi',
                                    driver=self.__driver,
                                    opsin=self.__results.opsin[1],
                                    compare_by=compare_by,
                                    group_by=group_by,
                                    color_by=color_by,
                                    distance_cutoff=distance_cutoff,
                                    cutoff_in_window=cutoff_in_window)
        return dabest_object
