#!/usr/bin/python
# -*-coding: utf-8 -*-
# Author: Joses Ho
# Email : joseshowh@gmail.com

"""
Track plotting for OSAR objects.

Options available:
    single_fly_track
    fly_tracks_by_group
    fly_tracks
    border_plots
"""



def single_fly_track(osar,
                    fly,
                    normalize_max_speed=None,
#                   colormap='magma',
#                   colorbar=True,
                    # arena_length_mm=55,
                    # zone_count=4,
                    track_width=3,
                    ax=None,
                    show_lighted_regions=True,
                    annotate=True,
                    add_title=False,
                    hide_axes=False,
                    ):
    '''
    Given a `fly` from an `osar` object, its x-position will be plotted
    over time.

    Keywords
    --------
        osar: an `osar` object

        fly: string

        normalize_max_velocity: float, default None
            The maximum velocity (in mm/s) to set as the maximum of the colorbar.
            If None, defaults to the maximum velocity across the entire `osar`
            experiment (aka `osar.max_velocity_mm_per_s`).

        track_width: integer, default 3
            The point size of the line used to plot the x-position.

        show_lighted_regions: boolean, default True

        annotate: boolean, default True
            If True, annotate the track with detected traversals and reversals.

        add_title: boolean, default False
            If True, a title will be added to the subplot in the format:
            <fly> <genotype> <status>
                  <originalcsv>

        hide_axes: boolean, default False
            If True, both the x and y axes will be hidden.
            Simply applies the command
            matplotlib.AxesSubplot.<axis>.set_visible(False) to both axes.

    Returns
    -------
        A matplotlib AxesSubplot object.
    '''
    import numpy as np

    import matplotlib.pyplot as plt
    from matplotlib import ticker as tks
    from matplotlib.patches import Rectangle
    from matplotlib.collections import LineCollection
    from matplotlib.colors import ListedColormap, BoundaryNorm
    plt.rcParams['svg.fonttype'] = 'none' # ensure all text renders as text.

    import seaborn as sns
    sns.set(style='ticks', context='poster')

    try:
        # Get track for fly of interest.
        track = osar.tracks[fly].copy()
    except KeyError:
        raise KeyError('{} cannot be found in the smoothed tracks. '
                       'Please check.'.format(fly))

    # Munge the data to segment the positions clearly
    ## The following if-elif loop ensures legacy compatibility.
    if 'patternState' in track.columns:
        ps_col = 'patternState'
    elif 'PatternState' in track.columns:
        ps_col = 'PatternState'
    track.loc[:, 'segments'] = track[ps_col]
    min_idx = track[track[ps_col]=='Pattern 01'].index.min()
    track.loc[0:min_idx, 'segments'] = 'EXPT_START'

    # Get metadata of fly as a dict.
    try:
        if 'flyid' in osar.flies.columns:
            metadata = osar.flies[osar.flies.flyid == fly].T
        else:
            metadata = osar.flies.loc[fly].T
    except KeyError:
        raise KeyError('We cannot find {} in `flies`!'.format(fly))
    # metadata.columns = ['fly-of-interest']
    metadata = metadata.to_dict()#['fly-of-interest']
    # Have to add 0 and 55 to start and end of the mz arrays,
    mz = [np.insert(np.append(np.array(x) * osar.PIXEL_LENGTH_MM,55),
                              0,0)
          for x in metadata['midzones']]

    # quadrantLength = arena_length_mm / zone_count

    # Initialise plot.
    if ax is None:
        fig, ax = plt.subplots( figsize = (4,12) )
    # the figure width of 4in. seems to give the x-axis to the right scale.

    for p in track.segments.unique():
        temp_track = track[track.segments==p]
        # Define x and y.
        ## Following loop is for legacy compatibility.
        if 'seconds' in track.columns:
            seconds_col = 'seconds'
            ps_col = 'patternState'
            es_col = 'experimentState'
        elif 'Seconds' in track.columns:
            seconds_col = 'Seconds'
            ps_col = 'PatternState'
            es_col = 'ExperimentState'
        y = temp_track[seconds_col]
        try:
            x = temp_track.cX_smoothed * osar.PIXEL_LENGTH_MM
        except AttributeError:
            x = temp_track.cX_smooth * osar.PIXEL_LENGTH_MM
    #     dy_dx = track.velocity.abs() * osar.PIXEL_LENGTH_MM

        # Creating line segments.
        # See http://scipy-cookbook.readthedocs.io/items/Matplotlib_MulticoloredLine.html
        # Make each xy point a line segment.
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

    #     # Create the line collection object, setting the colormapping parameters.
    #     # Have to set the actual values used for colormapping separately.
    #     if normalize_max_speed is None:
    #         max_velocity=osar.max_velocity_mm_per_s
    #     else:
    #         max_velocity=np.float(normalize_max_speed)

        if p == 'Pattern 01':
            cmap = ListedColormap(['black', metadata['light_color'],
                                   'black', metadata['light_color']])
            norm = BoundaryNorm(mz[0], cmap.N)
        elif p == 'Pattern 10':
            cmap = ListedColormap([metadata['light_color'], 'black',
                                   metadata['light_color'], 'black'])
            norm = BoundaryNorm(mz[1], cmap.N)
        else:
            cmap = ListedColormap(['silver'])
            norm = BoundaryNorm([0], cmap.N)

        lc = LineCollection(segments, cmap=cmap, norm=norm)
        lc.set_array(x)
#         lc.set_array(dy_dx) # set the array to color the LineCollection
        lc.set_linewidth(track_width)
        ax.add_collection(lc)

    if annotate:
        # Plot choice status.
        decs = track.loc[~np.isnan(track.decisionCode)].index
        for d in decs:
            ax.annotate(track.decision.loc[d],
                        fontsize = 13,
                        xycoords='data',
                        horizontalalignment='left',
                        xy=((track.cX.loc[d]* osar.PIXEL_LENGTH_MM)+2,
                            track[seconds_col].loc[d]))

    if show_lighted_regions:
        lane = fly.split('_')[-1]
        BORDERS = osar.borders.loc[lane].copy()  * osar.PIXEL_LENGTH_MM

        # Draw boxed regions to indicated lighted epochs and regions.
        for ix, p in enumerate(['Pattern 10', 'Pattern 01']):
            current_borders = BORDERS.loc[p]
            p_slice = track[(track[ps_col] == p) &
                            (track[es_col] != 'Light off')]
            t_start = p_slice[seconds_col].min()
            t_end = p_slice[seconds_col].max()
            t_span = np.int(np.round(t_end - t_start))

            rect1_width = current_borders['midzone2'] - current_borders['midzone1']
            rect_y_adj = rect1_width * ix
            rect = Rectangle(# coordinates of bottom left corner
                            (0+rect_y_adj, t_start),
                            # width of rectangle, and its height, in data coords.
                            rect1_width, t_span,
                            edgecolor='w', lw=0,
                            fill=True, alpha=0.25,
                            facecolor=metadata['light_color'])
            ax.add_patch(rect)

            rect2_width = current_borders['midzone3'] - current_borders['midzone2']
            rect_y_adj = rect2_width * ix
            rect = Rectangle((current_borders['midzone2']+rect_y_adj, t_start),
                             rect2_width, t_span,
                             edgecolor='w',
                             lw=1,
                             fill=True, alpha=0.25,
                             facecolor=metadata['light_color'])
            ax.add_patch(rect)

    # Aesthetic adjustments.
    ax.set_ylim(161, 0)
    ax.set_xlim(0, 57)

    ax.set_xticks([0,5,15,25,35,45,55])
    ax.xaxis.set_minor_locator(tks.MultipleLocator(base=10.0))

    ax.yaxis.set_major_locator(tks.MultipleLocator(base=20.0))
    ax.yaxis.set_minor_locator(tks.MultipleLocator(base=5.0))

    # Label the plot.
    ax.set_ylabel('Time (s)')
    ax.set_xlabel('cX (mm)')

    sns.despine(ax = ax,
               offset = 10,
               trim = True,
               left=hide_axes,
               bottom=hide_axes)
    if hide_axes:
        ax.yaxis.set_visible(False)
        ax.xaxis.set_visible(False)

    if add_title:
        # Create title.
        originalcsv = metadata['original_count_log']
        genotype = metadata['genotype']
        status = metadata['status']
        title = fly+' '+genotype+' '+status+'\n'+originalcsv
        ax.set_title(title, fontsize = 10, y=1.025)

#     if colorbar:
#         # Create color bar
#         axBar = plt.colorbar(lc, orientation = 'horizontal',
#                              ax=ax,
#                              pad = 0.1,
#                              aspect = 15,
#                              ticks = tks.MultipleLocator(base=5.0))
#         axBar.set_label('Speed (mm/s)')

    return ax



def fly_tracks(osar, status, light_intensity,
               ncols=5, panel_width=1.5, panel_height=4,
               annotate=True, verbose=True,
               track_width=3,
               ax=None,
               show_lighted_regions=True,
               **subplot_kwargs):
    '''
    Given an `osar` object, and the relevant group of flies (via `status` and
    `light_intensity`), each fly's x-position will be
    plotted over time, and numbered

    Keywords
    --------
        osar: an `osar` object

        status: string, 'Offspring' or 'Sibling'

        light_intensity: string, or list.
            Any of these values: 'Eighth', 'Quarter', 'Half', 'Full'

        ncols: int, default 5
            The number of columns in the resulting grid of AxesSubplots.

        panel_width, panel_width: float, default 1.75 and 4 respectively.
            The size, in INCHES, of each panel in the grid.

        annotate: boolean, default True
            If True, will place a number above the track plot corresponding to
            the location in `osar_object.flies.index`. This is for easily pulling
            out interesting or problematic track plots for closer inspection.
        
        track_width: integer, default 3
            The point size of the line used to plot the x-position.

        show_lighted_regions: boolean, default True

        verbose: boolean, default True
            If True, will report the status of plotting.

        subplot_kwargs: any keyword arguments that can be passed to `plt.subplots`.

    Returns
    -------
        An array of matplotlib AxesSubplot objects.
    '''
    import sys
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    plt.rcParams['svg.fonttype'] = 'none' # ensure all text renders as text.

    import seaborn as sns
    sns.set(style='ticks',context='poster')

    if isinstance(light_intensity, str):
        list_of_fly_ids = osar.flies[(osar.flies.light_intensity == light_intensity) &
                                    (osar.flies.status == status)]\
                                    .index.tolist()

    elif isinstance(light_intensity, list):
        ff = osar.flies[osar.flies.light_intensity.isin(light_intensity)]
        list_of_fly_ids = ff[ff.status == status].index.tolist()

    total_axes_count = len(list_of_fly_ids)

    row_count = np.int(np.ceil(total_axes_count / ncols))

    f, a = plt.subplots(ncols=ncols,
                        nrows=row_count,
                        figsize=(ncols * panel_width,
                                 row_count * panel_height),
                        **subplot_kwargs
                        )

    for j, fly in enumerate(list_of_fly_ids):
        if verbose:
            sys.stdout.write('\r{} fly {} of {}'.format(fly,
                                                      j+1, total_axes_count))
        single_fly_track(osar, fly,
                         track_width=track_width,
                         show_lighted_regions=show_lighted_regions,
                         annotate=annotate,
                         add_title=False,
                         hide_axes=True,
                         ax=a.flat[j])
        a.flat[j].set_title(j)
        sys.stdout.flush()

    axes_to_turn_off = ncols * row_count - total_axes_count
    if axes_to_turn_off > 0:
        for ax in a.flat[-axes_to_turn_off::]:
            ax.axis('off')
    plt.tight_layout()

    driver = osar.driver
    f.suptitle('{}; {} flies, light intensity {}'.format(driver, status,
                                                         light_intensity),
               fontsize=20,
               fontweight='bold',
               y=1.01)

    return f



def __resample_for_border_plot(trajectory, resample_by='0.04S'):
    """
    Convenience function to produce a resampled timeseries for plotting.
    """
    import pandas as pd

    tt = trajectory.copy()
    tt.index = tt.index + 10
    tt.index = pd.to_datetime(tt.index, unit='s')
    tt_resamp = tt.resample(resample_by).mean()
    tt_resamp.fillna(method='bfill', inplace=True)
    tt_resamp.dropna(inplace=True)

    secs = tt_resamp.index.second
    msecs = tt_resamp.index.microsecond
    new_index = secs + (msecs*1e-6) - 10
    tt_resamp.index = new_index

    return tt_resamp



def __center_trajectory(track, PIXEL_LENGTH_MM,
                        light_border_mm, border_width,
                        x_direction, entry_frame,
                        fly, epoch, tx_id):
    import numpy as np
    tx = track.copy()

    tx.loc[:, 'Seconds_centered'] = tx.Seconds - tx.Seconds[entry_frame]

    tx['cX_smoothed_mm'] = tx.cX_smoothed * PIXEL_LENGTH_MM
    tx['cX_smoothed_mm_centered'] =  tx.cX_smoothed_mm - light_border_mm

    if x_direction == -1:
        tx.loc[:,'cX_smoothed_mm_centered'] =  -tx.cX_smoothed_mm_centered

    # For some reason, you have to filter out entries with cX < 10mmm,
    # and then recentre them.
    if tx.cX_smoothed_mm_centered[entry_frame] < -10:
        tx.loc[:,'cX_smoothed_mm_centered'] = tx.cX_smoothed_mm_centered + 10

    tx = tx[['Seconds_centered','cX_smoothed_mm_centered']].copy()
    tx.rename(columns={'cX_smoothed_mm_centered':
                        '{}_{}_tx{}'.format(fly, epoch, tx_id)},
                        inplace=True)
    tx.set_index('Seconds_centered',inplace=True)

    return tx


def interpolate_dropna(df):
    import pandas as pd

    d = df.copy()

    d.index = pd.to_datetime(d.index, unit='s')
    d_interp = d.interpolate(method='time')
    d_interp.index = d_interp.index.values.astype(int) / 1e9
    d_interp.dropna(inplace=True)

    return d_interp



def border_plots(to_plot, epoch,
                 light_intensities=['Half', 'Full'],
                 show_plot=True,
                 print_subplot_title=True,
                 time_before_zone_entry=1,
                 time_after_zone_entry=2,
                 distance_cutoff=5,
                 minimum_pixel_velocity=5,
                 fig_size=(25,12),
                 show_verge=True,
                 single_path_lw=0.7,
                 single_path_alpha=0.2,
                 single_path_color='black',
                 summary_func='mean',
                 show_summary_path=True,
                 summary_path_color=None,
                 summary_path_lw=2,
                 summary_path_alpha=0.9,
                 lighted_region_alpha=0.1,
                 distance_xlim=None,
                 verbose=True,
                 return_plotted_trajectories=False):
    """
    Plots every trajectory approaching the light-dark (and dark-light) border.

    Keywords
    --------
    to_plot: OSAR object.

    epoch: string.
        The illumination epoch to plot. Either 'Pattern 01' or 'Pattern 10'.

    light_intensities: list or string, default ['Half', 'Full'].
        The light intensities to include. If 'all', all light intensities are
        included.

    show_plot: boolean, default True
        If False, no plot will be produced.
        
    print_subplot_title: boolean, default True
        Whether to attach an info-rich title to each subplot.

    time_before_zone_entry, time_before_zone_entry: float or int,
                                        default 1 and 2 respectively.
        The time window to extract and plot, before and after the choice
        zone entry respecitvely.

    distance_cutoff: float, default 5
        Any flies that travelled less than this amount in millimeters across the
        illuminated epoch(s) will not be plotted and included in the analysis.

    minimum_pixel_velocity: float, default 5
        Any choice zone entries with a pixel velocity (pixels/s) less than this
        will be excluded from the analysis.

    show_verge: boolean, True.
        If True, the verge will be indicated with dotted lines.

    summary_func: 'mean' or 'median'

    show_summary_path: boolean, default True

    single_path_lw, single_path_alpha, summary_path_lw,
    summary_path_alpha: floats.
        Line width and alpha settings for the single paths and the mean paths.

    lighted_region_alpha: float, default 0.1
        The alpha (transparency) of the lighted region.

    distance_xlim: tuple, default None.
        The x-limits to plot. If None, the x-limits will be normalized to the
        minimum and maximum x-coordinates across the entire quartet.

    verbose: boolean, default True

    return_plotted_trajectories: boolean, False
        If True, returns the actual plotted trajectories as a dictionary of
        pandas DataFrames.

    Returns
    -------
    A matplotlib Figure, a dict, and a pandas DataFrame.

    The matplotlib Figure consists of a quartet of AxesSubplots. The first row
    depicts light-dark border approaches, and the second row shows dark-light
    border approaches.

    The dictionary returns the raw trajectories for each status (Sibling vs.
    Offspring) and verge type (Dark or Light) combination.

    The pandas DataFrame returns of Proportion of Exits to the Opposite Side.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import sys

    from .bootstrap_utils import summary_ci_1group
    from ..munger import is_on_other_side

    # Initialise mappers.
    light_color_opsin_mapper = {'Chrimson': 'red',
                                'GtACR1':'green'}
    summary_path_color_mapper = {'Chrimson': 'deeppink',
                                 'GtACR1':'lightseagreen'}
    verge_mapper = {'dark':0, 'light':1}

    EPOCH = epoch
    accepted_epochs = ['Pattern 01', 'Pattern 10']
    if isinstance(EPOCH, str):
        EPOCH = [EPOCH]
    for e in EPOCH:
        if e not in accepted_epochs:
            raise ValueError("{} is not recognized as one of {}".format(EPOCH,
                                                                accepted_epochs))

    # Properly handle the experiment(s).
    if isinstance(to_plot.osar_version, str) is False:
        raise ValueError("The object supplied to `to_plot` is not an OSAR object.")
    METADATA = to_plot.flies
    DRIVER = to_plot.driver
    RESULTS = to_plot.results.copy()

    # Handle light intensities.
    ACCEPTED_INTENSITIES = METADATA.light_intensity.cat.categories.tolist()
    if isinstance(light_intensities, str):
        if light_intensities in ACCEPTED_INTENSITIES:
            light_intensities = [light_intensities]
        else:
            raise ValueError('{} is not found in the metadata'.format(light_intensities))
    elif isinstance(light_intensities, list):
        if [a for a in light_intensities
            if a in ACCEPTED_INTENSITIES] != light_intensities:
            errstr = 'One or more intensities in {}'.format(light_intensities) +\
                'is not found in the metadata.'
            raise ValueError(errstr)

    FRAMEWINSIZE_LOW = int(time_before_zone_entry * to_plot.FRAME_RATE) + 1
    FRAMEWINSIZE_HIGH = int(time_after_zone_entry * to_plot.FRAME_RATE)
    OPSIN = to_plot.results.opsin.unique()[0]
    ROI_LIGHT_COLOR = light_color_opsin_mapper[OPSIN]
    
    if summary_path_color is None:
        SUMMARY_PATH_COLOR = summary_path_color_mapper[OPSIN]
    else:
        SUMMARY_PATH_COLOR = summary_path_color
        
    PIXEL_LENGTH_MM = to_plot.PIXEL_LENGTH_MM
    BORDER_WIDTH = to_plot.CHOICEZONE_WIDTH_MM / 2

    trajectories_dict = {}
    plot_dict = {}
    fly_count_dict = {}
    summary_list = []
    SUMMARY_COLUMNS = ['driver', 'opsin', 'genotype', 'light_intensity',
                       'status', 'expt_date', 'expt_time', 'border',
                       'flyid', 'epoch','trajectory_id']
    # To normalize the x-limits.
    # xmax = []
    # xmin = []

    ## START PLOT ##
    if show_plot is True:
        # sns.set(style='ticks', font_scale=1.5)
        f, axx = plt.subplots(nrows=2,ncols=2,
                              figsize=fig_size,
                              gridspec_kw={'hspace':0.6,
                                           'wspace':0.25})
    else:
        show_summary_path = False

    for rownum, vtype in enumerate(verge_mapper.keys()):
        row_facets = []
        VERGE_TYPE = vtype
        VERGE_TYPE_CODE = verge_mapper[VERGE_TYPE]


        for colnum, status in enumerate(METADATA.status.unique()):

            # Select for light intensity and status.
            flies_by_light = METADATA[
                            (METADATA.light_intensity.isin(light_intensities)) &
                            (METADATA.status == status)
                            ].index.tolist()

            all_tracks = []

            if show_plot is True:
                ax = axx[rownum, colnum]

                for_summ_plot = []

                # Normalize the x-limits
                if distance_xlim is None:
                    # xmin_arr = np.array(xmin)
                    # xmax_arr = np.array(xmax)
                    final_xmin = -26
                    final_xmax = 6
                else:
                    final_xmin = distance_xlim[0]
                    final_xmax = distance_xlim[1]

                ax.set_ylim(-time_before_zone_entry, time_after_zone_entry)
                ax.set_xlim(final_xmin, final_xmax)

            for E in EPOCH:
                # Filter by distance moved within epoch.
                total_dist_epoch = 'total_distance_{}'.format(E)
                RESULTS[total_dist_epoch] = RESULTS['distance_in_dark_{}'.format(E)] + \
                                            RESULTS['distance_in_light_{}'.format(E)]

                flies_by_locomotion = RESULTS[RESULTS[total_dist_epoch] >
                                                distance_cutoff].index.tolist()

                flies_to_plot = list(set(flies_by_locomotion)\
                                    .intersection(set(flies_by_light)))

                TOTAL_FLY_COUNT = len(flies_to_plot)

                if verbose:
                    print('\nPlotting {} {} {}, Total flies={}'.format(E, VERGE_TYPE,
                                                                       status, len(flies_to_plot)))

                for k, fly in enumerate(flies_to_plot):
                    if verbose:
                        sys.stdout.write('\rFly {} of {}'.format(k+1, TOTAL_FLY_COUNT))

                    # Get track.
                    tx = to_plot.tracks[fly]
                    # Get the FRAME corresponding to the END OF THE EPOCH.
                    last_frame = tx[tx.PatternState == E].iloc[-1].name
                    # Get flyid
                    flyid = fly.split('_')[-1]
                    # Extract all choice zone entries in EPOCH,
                    # with an x-velocity in pixels exceeding minimum_pixel_velocity.

                    # Extract all choice zone entries in EPOCH,
                    # with an x-velocity in pixels exceeding minimum_pixel_velocity.
                    all_choice_zone_entries = tx[(tx.choiceZoneStatus ==1 ) &
                                                 (tx.PatternState == E) &
                                                 (tx.cX_velocity.abs() > minimum_pixel_velocity)
                                                 ].copy()
                    # Move on to next fly if the there are no choice zone entries.
                    if len(all_choice_zone_entries) == 0:
                        TOTAL_FLY_COUNT -= 1
                    else:
                        # Now get all exits.
                        all_choice_zone_exits = tx[(tx.choiceZoneStatus == -1) &
                                                   (tx.PatternState == E)].copy()
                        # Make sure they line up.
                        if len(all_choice_zone_exits) < len(all_choice_zone_entries):
                            all_choice_zone_exits = all_choice_zone_exits.append(tx.loc[last_frame])

                        for df in [all_choice_zone_entries, all_choice_zone_exits]:
                            df.reset_index(inplace=True)
                            df.rename(columns={'index':'frame'},inplace=True)

                        verge_type_choice_zone_entries = \
                            all_choice_zone_entries[all_choice_zone_entries.lightStatus_new==VERGE_TYPE_CODE]
                        verge_type_choice_zone_exits = \
                            all_choice_zone_exits.loc[verge_type_choice_zone_entries.index]

                        # Loop through each choice zone entry, pull out
                        # the desired window after the choice zone entry from dark.
                        for l, idx in enumerate(verge_type_choice_zone_entries.index):
                            entry_frame = verge_type_choice_zone_entries.loc[idx, 'frame']
                            exit_frame = verge_type_choice_zone_exits.loc[idx, 'frame']
                            current_track = tx[entry_frame-FRAMEWINSIZE_LOW:
                                               exit_frame].copy()

                            if exit_frame - entry_frame > 2:
                                if show_plot:
                                    plot_track = tx[entry_frame - FRAMEWINSIZE_LOW:
                                                    entry_frame + FRAMEWINSIZE_HIGH]
                                # Figure out which border we need to pull out.
                                x_direction = verge_type_choice_zone_entries.loc[idx, 'cX_direction']

                                if x_direction == -1:
                                    zone = int(verge_type_choice_zone_entries.zone_new[idx])
                                else:
                                    zone = int(verge_type_choice_zone_entries.zone_new[idx]) + 1

                                if zone in [1,2,3]:
                                    light_border_px = to_plot.borders.loc[flyid]\
                                                        .loc[E]\
                                                        .loc['midzone{}'.format(zone)]
                                    light_border_mm = light_border_px * PIXEL_LENGTH_MM

                                    tx_id = l + 1
                                    centered_track = __center_trajectory(
                                                    current_track, PIXEL_LENGTH_MM,
                                                    light_border_mm, BORDER_WIDTH,
                                                    x_direction, entry_frame,
                                                    fly, E, tx_id)
                                    all_tracks.append(centered_track)

                                    if show_plot is True or show_summary_path is True:
                                        centered_plot = __center_trajectory(
                                                        plot_track, PIXEL_LENGTH_MM,
                                                        light_border_mm, BORDER_WIDTH,
                                                        x_direction, entry_frame,
                                                        fly, E, tx_id)
                                        for_summ_plot.append(centered_plot)

                                    if show_plot is True:
                                        p = centered_plot[(centered_plot.index > -time_before_zone_entry) &
                                                          (centered_plot.index < time_after_zone_entry)].copy()

                                        p = p[(p.values > final_xmin) &
                                              (p.values < final_xmax)]

                                        ax.plot(p.values, p.index,
                                                lw=single_path_lw,
                                                alpha=single_path_alpha,
                                                c=single_path_color)

            tracks = pd.concat(all_tracks, sort=False, axis='columns')

            if show_summary_path is True:
                tx = pd.concat(for_summ_plot, sort=False, axis='columns')

                tx_interp_dropna = interpolate_dropna(tx)

                if summary_func == 'mean':
                    summ_plot_line = tx_interp_dropna.mean(axis=1)

                elif summary_func == 'median':
                    summ_plot_line = tx_interp_dropna.median(axis=1)

                else:
                    errstr = '{} is not a recognized summary function.'\
                             .format(summary_func)
                    raise ValueError(errstr)

                # summ_plot_line.index = summ_plot_line.index.values.astype(int) / 1e9

                ax.plot(summ_plot_line.values, summ_plot_line.index,
                        lw=summary_path_lw, alpha=summary_path_alpha,
                        c=SUMMARY_PATH_COLOR)

            current_facet = "{} {} {}".format(VERGE_TYPE, DRIVER, status)
            row_facets.append(current_facet)

            # Save the plotted and centred tracks.
            trajectories_dict[current_facet] = interpolate_dropna(tracks)
            fly_count_dict[current_facet] = TOTAL_FLY_COUNT
            if show_plot:
                plot_dict[current_facet] = for_summ_plot
            # # Get max and min for this set of tracks.
            # xmax.append(tracks.max().max())
            # xmin.append(tracks.min().min())

            # Create summary of current tracks
            txx = tracks.copy()
            txx.columns = pd.MultiIndex.from_tuples(
                            txx.columns.str.split('_').tolist(),
                            names=SUMMARY_COLUMNS).droplevel(level='border')
            get_final_pos = txx.apply(is_on_other_side, axis='rows').astype(int)
            get_final_pos.name = 'exit_to_other_side'
            final_pos = get_final_pos.reset_index()
            final_pos['verge_type'] = np.repeat(VERGE_TYPE, len(final_pos))
            summary_list.append(final_pos)

            if show_plot:
                if show_verge:
                    for b in [-BORDER_WIDTH, BORDER_WIDTH]:
                        ax.axvline(x=b, lw=1, color='k', ls='dashed')
                ax.invert_yaxis()

                ax.set_ylabel('Time (s)')
                ax.set_xlabel('Distance from border (mm)')

                final_states = tracks.apply(is_on_other_side, axis='rows')
                prop_exits = summary_ci_1group(final_states.astype(int).values,
                                               np.mean)
                prop_string = "Prop. Exits to Other Side={} [{}; {}]"\
                                .format(np.round(prop_exits['summary'], 3),
                                        np.round(prop_exits['bca_ci_low'], 3),
                                        np.round(prop_exits['bca_ci_high'], 3))

                light_intensities_title = ', '.join(light_intensities)
                title = "{}; N={}, Trajectories={};\nEpoch={}; Intensities={}; {}"\
                            .format(current_facet, TOTAL_FLY_COUNT,
                                    len(tracks.columns),
                                    EPOCH, light_intensities_title,
                                    prop_string)
                if print_subplot_title:
                    ax.set_title(title, pad=10, fontsize=14, fontweight='bold')

    if show_plot:
        # Normalize the x-limits
        if distance_xlim is not None:
            final_xmin = distance_xlim[0]
            final_xmax = distance_xlim[1]

        axvspan_kwargs = {'alpha':lighted_region_alpha,
                          'color':ROI_LIGHT_COLOR}
        for n, ax in enumerate(axx.flatten()):
            # ax.set_xlim(final_xmin, final_xmax)
            if n < 2:
                ax.axvspan(xmin=0, xmax=final_xmax, **axvspan_kwargs)
                ax.axvspan(xmin=final_xmin, xmax=-13.75, **axvspan_kwargs)
                # if final_xmin < -13.75:
                #     ax.axvspan(xmin=-13.75*2, xmax=-13.75, **axvspan_kwargs)
            else:
                ax.axvspan(xmin=-13.75, xmax=0, **axvspan_kwargs)
                # if final_xmin < -13.75:
                #     ax.axvspan(xmin=-13.75, xmax=0, **axvspan_kwargs)
                # else:
                #     ax.axvspan(xmin=final_xmin, xmax=0, **axvspan_kwargs)
                # ax.axvspan(xmin=13.75, xmax=13.75*2, **axvspan_kwargs)
            sns.despine(ax=ax, trim=True, offset={'bottom':10, 'left':5})

    # Munge Summary
    summary = pd.concat(summary_list).reset_index()
    summary.loc[:,'flyid'] = summary.expt_date.str.cat(
                                summary.expt_time.str.cat(summary.flyid,
                                                          sep=' '),
                                sep=' ')
    summary.drop(['index', 'expt_date', 'expt_time'], axis=1, inplace=True)

    if show_plot:
        if return_plotted_trajectories:
            return f, trajectories_dict, summary, plot_dict
        else:
            return f, trajectories_dict, summary
    else:
        return trajectories_dict, summary




# def __get_len_post_dropna(s):
#     return len(s.dropna())
#
# def __normalize_trajectory_indexes(df, index_type):
#     """
#     Convenience function to equalize all the indexes for a trajectory dataframe.
#
#     index type: string
#         Either 'shortest' or 'first'.
#     """
#     import numpy as np
#     import pandas as pd
#
#     # Get trajectory with greatest loss after dropna
#     benchmark_col = df.apply(__get_len_post_dropna, axis=0).idxmin()
#
#     if index_type == 'shortest':
#         first_idx = df.loc[:,benchmark_col].dropna().index
#     elif index_type == 'first':
#         first_idx = df.loc[:,df.columns[0]].dropna().index
#
#     round_idx = np.round(first_idx, 2)
#
#     t = []
#     for i, c in enumerate(df.columns):
#         a = df.loc[:, c].dropna()
#         a.index = round_idx
#         t.append(a)
#
#     out = pd.concat(t, axis=1)
#     return out
