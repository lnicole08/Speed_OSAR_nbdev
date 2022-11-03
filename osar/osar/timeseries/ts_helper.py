#!/usr/bin/python
# -*-coding: utf-8 -*-
# Author: Joses Ho
# Email : joseshowh@gmail.com
"""
Convenience Functions for timeseries plotting.
"""



def time_series_generic(timeseries, ax=None, ci=95, mean_axis=1,
                        error_ribbon_alpha=0.4, draw_zero_line=False,
                        **line_kwargs):
    """Convenience function to plot the mean of a timeseries with error ribbon.

    Keywords
    --------
    timeseries: pandas timeseries-like

    ax: matplotlib Axes object, default none

    ci: float, default 95
        The desired confidence interval that will be plotted as an error ribbon.

    mean_axis: 1 or 0
        The axis upon which to compute the mean.
        0 -- Compute column means
        1 -- Compute row means

    error_ribbon_alpha: float, default 0.4
        The alpha of the error ribbon.

    line_kwargs: matplotlib kwargs for Axes.plot() used to create the line
        timeseries.

    Returns
    -------
    matplotlib Axes
    """
    import matplotlib.pyplot as _plt
    from scipy.stats import norm as _norm

    if ci < 0.00001 or ci > 99.99999:
        err = '`ci` must between 1 and 99. The CI you entered ({}) is invalid.'.format(ci)
        raise ValueError(err)
    outer_bound = (1-(ci/100)) / 2
    zscore = ci/100 + outer_bound
    multiplier = _norm.ppf(zscore)

    if ax is None:
        ax = _plt.gca()

    mean = timeseries.mean(axis=mean_axis)
    half_ci = timeseries.sem(axis=mean_axis) * multiplier

    lower = mean - half_ci
    upper = mean + half_ci

    if 'lw' not in line_kwargs.keys() or 'linewidth' not in line_kwargs.keys():
        line_kwargs['linewidth'] = 1
    line_kwargs['zorder'] = 2

    ribbon_kwargs = {}
    if 'color' in line_kwargs.keys():
        ribbon_kwargs['color'] = line_kwargs['color']
    ribbon_kwargs['zorder'] = 1
    ribbon_kwargs['alpha'] = error_ribbon_alpha

    ax.plot(mean.index, mean, **line_kwargs)

    ax.fill_between(x=mean.index,
                    y1=lower,
                    y2=upper,**ribbon_kwargs)

    if draw_zero_line:
        ax.axhline(y=0, linewidth=0.25,
        color='k', zorder=0)

    return ax




def timeseries_plot(driver, meta, tracks,
                    framerate, pixel_length_mm,
                    rolling_winsize_seconds=5,
                    resample_winsize_seconds=0.1,
                    # smoothing_winsize_seconds=5,
                    show_illumination_epochs=True,
                    return_raw_data=False,
                    verbose=False):
    """
    Create timeseries plot.

    Keywords
    --------
    meta: DataFrame

    tracks: list of DataFrames

    framerate: int

    pixel_length_mm: float

    rolling_winsize_seconds: int, default 10
        The time window used to calculate the PI. By default, for each time-
        point, the PI for the last 10 seconds will be computed, for each fly.

    resample_winsize_seconds: float, default 0.1
        The new time interval to resample each fly's timeseries by. This is to
        allow proper calculation of the mean.

    show_illumination_epochs: boolean, default True
        If true, plot the illumination epochs as rectangles at the top of
        the PI axes.

    return_raw_data: boolean, default False
        If True, will return a pandas DataFrame of the rolling PI,
        and a pandas DataFrame of the rolling x-velocity.

    verbose: boolean, default False
        If True, the subplot titles will be printed to the console.

    Returns
    -------
    An array of matplotlib Axes objects, (and a pandas DataFrame of the rolling PI,
    and a pandas DataFrame of the rolling x-velocity, if return_raw_data is True).
    """

    import sys
    sys.path.append('..')

    import numpy as np
    import pandas as pd

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

    import seaborn as sns
    sns.set(style='ticks',context='poster')

    from ..munger import pi_array
    from ..plot_helpers import make_categorial_palette

    winsize_frames = framerate * rolling_winsize_seconds
    resample_period = '{}S'.format(resample_winsize_seconds)

    pi = []
    velocity = []

    if 'flyid' in meta.columns:
        meta.set_index('flyid', inplace=True, drop=True)

    for k, fly in enumerate(meta.index):
        sys.stdout.write('\rFly {} of {}'.format(k+1, len(meta)))

        # Isolate track of interest.
        tx = tracks[fly].copy()
        # change columns for backward compatibility.
        tx.columns = [a.lower() for a in tx.columns]

        # Get rolling window time spent (PI).
        rolling_pi = tx.lightstatus_new.rolling(window=winsize_frames)\
                                       .apply(pi_array)
        rolling_pi.index = pd.to_datetime(tx.seconds, unit='s')

        # Compute rolling window
        rolling_dist = tx.cx_smoothed_diff.abs().rolling(window=winsize_frames)\
                                                .apply(np.sum)
        rolling_dist.index = pd.to_datetime(tx.seconds, unit='s')

        # # Get velocity
        # velocity_abs = tx.cx_velocity.abs()
        # velocity_abs.index = pd.to_datetime(tx.seconds, unit='s')
        # velocity.append(velocity_abs)

        # Resample both PI and distance.
        resamp_kwargs = dict(rule='{}S'.format(resample_winsize_seconds),
                             closed='left',
                             label='right')
        pi_resamp = rolling_pi.resample(**resamp_kwargs).mean()
        dist_resamp = rolling_dist.resample(**resamp_kwargs).mean()
        speed_resamp = (dist_resamp * pixel_length_mm) / rolling_winsize_seconds

        pi.append(pi_resamp)
        velocity.append(speed_resamp)
        sys.stdout.flush()

    # Take the last tx DataFrame as indicative of light on and light off
    # timings. Analysis of various timelogs indicate that light on/off timings
    # have a tight spread of 0.2 seconds. So we can safely assume
    # (for visualization purposes) that the last `tx` is representative.
    lightstart = []
    lightend = []
    for e in [a for a in tx.experimentstate.unique()
             if a.startswith('Test ')]:
        lightstart.append(tx[tx.experimentstate == e].seconds.min())
        lightend.append(tx[tx.experimentstate == e].seconds.max())

    # Concatenate both lists into DataFrames.
    pi = pd.concat(pi, axis=1)
    velocity = pd.concat(velocity, axis=1)

    # Assign MultiIndex columns.
    cols = [meta.light_intensity, meta.status,
            meta.opsin, meta.driver,
            meta.light_color, meta.index]
    pi.columns = cols
    velocity.columns = cols

    # # Resample PI
    # pi_resamp = pi.resample(resample_period,
    #                         closed='left',
    #                         label='right').mean()
    #
    # # Resample velocity
    # velocity_resamp = velocity.resample(resample_period,
    #                                     closed='left',
    #                                     label='right').mean()

    groups = [(l,s) for l in meta.light_intensity.cat.categories
                    for s in meta.status.cat.categories]
    rowcount = len(meta.light_intensity.cat.categories)

    f, ax = plt.subplots(nrows=rowcount,
                         figsize=(11, 4.5*rowcount),
                         sharex=True,
                         # gridspec_kw={'hspace':0.5}
                         )

    # Create color palette.
    pal = make_categorial_palette(meta, group_by='status', pal='tab10')
    pal_dict = dict(zip(meta.status.cat.categories, pal))

    # Create axes.
    axes_pi = []
    axes_speed = []
    for k, a in enumerate(ax):
        # create new axes for plotting speed.
        ax_pi = ax[k]
        divider = make_axes_locatable(ax_pi)
        ax_speed = divider.append_axes("bottom",
                                       size="100%",
                                       pad="15%")
        axes_pi.append(ax_pi)
        axes_speed.append(ax_speed)
    ax_pi_dict = dict(zip(meta.light_intensity.cat.categories, axes_pi))
    ax_speed_dict = dict(zip(meta.light_intensity.cat.categories, axes_speed))

    # Create timeseries x-axis formatter.
    min_sec_fmt = mdates.DateFormatter('%M:%S')

    for k, g in enumerate(groups):
        # Get info from groups.
        light_intensity = g[0]
        status = g[1]
        label = "{} {}".format(g[0], g[1])
        if verbose:
            print(label)

        # Determine the correct axes to plot on.
        ax_pi = ax_pi_dict[light_intensity]
        ax_speed = ax_speed_dict[light_intensity]

        # # Smooth with rolling window.
        # smooth_win_size = int(smoothing_winsize_seconds/resample_winsize_seconds)
        # pi_rolling = pi_resamp[g].rolling(window=smooth_win_size).mean()
        # speed_rolling = velocity_resamp[g].rolling(window=smooth_win_size).mean()

        # Plot the timeseries.
        time_series_generic(pi[g], ax=ax_pi,
                            draw_zero_line=True,
                            color=pal_dict[status],
                            label=g[1])
        time_series_generic(velocity[g], ax=ax_speed,
                            draw_zero_line=False,
                            color=pal_dict[status],
                            label=g[1])

        # Plot the illumnation epochs.
        if show_illumination_epochs is True:
            light_color = meta.light_color[0]
            for i in range(len(lightstart)):
                s = pd.to_datetime(lightstart[i], unit='s')
                e = pd.to_datetime(lightend[i], unit='s')
                ax_pi.axvspan(xmin=s, xmax=e,
                              ymin=0.99, ymax=1.04,
                              color=light_color,
                              clip_on=False)


        # Set xlim for both axes.
        for a in [ax_pi, ax_speed]:
            a.set_xlim(pd.to_datetime('1970-01-01 00:00:00'),
                       pd.to_datetime('1970-01-01 00:02:40'))

        # Create legend.
        l = ax_pi.legend(loc='upper left',
                         bbox_to_anchor=(1., 0.99))
        # set the linewidth of each legend object
        for legobj in l.legendHandles:
            legobj.set_linewidth(1.5)

        # ax_pi aesthetics
        ax_pi.set_ylim(-1.05, 1.05)
        ax_pi.set_title('{}\n'.format(light_intensity))
        ax_pi.set_ylabel('PI')
        # Turn off the x-axis for `ax_pi`.
        plt.setp(ax_pi.get_xticklabels(), visible=False)
        ax_pi.xaxis.set_tick_params(size=0)
        sns.despine(ax=ax_pi, bottom=True, trim=True)

        # ax_speed aesthetics
        ax_speed.set_ylim(0, 8)
        ax_speed.set_yticks([0, 2, 4, 6, 8])
        ax_speed.set_ylabel('Speed\n(mm/s)')
        ax_speed.set_xlabel('Time (min:ss)')
        # Format the ax_speed xaxis properly.
        ax_speed.xaxis.set_major_locator(mdates.MinuteLocator())
        ax_speed.xaxis.set_major_formatter(min_sec_fmt)
        ax_speed.xaxis.set_minor_locator(mdates.SecondLocator(bysecond = [30]))
        ax_speed.xaxis.set_minor_formatter(min_sec_fmt)
        ax_speed.xaxis_date()
        sns.despine(ax=ax_speed, offset={'bottom':2})

    fig_title = '{}; window size = {}s, resampled to {}fps'.format(driver,
                                                rolling_winsize_seconds,
                                                1/resample_winsize_seconds)

    f.suptitle(fig_title, fontweight='bold', y=1.01)
    f.tight_layout()

    if return_raw_data:
        return ax, pi, velocity
    else:
        return ax
