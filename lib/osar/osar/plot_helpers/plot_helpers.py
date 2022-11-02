#!/usr/bin/python
# -*-coding: utf-8 -*-
# Author: Joses Ho
# Email : joseshowh@gmail.com




def normalize_ylims(ax_arr,include_zero=False,draw_zero_line=False):
    """Custom function to normalize ylims for an array of axes."""
    ymins=list()
    ymaxs=list()
    for ax in ax_arr:
        ymin=ax.get_ylim()[0]
        ymax=ax.get_ylim()[1]
        ymins.append(ymin)
        ymaxs.append(ymax)
    new_min=np.min(ymins)
    new_max=np.max(ymaxs)
    if include_zero:
        if new_max<0:
            new_max=0
        if new_min>0:
            new_min=0
    for ax in ax_arr:
        ax.set_ylim(new_min,new_max)
    if draw_zero_line:
        for ax in ax_arr:
            ax.axhline(y=0,linestyle='solid',linewidth=0.5,color='k')




def meanci(mean, cilow, cihigh, idx, ax, alpha=0.8, marker= 'o',color='black', size=8, ls='solid',lw=1.2):
    """Custom function to normalize plot the mean and CI as a dot and a vertical line, respectively."""
    # Plot the summary measure.
    ax.plot(idx, mean,
             marker=marker,
             markerfacecolor=color,
             markersize=size,
             alpha=alpha
            )
    # Plot the CI.
    ax.plot([idx, idx],
             [cilow, cihigh],
             color=color,
             alpha=alpha,
             linestyle=ls,
            linewidth=lw
            )




# Define function for string formatting of scientific notation.
def sci_nota(num, decimal_digits=2, precision=None, exponent=None):
    """
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.

    Found on https://stackoverflow.com/questions/21226868/superscript-in-python-plots
    """
    import numpy as np

    if not exponent:
        exponent = int(np.floor(np.log10(abs(num))))
    coeff = round(num / float(10**exponent), decimal_digits)
    if not precision:
        precision = decimal_digits

    return r"${0:.{2}f}\times10^{{{1:d}}}$".format(coeff, exponent, precision)





def make_categorial_palette(df, group_by, pal='Set1'):
    """
    Create a categorical color palette.
    Pass a pandas DataFrame and the column to group by.
    """
    import seaborn as sns

    _cat_palette = sns.color_palette(pal, n_colors=len(df[group_by].unique()))
    return _cat_palette




def make_sequential_palette(df, group_by):
    """
    Create a sequential color palette.
    Pass a pandas DataFrame and the column to group by.
    """
    import seaborn as sns

    _seq_palette = sns.cubehelix_palette(n_colors=len(df[group_by].unique()))
    return _seq_palette





def format_timecourse_xaxis(ax):
    """
    Convenience function to format a timecourse plot's x-axis.
    """
    import matplotlib.ticker as _tk

    ax.set_xlim(0,21600)
    ax.xaxis.set_ticks(range(0,25200,3600))
    ax.xaxis.set_minor_locator( _tk.MultipleLocator(base=1800) )
    ax.set_xlabel('Time (h)',fontsize=17)
    newlabels = [str(int(t/3600)) for t in ax.xaxis.get_ticklocs(minor=False)]
    ax.set_xticklabels(newlabels)
    ax.tick_params(axis='x', which='major',length=10)
    ax.tick_params(axis='x', which='minor',length=6)




def merge_two_dicts(x, y):
    """
    Given two dicts, merge them into a new dict as a shallow copy.
    Any overlapping keys in `y` will override the values in `x`.

    Taken from https://stackoverflow.com/questions/38987/
    how-to-merge-two-python-dictionaries-in-a-single-expression

    Keywords:
        x, y: dicts

    Returns:
        A dictionary containing a union of all keys in both original dicts.
    """
    z = x.copy()
    z.update(y)
    return z




def generic_dabest_parser(plot_df, yvar, driver, opsin,
                         compare_by=None,
                         group_by=None,
                         color_by=None,
                         distance_cutoff=15,
                         cutoff_in_window=False,
                         clean=False):
    """
    Generic dabest parser for osar objects.
    """
    import numpy as np

    import matplotlib.pyplot as plt
    plt.rcParams['svg.fonttype'] = 'none' # ensure all text renders as text.
    import dabest

    import sys as sys
    sys.path.append('..')
    from .. import munger as munger

    # Handle kwargs for OSAR plotting.
    if compare_by is None:
        compare_by = 'status'
    if group_by is None:
        group_by = 'light_intensity'
    if color_by is None:
        color_by = 'genotype'

    # Now, create 'groups' column as categorical.
    if compare_by == group_by:
        grp_contrast_cols = [group_by]
    else:
        grp_contrast_cols = [group_by, compare_by]

    # drop flies that walked a total distance less than distance_cutoff.
    if cutoff_in_window is True:
        if 'total_distance_window_mm' in plot_df.columns:
            plot_col = 'total_distance_window_mm'
        else:
            plot_col = 'total_distance_window'
    else:
        if 'total_distance_assay_mm' in plot_df.columns:
            plot_col = 'total_distance_assay_mm'
        else:
            plot_col = 'total_distance_assay'
    plot_df = plot_df[plot_df[plot_col] > distance_cutoff].copy()

    plot_df.sort_values(grp_contrast_cols,inplace=True)
    plot_df['groups_with_contrast'] = munger.join_cols(plot_df,
                                                         grp_contrast_cols)
    plot_df.loc[:,'groups_with_contrast'] = plot_df.groups_with_contrast.astype('category')
    plot_df.loc[:,'groups_with_contrast'] = plot_df.groups_with_contrast.cat.set_categories(
                                            [a for a in plot_df.groups_with_contrast.unique()],
                                            ordered=True)
    # Drop unused categories.
    categoricals = plot_df.dtypes[plot_df.dtypes == 'category'].index.tolist()
    for c in categoricals:
        plot_df[c].cat.remove_unused_categories(inplace=True)

    # # Handle palette.
    # if 'custom_palette' not in other_kwargs.keys():
    #     if palette_type == 'categorical':
    #         color_palette = make_categorial_palette(plot_df, color_by)
    #     elif palette_type == 'sequential':
    #         color_palette = make_sequential_palette(plot_df, color_by)
    #     else:
    #         raise ValueError('{} is not a recognized palette_type. Only `categorical` and `sequential` are recognized.'.format(palette_type))
    #     custom_palette = dict( zip( [a for a in plot_df[color_by].unique()],
    #                             color_palette )
    #                       )
    # 
    # else:
    #     custom_palette = other_kwargs['custom_palette']

    # Sort out the idx.
    idx=[ tuple(i) for i in np.array_split(
                    [a for a in plot_df.groups_with_contrast.unique()],
                    len(plot_df[group_by].unique())
            )
        ]
        
    # Create dabest object for plotting.
    cols_of_interest = [group_by, compare_by, color_by, yvar, 
                        "groups_with_contrast"]
    
    dabest_data = dabest.load(plot_df[cols_of_interest], idx=idx, 
                              x="groups_with_contrast", y=yvar)
                              
    return dabest_data
                              

    # es = getattr(dabest_data, effect_size)
    # 
    # # Assign default contrastplot keyword arguments.
    # default_plot_kwargs = dict( swarm_ylim=(-1.05,1.05),
    #                            fig_size=(15,7),
    #                            raw_marker_size=4,
    #                            swarm_label=yvar,
    #                            custom_palette=custom_palette )
    # 
    # default_plot_kwargs = merge_two_dicts(default_plot_kwargs, other_kwargs)
    # 
    # # Make the plot and compute contrasts.
    # fig, con = dabest.plot(**all_kwargs)
    # # Make the title.
    # fig.suptitle('{} {}'.format(driver, opsin),
    #             fontsize=20,
    #             fontweight='bold')
    # 
    # return fig, con


def r2_and_slope(x,y):
    '''
    Custom functions to produce r2 from linear regression stats.
    For computing 95 CIs.
    '''
    import scipy as sp
    reg = sp.stats.linregress(x,y)
    return reg.rvalue**2, reg.slope
