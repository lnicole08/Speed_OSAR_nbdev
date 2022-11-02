#!/usr/bin/python
# -*-coding: utf-8 -*-
# Author: Joses Ho
# Email : joseshowh@gmail.com

def generic_heatmap_plotter(osar, colormap = 'magma'):
    '''generic docstring'''

    import sys

    import numpy as np
    import scipy as sp
    import scipy.interpolate as ipl

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    plt.rcParams['svg.fonttype'] = 'none' # ensure all text renders as text.

    import pandas as pd
    import seaborn as sns
    sns.set_style('ticks')

    # Set matplotlib aesthetic settings.
    SMALL_SIZE = 17
    MEDIUM_SIZE = 20
    BIGGER_SIZE = 23

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels

    plt.rc('xtick', direction='out')
    plt.rc('ytick', direction='out')

    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    plt.rc('ytick.major', size=12)            # ytick length in points.
    plt.rc('xtick.major', size=12)            # ytick length in points.

    plt.rcParams['svg.fonttype'] = 'none'
    mpl.rcParams['pdf.fonttype'] = 42

    # Edges to bin the histograms against.
    edges = np.round( np.array(np.arange(-5, 5.2, 0.1)),
                      decimals = 1 ) # we round  up the edges to get nice labels.
    binlength = np.diff(edges).mean().round(1) # calculate the bin length.

    # Create dicts of args used across all plots.
    histo_args = {'density':True,
                  'bins':edges }
    heatmap_cx_args = { 'vmin':0, 'vmax':1,
                      'cmap': colormap, # choose between viridis, plasma, inferno, and magma.
                      'yticklabels':False,
                      'xticklabels':10 }

    heatmap_speed_args={ 'vmin':0, 'vmax':25,
                        'cmap': colormap, # choose between viridis, plasma, inferno, and magma.
                        'yticklabels':False,
                        'xticklabels':10 }
    ## MAINPLOT ARGS
    gs_args={'wspace':0.2,
             'hspace':0.05,
             'width_ratios':[8,8,0.5], 'height_ratios':[9,1.5]}
    mainplot_args={'nrows':2, 'ncols':3,
                   'figsize':(20,10), 'gridspec_kw':gs_args}

    # args for summary histogram and the subsequent plot.
    summ_hist_args={'bins':edges,'density':True}

    summ_plot_cx_args={'x':edges[:-1], 'alpha':0.5, 'color':'black',
                     'linewidth':0.5,'antialiased':True}

    summ_plot_speed_args={'style':'-', 'color':'black',
                          'linewidth':0.5, 'antialiased':True}

    # Assign color appropriately.
    light_color = osar.results.light_color.unique()[0].lower()
    hspan_args={'facecolor':light_color,
                'edgecolor':light_color,
                'clip_on':False}

    # plot by light intensities
    all_intensities = [ l for l in osar.results.light_intensity.unique() ]

    out = []

    for intensity in all_intensities:
        sibiling_list = [k for k in osar.tracks.keys() if 'Sibling' in k]
        offspring_list = [k for k in osar.tracks.keys() if 'Offspring' in k]

        # Compile the lists  into a dict.
        all_hists={'sibling' : sibiling_list,
                   'offspring' : offspring_list}
        #### MAIN PLOTTING CODE ####
        problems=list()
        # all_df=list()
        all_hist=list()

        for status, l in all_hists.items(): # iterates over each key-pair in the dict `all_hists`.'
            sys.stdout.write('\rprocessing '+status)
            sys.stdout.flush()

            # Create lists to compute summary heatmaps.
            start_light_cx_summary=list()
            start_dark_cx_summary=list()
            start_light_speed_summary=list()
            start_dark_speed_summary=list()
            # Create lists to store all the per_fly histograms.
            start_light_cx_all=list()
            start_dark_cx_all=list()
            start_light_speed_all=list()
            start_dark_speed_all=list()

            for j, track in enumerate(l): # loop thru the CSVs in each csvlist (either sibling or offspring)
                sys.stdout.write('\rprocessing track '+track+' '+str(j+1)+' of '+str(len(l))+'\t')
                sys.stdout.flush()
                # Read the current track.
                ## track=full_offspring_csvlist[10] ## FOR TESTING ONLY
                df = osar.tracks[track]
                # df=pd.read_csv(os.path.join(dataDir,track))
                # Identify where entries and exits occur.
                entries=df[df.choiceZoneStatus==1].index
                exits=df[df.choiceZoneStatus==-1].index
                # These two indexes MUST and should have the same length.
                if len(entries)!=len(exits):
                    problems.append(track)
                    continue # If the entries and exits don't match, take note of it and skip to next track.
                # create lists to store cX and speed from when choice zone is occupied.
                start_light_cx_per_fly=list()
                start_light_speed_per_fly=list() # when the fly entered choice zone from the light.
                start_dark_cx_per_fly=list()
                start_dark_speed_per_fly=list() # when the fly entered choice zone from the dark.
                for idx, en in enumerate(entries):
                    # Loop through each choice zone.
                    zoomdf=df.ix[entries[idx]:exits[idx],:]
                    # Identify the x-coordinate of crossing the light-dark border.
                    crossing=zoomdf[~np.isnan(zoomdf.crossings)].cX_smoothed
                    if len(crossing)>0:
                        # Mark this point as zero, everything else centred around it.
                        cxx=( zoomdf.cX_smoothed-crossing.mean() ) * osar.PIXEL_LENGTH_MM
                        # If the fly is walking backwards into the choice zone, reverse it.
                        if zoomdf.cX_direction.iloc[0]==-1:
                            cxx=-cxx
                        ## Process speed! Incl. interpolation.
                        speed=np.abs(zoomdf.cX_velocity) * osar.PIXEL_LENGTH_MM, # The absolute speed, in mm/sec.
                        speed=speed[0] # Not sure why the above gives a tuple, of which the first element is what we want.
                        # Create interpolating function for speed.
                        f=ipl.interp1d(x=cxx,y=speed,kind='slinear',bounds_error=False)
                        # Here, we compute the range of interpolation from the max and min of the cxx.
                        minx=np.round(cxx.min(),decimals=1)
                        maxx=np.round(cxx.max(),decimals=1)
                        plotrange=np.round( np.arange(minx,maxx,0.1), decimals=1 )
                        # Perform interpolation.
                        interpol_speed=f(plotrange)
                        # Convert the interpolation to a pandas series with the index as the plotrange.
                        interpol_speed=pd.Series(interpol_speed,index=plotrange)
                    else:
                        continue # ignore this choice zone entry if no border crossings were made.
                    # Append to the correct list.
                    if zoomdf.lightStatus_new.iloc[0]==0:
                        start_dark_cx_per_fly.append(cxx)
                        start_dark_speed_per_fly.append(interpol_speed)
                    elif zoomdf.lightStatus_new.iloc[0]==1:
                        start_light_cx_per_fly.append(cxx)
                        start_light_speed_per_fly.append(interpol_speed)

                ### PROCESSING CX
                # flatten out cx lists, convert to pandas series, and dropnans.
                start_light_cx_per_fly=pd.Series( [cx for series in start_light_cx_per_fly for cx in series] ).dropna()
                start_dark_cx_per_fly=pd.Series( [cx for series in start_dark_cx_per_fly for cx in series] ).dropna()
                # Save these flattened lists to the respective summary lists.
                start_light_cx_summary.append(start_light_cx_per_fly)
                start_dark_cx_summary.append(start_dark_cx_per_fly)
                # Compute the histograms.
                start_light_hist,_=np.histogram(pd.Series(start_light_cx_per_fly),
                                                **histo_args)
                start_dark_hist,_=np.histogram(start_dark_cx_per_fly,
                                               **histo_args)
                # Divide by number of 1/bin length so each histogram sums to 1.
                start_light_hist=start_light_hist/(1/binlength)
                start_dark_hist=start_dark_hist/(1/binlength)
                # Save each histogram to the respective list.
                start_light_cx_all.append(start_light_hist)
                start_dark_cx_all.append(start_dark_hist)

                ### PROCESSING SPEED.
                # Compute mean speed profile for each fly, save to the respective list.
                start_light_speed_all.append(pd.DataFrame(start_light_speed_per_fly).mean(axis=0))
                start_dark_speed_all.append(pd.DataFrame(start_dark_speed_per_fly).mean(axis=0))

                # Indicate that we are done with this track.
                sys.stdout.write('\rprocessed track '+track+' '+str(j+1)+' of '+str(len(l))+'\t')
                sys.stdout.flush()

            # flatten out summary lists, save as pandas series for easy plotting.
            start_light_cx_summary=pd.Series( [cx for series in start_light_cx_summary for cx in series] )
            start_dark_cx_summary=pd.Series( [cx for series in start_dark_cx_summary for cx in series] )
            # convert lists to pandas DataFrame
            ## cx
            start_light_cx_all=pd.DataFrame(start_light_cx_all)
            start_dark_cx_all=pd.DataFrame(start_dark_cx_all)
            ## speed
            start_light_speed_all=pd.DataFrame(start_light_speed_all)
            start_dark_speed_all=pd.DataFrame(start_dark_speed_all)
            ## print max speed so I can tweak the colormap....
            print( '\nmax light speed (mm/s):', start_light_speed_all.max().max().round(1) )
            print( 'max dark speed (mm/s):', start_dark_speed_all.max().max().round(1) )
            # sort each df by the 1cm bin after the border, NaNs at the end.
            for df in [start_light_cx_all,start_dark_cx_all,
                       start_light_speed_all, start_dark_speed_all
                      ]:
                df.sort_values(axis=0, by=[0.],
                               ascending=False,
                               inplace=True,
                               na_position='last')
            # Assign appropriate column names for proper plotting.
            for df in [start_light_cx_all,start_dark_cx_all]:
                df.columns=edges[:-1] # drop the last edge so the bin-labels and column count matches.


            ##### Create the plot for cx occupancy.
            summ_hist_cx_axes=list()
            fig1,((axlight,axdark,axcolbar),
                 (axlightsum,axdarksum,xtra))=plt.subplots(**mainplot_args)
            xtra.axis('off') # We have created an extra axes for spacing; hide it!
            # Plot heatmaps for border approaches from light.
            sns.heatmap(start_light_cx_all,
                        cbar=False,ax=axlight,**heatmap_cx_args)
            axlight.xaxis.set_visible(False) # Hide xaxis for aesthetics and clarity.
            # Set the title
            axlight.text(.5,1.04,
                         'From Light',
                         fontsize=18,
                         horizontalalignment='center',
                         transform=axlight.transAxes)
            # Compute histogram for summary for border approaches from light.
            hlight,_=np.histogram(start_light_cx_summary,**summ_hist_args)
            hlightplot=pd.Series(hlight/(1/binlength),
                                 index=edges[:-1])
            # Plot the summary histogram for light approaches.
            axlightsum.fill_between(y1=hlightplot,**summ_plot_cx_args)
            axlightsum.set_xlim(start_light_cx_all.columns[0],
                               start_light_cx_all.columns[-1])
            sns.despine(ax=axlightsum,left=True)
        #     axlightsum.yaxis.set_visible(False)
            axlightsum.set_ylabel('\t\t')
            axlightsum.set_xlabel('Distance from\nlight-dark border (mm)')
            summ_hist_cx_axes.append(axlightsum) # Save axes into list for future ylim equalization.
            # Plot heatmaps for border approaches from dark.
            sns.heatmap(start_dark_cx_all,
                        cbar_ax=axcolbar,ax=axdark,**heatmap_cx_args)
            axcolbar.set_ylabel('Normalised Occupancy',labelpad=20,rotation=270) # Add title to our colorbar.
            axdark.xaxis.set_visible(False)
            axdark.text(.5,1.04,
                        'From Dark',
                        fontsize=18,
                        horizontalalignment='center',
                        transform=axdark.transAxes)
            # Compute histogram for summary for border approaches from dark.
            hdark,_=np.histogram(start_dark_cx_summary,**summ_hist_args)
            hdarkplot=pd.Series(hdark/(1/binlength),
                                 index=edges[:-1])
            # Plot the summary histogram for dark approaches.
            axdarksum.fill_between(y1=hdarkplot,**summ_plot_cx_args)
            axdarksum.set_xlim(start_dark_cx_all.columns[0],
                               start_dark_cx_all.columns[-1])
            sns.despine(ax=axdarksum,left=True)
        #     axdarksum.yaxis.set_visible(False)
            axdarksum.set_xlabel('Distance from\nlight-dark border (mm)')
            summ_hist_cx_axes.append(axdarksum)
            # Aesthetics commmon to both axes.
            for a in [axlight, axdark]:
                # Draw reference line.
                a.axvline(x=np.where(edges[:-1]==0.)[0][0],
                          linewidth=0.55,
                          color='white')
            # Draw reference hbar to indicate lighted region.
            axlightymax=axlight.get_ylim()[1]
            axlight.axhspan(xmax=0.5,
                            ymin=axlightymax+0.5,ymax=axlightymax+1,
                            **hspan_args)
            axdarkymax=axdark.get_ylim()[1] # The ymax for both axes should be the same, but do this just in case.
            axdark.axhspan(xmin=0.5,
                           ymin=axdarkymax+0.5,ymax=axdarkymax+1,
                           **hspan_args)
            ### Get summ histo plots ylim; compile and set all the ylims to same max limit.
            for l in summ_hist_cx_axes:
                maxy=list()
                maxy.append(l.get_ylim()[1])
            maxy=np.array(maxy)
            for l in summ_hist_cx_axes:
                l.set_ylim(0,maxy.max())
                l.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=4))
            # Add title.
            fig1.suptitle(status+' '+intensity+' '+osar.driver,
                         y=0.96,fontsize=22)
            out.append(fig1)
            # # Save figure.
            # fig1.savefig('cX_'+status+'_'+intensity+'_'+osar.driver+'.png',**savefig_args)


            ##### Create the plot for speed.
            summ_hist_speed_axes=list()
            fig2,((axlight,axdark,axcolbar),
                 (axlightsum,axdarksum,xtra))=plt.subplots(**mainplot_args)
            xtra.axis('off') # We have created an extra axes for spacing; hide it!

            # Plot heatmaps for border approaches from light.
            start_light_speed_all=start_light_speed_all.loc[:,-5.:5.]
            sns.heatmap(start_light_speed_all,
                        cbar=False,ax=axlight,**heatmap_speed_args)
            axlight.xaxis.set_visible(False) # Hide xaxis for aesthetics and clarity.
            # Set the title
            axlight.text(.5,1.04,
                         'From Light',
                         fontsize=18,
                         horizontalalignment='center',
                         transform=axlight.transAxes)
            # Plot the summary histogram for light approaches.
            mean=start_light_speed_all.mean()
            halfci=1.96*start_light_speed_all.sem()
            start_light_speed_all.mean().plot(color='black',linewidth=0.5,ax=axlightsum)
            axlightsum.fill_between(x=mean.index,
                                    y1=mean-halfci, y2=mean+halfci,
                                    alpha=0.25, color=light_color)
            axlightsum.set_xlim(start_light_speed_all.columns[0],
                               start_light_speed_all.columns[-1])
            sns.despine(ax=axlightsum,left=True)
            axlightsum.set_xlabel('Distance from\nlight-dark border (mm)')
            summ_hist_speed_axes.append(axlightsum) # Save axes into list for future ylim equalization.

            # Plot heatmaps for border approaches from dark.
            start_dark_speed_all=start_dark_speed_all.loc[:,-5.:5.]
            sns.heatmap(start_dark_speed_all,
                        cbar_ax=axcolbar,ax=axdark,**heatmap_speed_args)
            axcolbar.set_ylabel('Speed (mm/s)',labelpad=20,rotation=270) # Add title to our colorbar.
            axdark.xaxis.set_visible(False)
            axdark.text(.5,1.04,
                        'From Dark',
                        fontsize=18,
                        horizontalalignment='center',
                        transform=axdark.transAxes)
            # Plot the summary histogram for dark approaches.
            mean=start_dark_speed_all.mean()
            halfci=1.96*start_dark_speed_all.sem()
            start_dark_speed_all.mean().plot(color='black',linewidth=0.5,ax=axdarksum)
            axdarksum.fill_between(x=mean.index,y1=mean-halfci,y2=mean+halfci,alpha=0.25,color='grey')
            axdarksum.set_xlim(start_dark_speed_all.columns[0],
                               start_dark_speed_all.columns[-1])
            sns.despine(ax=axdarksum,left=True)
            axdarksum.set_xlabel('Distance from\nlight-dark border (mm)')
            summ_hist_speed_axes.append(axdarksum)
            # Aesthetics commmon to both axes.
            for a in [axlight, axdark]:
                # Draw reference line.
                a.axvline(x=np.where(edges[:-1]==0.)[0][0],
                          linewidth=0.55,
                          color='white')
            # Draw reference hbar to indicate lighted region.
            axlightymax=axlight.get_ylim()[1]
            axlight.axhspan(xmax=0.5,
                            ymin=axlightymax+0.5,ymax=axlightymax+1,
                            **hspan_args)
            axdarkymax=axdark.get_ylim()[1] # The ymax for both axes should be the same, but do this just in case.
            axdark.axhspan(xmin=0.5,
                           ymin=axdarkymax+0.5,ymax=axdarkymax+1,
                           **hspan_args)
            for l in summ_hist_speed_axes:
                maxy=list()
                maxy.append(l.get_ylim()[1])
            maxy=np.array(maxy)
            for l in summ_hist_speed_axes:
                l.set_ylim(0,maxy.max()+1) # Adding 1 to ylim max to give more buffer?
                l.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=4))
            # Add title.
            fig2.suptitle(status+' '+intensity+' '+osar.driver,
                         y=0.96,fontsize=22)
            out.append(fig2)
            # # Save figure.
            # fig2.savefig('speed_'+status+'_'+intensity+'_'+osar.driver+'.png',**savefig_args)

    return out
