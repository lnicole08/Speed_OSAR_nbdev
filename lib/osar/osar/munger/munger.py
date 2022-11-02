#!/usr/bin/python
# -*-coding: utf-8 -*-
# Author: Joses Ho
# Email : joseshowh@gmail.com

"""
Convenience functions for munging of countlogs.
"""
import numba

def create_borders(path_to_data,border_folder):
    '''This function creates a DataFrame with all the borders for each fly,
    for each patternState.'''

    import os as os
    import numpy as np
    import pandas as pd

    borders_dir = os.path.join(path_to_data,border_folder)
    bordercsv = [f for f in os.listdir(borders_dir) if '.csv' in f]
    border_list = list()

    for n, border_csv in enumerate(np.sort(bordercsv)):
        # filename is the pattern name.
        p_name = border_csv.split('.')[0]
        temp = pd.read_csv(os.path.join(borders_dir, border_csv))
        # only get the borders in pixels.
        temp = temp.filter(regex='pix').T.iloc[:, 1:len(temp)]

        temp['flyID'] = range(1, len(temp)+1)
        temp['flyID'] = 'fly' + temp['flyID'].astype(str)
        temp['patternState'] = np.repeat(p_name, len(temp))
        temp.index = [temp.flyID, temp.patternState]
        temp.drop(['flyID', 'patternState'], axis = 1, inplace = True)

        border_list.append(temp)
    # Turn list into DataFrame
    out = pd.concat(border_list)
    # Add appropriate columns.
    out.columns = 'midzone' + pd.Series(range(1, len(temp.columns)+1))\
                                .astype(str)
    return out


def join_cols(df,cols,sep='; '):
    """
    Convenience function to concatenate all the columns found in
    the list `cols`, with `sep` as the delimiter.

    Keywords
    --------

    df: a pandas DataFrame

    cols: a list of column names in `df`.

    sep: str, default '; '
        The delimiter used to seperate the concatenated columns.
    """
    try:
        base_col = df[ cols[0] ].astype(str).copy()
        if len(cols)>1: # if more than one column...
            col_list = cols[1:].copy()
            try:
                for j,col in enumerate(col_list):
                    out_col = base_col + sep + df[ col ].astype(str)
                    base_col = out_col
                return out_col
            except KeyError:
                    print('`{}` is not found in the feeds. Please check.'.format(col))
        else:
            # only one column in list; so do nothing.
            return base_col
    except KeyError:
        print('`{}` is not found in the feeds. Please check.'.format(cols[0]))


def cat_categorical_columns(df, group_by, compare_by):
   """
   Convenience function to concatenate categorical columns for
   contrast plotting purposes.
   """
   df_out - df.copy()
   # add compare_by to the group_by list.
   gby = group_by.copy()
   gby.append(compare_by)
   # create new categorical column.
   df_out['plot_groups']=join_cols(df_out, group_by)
   # Create another categorical column.
   df_out['plot_groups_with_contrast']=join_cols(df_out, gby)
   return df_out


@numba.jit(nopython=True)
def pi_for_against(votes_for, votes_against):
    """
    Given two integers `votes_for` and `votes_against`,
    computes the PI.
    """
    diff_votes = votes_for - votes_against
    total_votes = votes_for + votes_against

    return diff_votes / total_votes


@numba.jit(nopython=True)
def pi_array(array):
    """
    Calculate the preference index for an array of 1s and 0s, where 1 denotes
    a choice for the object of interest.
    """
    total_votes = len(array)
    votes_for = array.sum()
    votes_against = total_votes - votes_for
    votes_diff = votes_for - votes_against
    pi = votes_diff / total_votes

    return pi


# # @numba.jit(nopython=True)
# def log2_ratio(numer, denom):
#     import numpy as np
#
#     numer = numer.copy()
#     denom = denom.copy()
#     numer[numer == 0] = np.nan
#     denom[denom == 0] = np.nan
#     speed_ratio = np.divide(numer, denom)
#     speed_ratio[speed_ratio == 0] = np.nan
#
#     return np.log2(speed_ratio)



def is_on_other_side(track):
    """
    Convenience function to compute if the final position is on the same side.
    Used on trajectory output of `osar.border_plots`.
    """
    last_position = track.dropna().iloc[-1]
    return last_position > 0


def summarize_trajectory_exits(df):
    """
    Convenience function to summarize the trajectory exit decisions.
    `df` should be the sumamry output from osar.border_plots().
    """
    import pandas as pd

    fly_count = df.groupby(['driver','status','verge_type','flyid']).size()\
                  .groupby(['driver','status','verge_type']).size()
    trajectory_count = df.groupby(['driver','status','verge_type']).size()
    props = df.groupby(['driver','status','verge_type']).mean()
    out_summary = pd.concat([fly_count, trajectory_count, props],
                            axis=1)
    out_summary.columns = ['Fly_count', 'Trajectory_count',
                           'Proportion_{}'.format(props.columns[0])]
    out_summary.reset_index()

    return out_summary


def get_consecutive_runs(data):
    """
    Subroutine accessed by get_chunk_speeds().
    Numpy only function to get consecutive runs from an array-like.
    See https://stackoverflow.com/a/46606745/6202321
    """
    import numpy as np
    sequences = np.split(data,
                         np.array(
                            np.where(np.diff(data) > 1)[0]
                                 ) + 1
                        )
    return sequences


def get_chunk_times(indexed_epoch):
    out = []
    chunks = get_consecutive_runs(indexed_epoch.Seconds)

    for chunk in chunks:
        idx = chunk.index.tolist()
        chunk_time = indexed_epoch.loc[idx].seconds_diff.sum()
        out.append(chunk_time)

    return out


def get_chunk_distances(indexed_epoch, pixel_length_mm=0.104):
    out = []
    chunks = get_consecutive_runs(indexed_epoch.Seconds)

    for chunk in chunks:
        idx = chunk.index.tolist()
        chunk_dist = indexed_epoch.loc[idx].cX_smoothed_diff.abs().sum() * pixel_length_mm
        out.append(chunk_dist)

    return out


def get_chunk_speeds(indexed_epoch, pixel_length_mm=0.104):
    """Subroutine accessed by get_chunk_speeds_by_illumination().
    Used to compute the speed of each consecutive chunk in an indexed epoch,
    which is a non-contiguous set of consecutive stretches of simliar illumination."""

    out = []
    chunks = get_consecutive_runs(indexed_epoch.Seconds)

    for chunk in chunks:
        idx = chunk.index.tolist()
        chunk_time = indexed_epoch.loc[idx].seconds_diff.sum()
        if chunk_time < 1:
            pass

        else:
            chunk_dist = indexed_epoch.loc[idx].cX_smoothed_diff.abs().sum() * pixel_length_mm
            chunk_speed = chunk_dist / chunk_time
            out.append(chunk_speed)

    return out


def get_chunk_speeds_by_illumination(full_epoch, pixel_length_mm=0.104):
    """Subroutione accessed by log2_speedratio().
    Returns a dict for the speeds in each of the dark stretches and
    light stretches.
    """
    from numpy import log2

    in_light = full_epoch[full_epoch.lightStatus_new == 1].reset_index()
    in_dark = full_epoch[full_epoch.lightStatus_new == 0].reset_index()

    speed_chunks_in_light = get_chunk_speeds(in_light,
                                            pixel_length_mm=pixel_length_mm)
    speed_chunks_in_dark = get_chunk_speeds(in_dark,
                                            pixel_length_mm=pixel_length_mm)

    return {'speed_chunks_in_light':speed_chunks_in_light,
           'speed_chunks_in_dark':speed_chunks_in_dark}


def get_speedratio_from_speed_chunks(speed_chunk_dict):
    """Subroutione accessed by log2_speedratio().
    Takes the mean of the speeds in dark and light, and returns the log2
    of the light speed - dark speed ratio.
    """
    from numpy import mean, log2

    mean_speed_light = mean(speed_chunk_dict['speed_chunks_in_light'])
    mean_speed_dark = mean(speed_chunk_dict['speed_chunks_in_dark'])

    return log2(mean_speed_light/mean_speed_dark)


def log2_speedratio(full_epoch, pixel_length_mm=0.104):
    speeds = get_chunk_speeds_by_illumination(full_epoch,
                                              pixel_length_mm=pixel_length_mm)
    return get_speedratio_from_speed_chunks(speeds)




def process_countlogs(borders,
                     driver,
                     path_to_data,
                     countlog_folder,
                     TRAVERSAL_WINDOW_IDX,
                     LIGHT_ON_TIME_START,
                     LIGHT_OFF_TIME_START,
                     PIXEL_LENGTH_MM,
                     CHOICEZONE_WIDTH_MM,
                     WINDOW_SIZE,
                     ZONE_WIDTH ):

    import os
    import sys
    import pandas as pd
    import numpy as np

    # To store the fly metadata.
    LIGHT_INTENSITIES = ['Eighth','Quarter','Half','Full']

    # To store metadata for each fly.
    metadata = []
    # To store each fly's fly's smoothed track as its own dataframe.
    munged_flies = {}
    # To store the stats for each fly's smoothed track.
    results = []

    # bug catching.
    weird_flies_out = dict()

    # Define how zones map to light status,
    # based on the current pattern.
    # In the future, you can add more pattern variations below.
    pattern_01_map = {0.: 0, 1.: 1, 2.: 0, 3.: 1}
    pattern_10_map = {0.: 1, 1.: 0, 2.: 1, 3.: 0}
    # `BASELINE` is doubled for dark baseline indices.
    pattern_dict = {'BASELINE': pattern_01_map,
                    'Pattern 01': pattern_01_map,
                    'INTER-EPOCH': pattern_10_map,
                    'Pattern 10': pattern_10_map}
    zone_bin_labels = [0,1,2,3]
    choice_zone_labels = [-1, 1, -2, 2, -3, 3, -4]
    # This dict is used to remap the choiceZoneStatus values.
    # 10 = fly is in choice zone
    # 1 = entry into choice zone
    # -1 = exit from choice zone
    # np.nan = Fly is not in choice zone.
    cz_dict = {True: 10, False: np.nan, 11: 1, -11: -1}
    # This dict is used to map the decision codes to text.
    decision_map = {1:'to_light', -1:'to_dark',
                   11:'reverse_to_light', -11:'reverse_to_dark',
                   111:'.TO_LIGHT.', -111:'.TO_DARK.',
                   1111:'_TO_LIGHT_', -1111:'_TO_DARK_'}

    countlog_path = os.path.join(path_to_data, countlog_folder)
    csvlist = [c for c in os.listdir(countlog_path)
                if c.startswith('CountLog') and c.endswith('.csv')]

    # csvlist = [csvlist[5]] # comment out; only for DEBUG.

    for z, csvname in enumerate(csvlist):
        # print(csvname) ## REMOVE AFTER DEBUG
        sys.stdout.write('\rProcessing CSV {} of {}'.format(z+1, len(csvlist)))
        # Load the countlog.
        # csvname = csvlist[1] # comment out; only for DEBUG.
        # split names to add to each dataframe.
        name_split = csvname.strip('.csv').split('_')
        date = name_split[1]
        timestart = name_split[2]
        genotype = name_split[3]

        if genotype.startswith('w1118') or genotype.startswith('W1118'):
            status = 'Sibling'
        elif genotype.startswith(driver):
            status = 'Offspring'
        else:
            raise ValueError('Cannot identify genotype for file {}'.format(csvname))

        light_intensity = name_split[4]
        if light_intensity not in LIGHT_INTENSITIES:
            raise KeyError('{} is not one of the defined light intensities {}'\
                            .format(light_intensity,LIGHT_INTENSITIES))
        light_color = name_split[5]
        if light_color =='Green':
            opsin = 'GtACR1'
        elif light_color == 'Red':
            opsin = 'Chrimson'
        else:
            err = 'Cannot identify the opsin given the light color {}.'.format(light_color)
            err = err + 'Currently, only "Green" and "Red" are recognized.'
            raise ValueError(err)

        # Load the current CSV.
        countlog = pd.read_csv(os.path.join(
                                path_to_data, countlog_folder, csvname))
                                
        # Drop IDLE logrows.
        countlog = countlog[countlog.ExperimentState != 'Idle'].copy()
        
        # Pull out metadata columns, ensuring it is _not_ a slice but a copy.
        meta_cols = countlog[['Seconds','Time','ExperimentState','PatternState']].copy()
        
        # Properly re-assign PatternState based on ExperimentState.
        state_pattern_map = dict(zip(meta_cols.ExperimentState.unique(), 
                                     meta_cols.PatternState.unique()
                                     )
                                )
        meta_cols.loc[:,'PatternState'] = meta_cols.ExperimentState.map(state_pattern_map)

        ## Identify start of first epoch (Pattern 01).
        first_epoch_start = meta_cols[meta_cols.PatternState == 'Pattern 01'].Seconds.min()
        ## Set the entire 'Pattern 00' period before this as BASELINE.
        meta_cols.loc[(meta_cols.PatternState == 'Pattern 00') &
                      (meta_cols.Seconds < first_epoch_start), 'PatternState'] = 'BASELINE'
        ## The rest of the PatternStates with 'Pattern 00' will,
        ## by definition, be 'INTER-EPOCH'.
        meta_cols.loc[(meta_cols.PatternState == 'Pattern 00'), 'PatternState'] = 'INTER-EPOCH'

        # Get time diff to calculate velocity later.
        time_diff = countlog.Seconds.diff()
        meta_cols.loc[:,'seconds_diff'] = time_diff
        # Get unique objects in CSV for later.
        unique_objs = [o.strip('cX') for o in countlog.filter(regex='cX').columns]
        # Smooth the countlog.
        rolling_args = dict(window=WINDOW_SIZE,
                            min_periods=1,
                            center=True,
                            win_type="triang")

        cX_only = countlog.filter(regex='cX').copy()
        cX_smoothed = cX_only.rolling(**rolling_args).mean()
        # Calculate the raw diff, so our velocity calc will produce a signed number.
        cX_smoothed_diff = cX_smoothed.diff()
        cX_smoothed_diff.columns = cX_smoothed_diff.columns + '_diff'

        cX_velocity = cX_smoothed_diff.div(time_diff, axis='index')
        cX_velocity.columns = cX_smoothed.columns + '_velocity'

        # Determine direction status.
        # 1 = Rightwards; -1 = Leftwards
        cX_direction = np.sign(np.round(cX_velocity))
        cX_direction.columns = cX_smoothed.columns + '_direction'

        # Check if there is cY.
        cY_present = False
        cY_only = countlog.filter(regex='cY').copy()
        if len(cY_only.columns) == len(cX_only.columns):
            cY_present = True

        if cY_present is True:
            cY_smoothed = cY_only.rolling(**rolling_args).mean()
            cY_smoothed.columns = cY_only.columns + '_smoothed'

            cY_smoothed_diff = cY_smoothed.diff()
            cY_smoothed_diff.columns = cY_smoothed_diff.columns + '_diff'

            cY_velocity = cY_smoothed_diff.div(time_diff, axis='index')
            cY_velocity.columns = cY_smoothed.columns + '_velocity'

        all_zones = []
        all_light_status = []
        all_cz_occupancy = []
        # Identify the rightmost extreme of the arena.
        max_smoothed_cX = cX_smoothed.max().max()

        for zz, obj in enumerate(unique_objs):
            fly_zones = []
            fly_light_status = []
            fly_cz_occupancy = []
            fly_id = obj.strip('_').replace('Obj','fly')
            for p in borders.index.levels[1]:
                # Identify which zone-to-lightStatus dict to use.
                current_pattern_map = pattern_dict[p]
                # Identify the idx matching the current pattern.
                pattern_index = meta_cols[meta_cols.PatternState == p].index

                # Create the bins for re-zoning.
                zone_bins = [b for b in borders.loc[fly_id, p]]
                # Create the bins for identifying choice zone occupancies.
                choice_zone_bins = [[z-ZONE_WIDTH, z+ZONE_WIDTH]
                                    for z in zone_bins]
                choice_zone_bins = [j for k in choice_zone_bins for j in k]
                # Add 0 and rightmost extreme to as bin bookends.
                for b in [zone_bins, choice_zone_bins]:
                    b.insert(0, 0)
                    b.append(max_smoothed_cX)

                cX_for_pattern = cX_smoothed.loc[pattern_index, '{}cX'.format(obj)]
                try:
                    new_zones = pd.cut(cX_for_pattern,
                                       bins=zone_bins,
                                       labels=zone_bin_labels)
                except ValueError:
                    raise ValueError('Cannot create border bins.'
                                    ' Check that the framerate is correct'
                                    ' or that you have included BASELINE.csv'
                                    ' in the `border` subfolder.')
                new_light_status = new_zones.map(current_pattern_map)

                cz_occupancy = pd.cut(cX_for_pattern,
                                      bins=choice_zone_bins,
                                      labels=choice_zone_labels)
                cz_occupancy = cz_occupancy.map(lambda x: x > 0)

                # convert back to integer to mitigate this error:
                # https://github.com/pandas-dev/pandas/issues/18646
                # (which should be resolved in pandas 0.23.0)
                new_zones = new_zones.astype('float')

                fly_zones.append(new_zones)
                fly_light_status.append(new_light_status)
                fly_cz_occupancy.append(cz_occupancy)

            fly_zones = pd.concat(fly_zones, axis=0)
            all_zones.append(fly_zones)

            fly_light_status = pd.concat(fly_light_status, axis=0)
            all_light_status.append(fly_light_status)

            fly_cz_occupancy = pd.concat(fly_cz_occupancy, axis=0)
            all_cz_occupancy.append(fly_cz_occupancy)

        # concatenate all the lists into dataframes.
        # and modify their column names as appropriate.
        all_zones = pd.concat(all_zones, axis=1)
        all_zones.columns = all_zones.columns.str.replace('_cX', '_zone_new')
        all_light_status = pd.concat(all_light_status, axis=1)
        all_light_status.columns = all_light_status.columns.str.replace('_cX', '_lightStatus_new')
        all_cz_occupancy = pd.concat(all_cz_occupancy, axis=1)
        all_cz_occupancy.columns = all_cz_occupancy.columns.str.replace('_cX', '_choiceZoneStatus')
        # slightly convoluted way of organising the cz_occupancy dataframe...
        # We need to make sure we get a complete view of the experiment duration
        # to properly assign exits and entries.
        null_df = pd.DataFrame(index=cX_only.index,
                               columns=all_cz_occupancy.columns)
        null_df.fillna(False, inplace=True)
        null_df.loc[all_cz_occupancy.index] = all_cz_occupancy
        all_cz_occupancy = null_df

        # Identify all crossings of borders.
        all_crossings = all_light_status.diff().abs()
        all_crossings.columns = all_crossings.columns.str.replace('lightStatus_new', 'crossings')
        all_crossings.replace(0, np.nan, inplace=True)

        # Identify choice zone entries and exits.
        # diff and multiply by integer, so we can differentiate
        # True from ones.
        all_cz_entries = all_cz_occupancy.diff() * 11
        # Find turns, and assign when the turns are
        # to the diffed one.
        turns = all_cz_entries.loc[1:,]!=0
        all_cz_occupancy[turns] = all_cz_entries[turns]

        # Identify the epoch start and end indicies.
        epoch_start_idx=[]
        epoch_end_idx=[]
        for p in borders.index.levels[1]:
            pstate = meta_cols[meta_cols.PatternState == p]
            epoch_start_idx.append(pstate.index[0])
            epoch_end_idx.append(pstate.index[-1])
        epoch_start_idx = np.sort(epoch_start_idx)
        epoch_end_idx = np.sort(epoch_end_idx)

        # Handle boundary conditions where
        # the fly is still in the choice zone
        # when the LAST epoch ends (ie very end of experiment).
        still_in_cz = all_cz_occupancy.loc[epoch_end_idx[-1]]
        still_in_cz[still_in_cz == True] = -11

        # merge all the dataframes together.
        cX_smoothed.columns = cX_only.columns + '_smoothed'
        joined = cX_smoothed.join(all_cz_occupancy, how='outer')

        for df in [all_light_status, all_crossings, all_zones,]:
            joined = joined.join(df, how='outer')

        all_munged = pd.concat([cX_only, joined,
                               cX_velocity, cX_direction],
                              axis=1)
        if cY_present is True:
            cY_munged = pd.concat([cY_only, cY_smoothed, cY_velocity],
                                  axis=1)
            all_munged = pd.concat([all_munged, cY_munged], axis=1)

        # select the columns belonging to the fly.
        for j, obj in enumerate(unique_objs):
            # Create the fly_key used to uniquely ID each fly.
            fly_id = obj.strip('_').replace('Obj','fly')
            fly_key = '_'.join([driver,opsin,genotype,
               light_intensity,status,
               date,timestart,
               'border'+str(CHOICEZONE_WIDTH_MM)+'mm',
               fly_id])
            # Grab the light-dark borders for each object.
            mz = list()
            for p in borders.index.levels[1]:
                midzones = borders.loc[fly_id, p].tolist()
                mz.append(midzones)
            # Add this fly to metadata.
            meta = {'flyid': fly_key,
                    'exptDate': date,
                    'driver' : driver,
                    'original_count_log': csvname,
                    'expt_time_start': timestart,
                    'genotype': genotype,
                    'opsin': opsin,
                    'light_intensity': light_intensity,
                    'light_color': light_color,
                    'midzones' : mz,
                    'status': status}
            metadata.append(meta)

            temp_fly_df = all_munged.filter(regex=obj).copy()
            # rename the columns by removing the preceding `objX_`.
            temp_fly_df.columns = temp_fly_df.columns.str.replace(obj, '')
            # Add metadata columns to temp_fly_df
            temp_fly_df = pd.concat([meta_cols, temp_fly_df], axis=1)

            # Identify the epoch start and end indicies.
            # epoch_start_idx_ = []
            # epoch_end_idx_ = []
            epoch_idx = dict()
            for p in borders.index.levels[1]:
                pstate = temp_fly_df[temp_fly_df.PatternState == p]
                # epoch_start_idx_.append(pstate.index[0])
                # epoch_end_idx_.append(pstate.index[-1])
                epoch_idx[p] = (pstate.index[0], pstate.index[-1])
            # print(epoch_idx) # FOR DEBUG
            # epoch_start_idx_ = np.sort(epoch_start_idx_)
            # epoch_end_idx_ = np.sort(epoch_end_idx_)

            # remap the choiceZoneStatus.
            temp_fly_df.loc[:,'choiceZoneStatus'] = temp_fly_df.choiceZoneStatus.map(cz_dict)

            for pattern in epoch_idx.keys():
                start = epoch_idx[pattern][0]
                end = epoch_idx[pattern][1]

                # Ensure that any flies still already in the choice zone
                # at the start of each epoch are deemed to have entered.
                cz_status_at_start = temp_fly_df.loc[start, 'choiceZoneStatus']
                if cz_status_at_start == 10: # mark presence as entry.
                    temp_fly_df.loc[start, 'choiceZoneStatus'] = 1
                elif cz_status_at_start == -1: # Disregard exit at start.
                    temp_fly_df.loc[start, 'choiceZoneStatus'] = np.nan

                # Ensure that any flies still in the choice zone
                # either last minute entry or lingering till the
                # very last frame, are deemed to have exited.
                cz_status_at_end = temp_fly_df.loc[end, 'choiceZoneStatus']
                if cz_status_at_end == 1: # last minute entry.
                    temp_fly_df.loc[end, 'choiceZoneStatus'] = np.nan
                elif cz_status_at_end == 10: # lingering.
                    temp_fly_df.loc[end, 'choiceZoneStatus'] = -1

            entries = temp_fly_df[temp_fly_df.choiceZoneStatus == 1].index
            exits = temp_fly_df[temp_fly_df.choiceZoneStatus == -1].index

            # If there is one more exit than entrance,
            # drop the first exit.
            if len(entries) == len(exits) - 1 :
                temp_fly_df.loc[exits[0], 'choiceZoneStatus'] = np.nan
                exits = exits[1:]

            # Possibly deprecated error check. Leave here just in case?
            if len(entries) != len(exits):
                if len(entries) > len(exits):
                    err = '{} has MORE entries than exits.'.format(fly_key)
                elif len(entries) < len(exits):
                    err = '{} has LESS entries than exits.'.format(fly_key)
                weird_flies_out[fly_key] = err

            # Create column for decision code.
            temp_fly_df['decisionCode'] = np.repeat(np.nan, len(temp_fly_df))

            # identify the light status at entries and exits.
            light_at_entry = temp_fly_df.loc[entries, 'lightStatus_new'].values
            light_at_exit = temp_fly_df.loc[exits, 'lightStatus_new'].values
            # Double check that no NaNs are in the expected entries and exits.
            for ls_array in [light_at_entry, light_at_entry]:
                if np.isnan(ls_array).any(): raise ValueError('nan detected at light exit/entry.')
            direction_at_entry = temp_fly_df.loc[entries, 'cX_direction'].values
            direction_at_exit = temp_fly_df.loc[exits, 'cX_direction'].values

            # Perform boolean checks to identify the kind of choice made at each exit point.
            light_check = (light_at_entry == light_at_exit)
            velocity_check = (direction_at_entry == direction_at_exit)

            good_traversals = ~light_check & velocity_check
            traversals_to_light = good_traversals & (light_at_exit == 1)
            traversals_to_dark = good_traversals & (light_at_exit == 0)
            temp_fly_df.loc[exits[traversals_to_light], 'decisionCode'] = 1
            temp_fly_df.loc[exits[traversals_to_dark], 'decisionCode'] = -1

            good_reversals = light_check & ~velocity_check
            reversals_to_light = good_reversals & (light_at_exit == 1)
            reversals_to_dark = good_reversals & (light_at_exit == 0)
            temp_fly_df.loc[exits[reversals_to_light], 'decisionCode'] = 11
            temp_fly_df.loc[exits[reversals_to_dark], 'decisionCode'] = -11

            # special boundary conditions. Either the fly does not
            # make a traversal and moves deeper into the zone it is in,
            # or it stays in the choice zone and makes multiple traversals
            # before finally exiting with the same velocity.
            boundary1 = light_check & velocity_check
            boundary1_to_light = boundary1 & (light_at_exit == 1)
            boundary1_to_dark = boundary1 & (light_at_exit == 0)
            temp_fly_df.loc[exits[boundary1_to_light], 'decisionCode'] = 111
            temp_fly_df.loc[exits[boundary1_to_dark], 'decisionCode'] = -111

            # special boundary condition again. The fly was (ostensibly)
            # making a traversal but the light got turned on/off.
            boundary2 = ~light_check & ~velocity_check
            boundary2_to_light = boundary2 & (light_at_exit == 1)
            boundary2_to_dark = boundary2 & (light_at_exit == 0)
            temp_fly_df.loc[exits[boundary2_to_light], 'decisionCode'] = 1111
            temp_fly_df.loc[exits[boundary2_to_dark], 'decisionCode'] = -1111

            # Create a human-readable column of `decisionCode`
            temp_fly_df.loc[:,'decision'] = temp_fly_df.decisionCode.map(decision_map)

            munged_flies[fly_key] = temp_fly_df

            # Now, get the stats for each fly!
            # Create a dictionary to store various calculations.
            res = {}
            res['flyid'] = fly_key
            # Select the desired timewindow.
            winix = list()
            # Create cX_smoothed_diff
            temp_fly_df.loc[:, 'cX_smoothed_diff'] = temp_fly_df.cX_smoothed.diff()

            for k in epoch_idx.keys():
                idxs = epoch_idx[k]
                pattern_df = temp_fly_df.loc[idxs[0]: idxs[1]].copy()
                # Get window start and window end RELATIVE to the dataset's Seconds column.
                if k.startswith('Pattern '):
                    lightOnStart = pattern_df.Seconds.min()
                    wStart = lightOnStart + LIGHT_ON_TIME_START
                    wEnd = lightOnStart + LIGHT_OFF_TIME_START
                    pattern_df_win = pattern_df[(pattern_df.Seconds > wStart) &
                                                (pattern_df.Seconds < wEnd)]
                    # Get the INDEX of this window
                    ind = pattern_df_win.index
                    # convert to list, then array, then append to the list `winix`.
                    winix.append(np.array(ind.tolist()))

                elif k == 'BASELINE' or k == 'INTER-EPOCH':
                    pattern_df_win = pattern_df
                    # Get the INDEX of this window
                    ind = pattern_df.index
                else:
                    raise ValueError('{} is not a recognised pattern'.format(k))

                # Light Attraction
                res['reversals_to_light_{}'.format(k)] = float(len(pattern_df_win[pattern_df_win.decisionCode == 11]))
                res['traversals_to_light_{}'.format(k)] = float(len(pattern_df_win[pattern_df_win.decisionCode == 1]))
                res['reversals_to_dark_{}'.format(k)] = float(len(pattern_df_win[pattern_df_win.decisionCode == -11]))
                res['traversals_to_dark_{}'.format(k)] = float(len(pattern_df_win[pattern_df_win.decisionCode == -1]))

                votes_for_light = res['reversals_to_light_{}'.format(k)] + res['traversals_to_light_{}'.format(k)]
                votes_for_dark = res['reversals_to_dark_{}'.format(k)] + res['traversals_to_dark_{}'.format(k)]
                if (votes_for_light + votes_for_dark) > 0:
                    res['light_attraction_index_{}'.format(k)] = pi_for_against(votes_for_light, votes_for_dark)
                else:
                    res['light_attraction_index_{}'.format(k)] = np.nan


                # PI
                res['pi_smoothed_{}'.format(k)] = pi_array(pattern_df_win.lightStatus_new.values)


                # Speed
                res['log2_speed_ratio_{}'.format(k)] = log2_speedratio(pattern_df_win,
                                                                           PIXEL_LENGTH_MM)

                win_in_light = pattern_df_win[pattern_df_win.lightStatus_new == 1].reset_index()
                win_in_dark = pattern_df_win[pattern_df_win.lightStatus_new == 0].reset_index()

                res['distance_in_light_{}'.format(k)] = np.sum(get_chunk_distances(
                                                        win_in_light,
                                                        pixel_length_mm=PIXEL_LENGTH_MM))
                res['time_in_light_{}'.format(k)] = np.sum(get_chunk_times(win_in_light))

                res['distance_in_dark_{}'.format(k)] = np.sum(get_chunk_distances(
                                                        win_in_dark,
                                                        pixel_length_mm=PIXEL_LENGTH_MM))
                res['time_in_dark_{}'.format(k)] = np.sum(get_chunk_times(win_in_dark))


                # Time spent in choice window.
                inzone = pattern_df_win[pattern_df_win.choiceZoneStatus == 10]
                res['time_in_choice_zone_{}'.format(k)] = inzone.seconds_diff.sum()


                # Get indexes where traversals happened across the entire assay duration.
                xing_dark_idx = list()
                xing_light_idx = list()
                xing = temp_fly_df[temp_fly_df.crossings == 1].index
                for cross in xing:
                    try:
                        if temp_fly_df.loc[cross-1,'lightStatus_new'] > temp_fly_df.loc[cross+1,'lightStatus_new']:
                            # moving from light to dark.
                            xing_dark_idx.append(cross)
                        elif temp_fly_df.loc[cross-1,'lightStatus_new'] < temp_fly_df.loc[cross+1,'lightStatus_new']:
                            # moving from dark to light.
                            xing_light_idx.append(cross)
                    except KeyError:
                        pass
                if len(xing_dark_idx) > 0:
                    before = []
                    after = []
                    for idx in xing_dark_idx:
                        tempWinBefore = temp_fly_df.loc[idx-TRAVERSAL_WINDOW_IDX : idx-1]
                        before.append(np.divide(tempWinBefore.cX_smoothed_diff.abs().sum(),
                                                tempWinBefore.seconds_diff.abs().sum()))
                        # and speed after crossing...
                        tempWinAfter = temp_fly_df.loc[idx+1 : idx+TRAVERSAL_WINDOW_IDX]
                        after.append(np.divide(tempWinAfter.cX_smoothed_diff.abs().sum(),
                                               tempWinAfter.seconds_diff.abs().sum()))
                    res['speed_before_cross_to_dark_{}'.format(k)] = np.mean(before)
                    res['speed_after_cross_to_dark_{}'.format(k)] = np.mean(after)
                else:
                    res['speed_before_cross_to_dark_{}'.format(k)] = np.nan
                    res['speed_after_cross_to_dark_{}'.format(k)] = np.nan
                if len(xing_light_idx) > 0:
                    before = []
                    after = []
                    for idx in xing_light_idx:
                        tempWinBefore = temp_fly_df.loc[idx-TRAVERSAL_WINDOW_IDX : idx-1]
                        before.append(np.divide(tempWinBefore.cX_smoothed_diff.abs().sum(),
                                                tempWinBefore.seconds_diff.abs().sum()))
                        # and speed after crossing...
                        tempWinAfter = temp_fly_df.loc[idx+1 : idx+TRAVERSAL_WINDOW_IDX]
                        after.append(np.divide(tempWinAfter.cX_smoothed_diff.abs().sum(),
                                               tempWinAfter.seconds_diff.abs().sum()))
                    res['speed_before_cross_to_light_{}'.format(k)] = np.mean(before)
                    res['speed_after_cross_to_light_{}'.format(k)] = np.mean(after)
                else:
                    res['speed_before_cross_to_light_{}'.format(k)] = np.nan
                    res['speed_after_cross_to_light_{}'.format(k)] = np.nan

            # concatenate the list into a NumPy array.
            winix = np.concatenate(winix)
            # slice out the timewindow into its own DataFrame.
            fly_win = temp_fly_df.iloc[winix]

            # Compute log2 speed ratios for all illumination epochs.
            illum_epochs = [a for a in temp_fly_df.PatternState.unique()
                            if a != 'BASELINE' and a != 'INTER-EPOCH']
            temp_fly_df_reindexed = temp_fly_df.set_index('PatternState')
            all_illum_data = temp_fly_df_reindexed.loc[illum_epochs]

            all_speeds = get_chunk_speeds_by_illumination(all_illum_data, PIXEL_LENGTH_MM)
            res['speed_in_light'] = np.mean(all_speeds['speed_chunks_in_light'])
            res['speed_in_dark'] = np.mean(all_speeds['speed_chunks_in_dark'])
            res['log2_speed_ratio'] = np.log2(res['speed_in_light'] / res['speed_in_dark'])


            # ## The line below ensures we do not capture the BASELINE state.
            # # fly_win_light_on = fly_win[fly_win.ExperimentState != 'Light off']
            # # win_in_light = fly_win_light_on[fly_win_light_on.lightStatus_new == 1]
            # # win_in_dark = fly_win_light_on[fly_win_light_on.lightStatus_new == 0]
            # # win_in_light = fly_win[fly_win.lightStatus_new == 1]
            # # win_in_dark = fly_win[fly_win.lightStatus_new == 0]
            #
            # fly_in_dark = temp_fly_df[
            #                  (temp_fly_df.ExperimentState == 'Light off') |
            #                  ((temp_fly_df.ExperimentState != 'Light off') & \
            #                   (temp_fly_df.lightStatus_new == 0))]
            # fly_in_light = temp_fly_df[
            #                   (temp_fly_df.ExperimentState != 'Light off') & \
            #                   (temp_fly_df.lightStatus_new == 1)]
            #
            # with np.errstate(divide='ignore',invalid='ignore'):
            #     res['speed_in_light'] = np.divide(fly_in_light.cX_smoothed_diff.abs().sum(),
            #                                       fly_in_light.seconds_diff.sum())
            #     res['speed_in_dark'] = np.divide(fly_in_dark.cX_smoothed_diff.abs().sum(),
            #                                      fly_in_dark.seconds_diff.sum())

            res['total_distance_window_mm'] = fly_win.cX_smoothed_diff.abs().sum()
            res['total_distance_assay_mm'] = temp_fly_df.cX_smoothed_diff.abs().sum()

            results.append(res)

        sys.stdout.flush()

    print('\nSummarising results for all flies...')


    # convert metadata to DataFrame
    results = pd.DataFrame(results)

    # Speed Ratios.
    # for p in borders.index.levels[1]:
    #     results['log2_speed_ratio_{}'.format(p)] = log2_ratio(results['speed_in_light_{}'.format(p)].values,
    #                                                           results['speed_in_dark_{}'.format(p)].values)

    for _metric in ['pi_smoothed',
                    'light_attraction_index',
                    'log2_speed_ratio']:
        baseline = results['{}_BASELINE'.format(_metric)]
        patt01 = results['{}_{}'.format(_metric, 'Pattern 01')]
        patt10 = results['{}_{}'.format(_metric, 'Pattern 10')]

        corrected_patt01 = '{}_{}_baseline_corrected'.format(_metric, 'Pattern 01')
        results[corrected_patt01] = patt01 - baseline
        corrected_patt10 = '{}_{}_baseline_corrected'.format(_metric, 'Pattern 10')
        results[corrected_patt10] = patt10 - (-baseline)

        for c in [corrected_patt01, corrected_patt10]:
            results.loc[:, c] = np.divide(results[c], 2)


    results['light_attraction_index_baseline_corrected'] = results\
        .filter(regex='light_attraction_index_Pattern [0-9][0-9]_baseline_corrected')\
        .mean(axis=1)
    results['pi_smoothed_baseline_corrected'] = results\
        .filter(regex='pi_smoothed_Pattern [0-9][0-9]_baseline_corrected')\
        .mean(axis=1)

    results['light_attraction_index_not_baseline_corrected'] = results\
        .filter(regex='light_attraction_index_Pattern [0-9][0-9]$')\
        .mean(axis=1)
    results['pi_smoothed_not_baseline_corrected'] = results\
        .filter(regex='pi_smoothed_Pattern [0-9][0-9]$')\
        .mean(axis=1)

    # results['log2_speed_ratio'] = log2_ratio(results.speed_in_light.values,
    #                                          results.speed_in_dark.values)

    # Reversal and Traversal PI calculations.
    results['reversals_to_light'] = results.filter(regex='reversals_to_light_Pattern').sum(axis=1)
    results['traversals_to_light'] = results.filter(regex='traversals_to_light_Pattern').sum(axis=1)
    results['reversals_to_dark'] = results.filter(regex='reversals_to_dark_Pattern').sum(axis=1)
    results['traversals_to_dark'] = results.filter(regex='traversals_to_dark_Pattern').sum(axis=1)

    results['reversal_light_pi'] = pi_for_against(results.reversals_to_light.values,
                                                  results.reversals_to_dark.values)
    results['traversal_light_pi'] = pi_for_against(results.traversals_to_light.values,
                                                 results.traversals_to_dark.values)



    # results['log2_speed_ratio_traversal_dark'] = log2_ratio(results.speed_after_cross_to_dark.values,
    #                 results.speed_before_cross_to_dark.values)
    # results['log2_speed_ratio_traversal_light'] = log2_ratio(results.speed_after_cross_to_light.values,
    #                 results.speed_before_cross_to_light.values)


    # convert distance and speed columns from pix and pix/sec to mm and mm/sec
    cols = np.concatenate([results.filter(regex='speed').columns.tolist(),
                           results.filter(regex='distance').columns.tolist()])
    #commented out as of 12 nov 2021
    #results.loc[:,cols] = results.loc[:,cols]
    

    # Convert all the metadata into a DataFrame.
    metadata = pd.DataFrame(metadata)
    # Make relevant cols in metadata Categorical.
    for c in ['driver', 'opsin', 'status',
              'light_intensity','genotype']:
        metadata.loc[:,c]=metadata[c].astype('category')
    # Set category order.
    metadata.loc[:,'status'] = metadata.status.cat.set_categories(['Sibling',
                                                                'Offspring'], ordered = True)
    metadata.loc[:,'light_intensity'] = metadata.light_intensity.cat.set_categories(LIGHT_INTENSITIES,
                                                                                    ordered = True)
    metadata.loc[:,'light_intensity'] = metadata.light_intensity.cat.remove_unused_categories()
    metadata.sort_values(by=['status', 'light_intensity',
                             'exptDate','expt_time_start'],
                         inplace=True)
    metadata.reset_index(drop=True, inplace=True)

    # based on the ordering of 'status', set 'genotypes' as categorical.
    genotypes = [a for a in metadata.genotype.unique()]
    metadata.loc[:,'genotype']=metadata.genotype.cat.set_categories(genotypes,
                                                                    ordered=True)

    # join metadata with results!
    results = metadata.set_index('flyid').join(results.set_index('flyid'))
    results.drop(labels=['midzones', 'original_count_log'], axis=1, inplace=True)


    print('All done.')
    return metadata.set_index('flyid'), munged_flies, results, weird_flies_out



def produce_delta_metrics(osar, driver=None,
                          light_intensities=['Half', 'Full'],
                          distance_cutoff=0,
                          include_peroi=True,
                          epoch='Pattern 01'):
    """Convenience function for producing metrics from an osar.results
    pandas DataFrame.

    osar: OSAR object.

    light_intensities: list of light intensities to include.

    epoch: the epoch to produce delta metrics for.
    """
    import pandas as pd
    import dabest
    from ..track_plots.track_plots import border_plots

    # Get desired results.
    df = osar.results[osar.results.light_intensity.isin(light_intensities)]
    df_indexed = df.set_index(['light_intensity','status'])

    # Compute delta PEROI.
    if include_peroi:
        t, summ = border_plots(osar, epoch,
                                light_intensities=light_intensities,
                                distance_cutoff=distance_cutoff,
                                show_plot=False, verbose=False)
        VERGE_TYPES = summ.verge_type.unique()
        summ_indexed = summ.set_index(['light_intensity',
                                        'verge_type', 'status'])

    if driver is None:
        DRIVER = osar.driver
    else:
        DRIVER = driver
    OPSIN = osar.results.opsin[0]

    out = []

    for light_intensity in light_intensities:
        results = {}
        results['light_intensity'] = light_intensity
        results['driver'] = DRIVER
        results['opsin'] = OPSIN

        sibling = df_indexed.loc[light_intensity].loc['Sibling']
        offspring = df_indexed.loc[light_intensity].loc['Offspring']

        delta_pi = dabest.bootstrap_tools.bootstrap(sibling.loc[:,'pi_smoothed_{}'.format(epoch)],
                                                    offspring.loc[:,'pi_smoothed_{}'.format(epoch)])

        results['delta_pi'] = delta_pi.summary
        results['delta_pi_ci_low'] = delta_pi.bca_ci_low
        results['delta_pi_ci_high'] = delta_pi.bca_ci_high


        delta_lai = dabest.bootstrap_tools.bootstrap(sibling.loc[:,'light_attraction_index_{}'.format(epoch)],
                                                     offspring.loc[:,'light_attraction_index_{}'.format(epoch)])
        results['delta_lai'] = delta_lai.summary
        results['delta_lai_ci_low'] = delta_lai.bca_ci_low
        results['delta_lai_ci_high'] = delta_lai.bca_ci_high


        delta_speedratio = dabest.bootstrap_tools.bootstrap(sibling.loc[:,'log2_speed_ratio_{}'.format(epoch)],
                                                            offspring.loc[:,'log2_speed_ratio_{}'.format(epoch)])
        results['delta_speedratio'] = delta_speedratio.summary
        results['delta_speedratio_ci_low'] = delta_speedratio.bca_ci_low
        results['delta_speedratio_ci_high'] = delta_speedratio.bca_ci_high

        if include_peroi:
            for verge in VERGE_TYPES:
                summ_current_slice = summ_indexed.loc[light_intensity].loc[verge]

                sibling = summ_current_slice.loc['Sibling'].groupby('flyid').mean()
                offspring = summ_current_slice.loc['Offspring'].groupby('flyid').mean()

                delta_peroi = dabest.bootstrap_tools.bootstrap(sibling.exit_to_other_side,
                                                               offspring.exit_to_other_side)
                MEASURE = 'delta_peroi_{}'.format(verge)
                results[MEASURE] = delta_peroi.summary
                results['{}_ci_low'.format(MEASURE)] = delta_peroi.bca_ci_low
                results['{}_ci_high'.format(MEASURE)] = delta_peroi.bca_ci_high

        out.append(results)

    return pd.DataFrame(out)
