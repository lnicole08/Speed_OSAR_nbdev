#!/usr/bin/python
# -*-coding: utf-8 -*-
# Author: Joses Ho
# Email : joseshowh@gmail.com

class timeseries_plotter(object):
    """
    Timeseries plot class for osar object.

    Plots available:
    ----------------
    • To produce contrast plots for PI, use the `pi` method.
    e.g `my_osar_object.timeseries.pi`

    Plots available:
        pi_and_speed

    """



    def __init__(self, osar): # pass along dataframe.
        # from ..plot_helpers import plot_helpers as __plot_helpers

        self.__meta = osar.flies
        self.__driver = osar.driver
        self.__tracks = osar.tracks
        self.__framerate = osar.FRAME_RATE
        self.__pixel_mm = osar.PIXEL_LENGTH_MM



    def pi_and_speed(self,
                     rolling_winsize_seconds=5,
                     resample_winsize_seconds=0.1,
                     # smoothing_winsize_seconds=5,
                     show_illumination_epochs=True,
                     return_raw_data=False,
                     verbose=False):
        """
        Create timeseries plot, with PI and Speed over time.

        Keywords
        --------
        rolling_winsize_seconds: int, default 10
            The time window used to calculate the PI. By default, for each time-
            point, the PI for the last 10 seconds will be computed, for each
            fly.

        resample_winsize_seconds: float, default 0.05
            The new time interval to resample each fly's timeseries by. This is
            to allow proper calculation of the mean.

        show_illumination_epochs: boolean, default True
            If true, plot the illumination epochs as rectangles at the top of
            the PI axes.

        return_raw_data: boolean, default False
            If True, will return a pandas DataFrame of the rolling PI,
            and a pandas DataFrame of the rolling x-velocity.

        Returns
        -------
        An array of matplotlib Axes objects, (and a pandas DataFrame of the
        rolling PI, and a pandas DataFrame of the rolling x-velocity, if
        return_raw_data is True).
        """
        from .ts_helper import time_series_generic, timeseries_plot

        out = timeseries_plot(self.__driver,
                             self.__meta.copy(),
                             self.__tracks,
                             self.__framerate,
                             self.__pixel_mm,
                             rolling_winsize_seconds=rolling_winsize_seconds,
                             resample_winsize_seconds=resample_winsize_seconds,
                             show_illumination_epochs=show_illumination_epochs,
                             return_raw_data=return_raw_data,
                             verbose=verbose)

        return out
