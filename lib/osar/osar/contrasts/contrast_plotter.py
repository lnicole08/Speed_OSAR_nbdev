#!/usr/bin/python
# -*-coding: utf-8 -*-
# Author: Joses Ho
# Email : joseshowh@gmail.com

"""
Plot functions for osar objects.
"""

#      # #####  #####    ##   #####  #   #    # #    # #####   ####  #####  #####
#      # #    # #    #  #  #  #    #  # #     # ##  ## #    # #    # #    #   #
#      # #####  #    # #    # #    #   #      # # ## # #    # #    # #    #   #
#      # #    # #####  ###### #####    #      # #    # #####  #    # #####    #
#      # #    # #   #  #    # #   #    #      # #    # #      #    # #   #    #
###### # #####  #    # #    # #    #   #      # #    # #       ####  #    #   #

# Add submodules below. The respective .py scripts
# should be in the same folder as espresso_plotter.py.
from . import pi as _pi
from . import log2_speed_ratio as _log2_speed_ratio

class contrast_plotter(object):
    """
    Contrast class for osar object.

    Plots available:
    ----------------
    • To produce contrast plots for PI, use the `pi` method.
    e.g `my_osar_object.contrasts.pi`

    Plots available:
        time_spent_in_light
        light_attraction_index
        traversals_towards_light
        reversals_towards_light

    • To produce contrast plots for speed ratios, use the `log2_speed_ratio`
    method, e.g `my_osar_object.contrasts.log2_speed_ratio`

    Plots available:
        light_speed_against_dark_speed
        cross_to_dark --> deprecated for now
        cross_to_light --> deprecated for now
    """

    #    #    #    #    #####
    #    ##   #    #      #
    #    # #  #    #      #
    #    #  # #    #      #
    #    #   ##    #      #
    #    #    #    #      #

    def __init__(self,osar): # pass along dataframe.
        self._results = osar.results.copy()
        self._driver = osar.driver
        # call obj.contrasts.xxx to access these methods.
        self.pi = _pi.pi_plotter(self)
        self.log2_speed_ratio = _log2_speed_ratio.log2_speed_ratio_plotter(self)
