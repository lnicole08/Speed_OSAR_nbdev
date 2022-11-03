#!/usr/bin/python
# -*-coding: utf-8 -*-
# Author: Joses Ho
# Email : joseshowh@gmail.com


def bootstrap_deltas(osar_results, epoch):
    """Convenience function to extract deltas from the results of an OSAR
    experiment."""

    import pandas as pd
    import dabest

    DRIVER = osar_results.driver[0]
    OPSIN = osar_results.opsin[0]
    df = osar_results.copy()
    df_indexed = df.set_index(['light_intensity','status'])

    out = []

    for light_intensity in df.light_intensity.unique():
        r = {}
        r['light_intensity'] = light_intensity
        r['driver'] = DRIVER
        r['opsin'] = OPSIN


        sibling = df_indexed.loc[light_intensity].loc['Sibling']
        offspring = df_indexed.loc[light_intensity].loc['Offspring']


        pi_metric = 'pi_smoothed_{}'.format(epoch)
        delta_pi = dabest.bootstrap_tools.bootstrap(sibling.loc[:, pi_metric],
                                                    offspring.loc[:, pi_metric])

        r['delta_pi'] = delta_pi.summary
        r['delta_pi_ci_low'] = delta_pi.bca_ci_low
        r['delta_pi_ci_high'] = delta_pi.bca_ci_high


        lai_metric = 'light_attraction_index_{}'.format(epoch)
        delta_lai = dabest.bootstrap_tools.bootstrap(sibling.loc[:, lai_metric],
                                                     offspring.loc[:, lai_metric])
        r['delta_lai'] = delta_lai.summary
        r['delta_lai_ci_low'] = delta_lai.bca_ci_low
        r['delta_lai_ci_high'] = delta_lai.bca_ci_high


        speed_metric = 'log2_speed_ratio_{}'.format(epoch)
        delta_speedratio = dabest.bootstrap_tools.bootstrap(sibling.loc[:, speed_metric],
                                                            offspring.loc[:, speed_metric])
        r['delta_speedratio'] = delta_speedratio.summary
        r['delta_speedratio_ci_low'] = delta_speedratio.bca_ci_low
        r['delta_speedratio_ci_high'] = delta_speedratio.bca_ci_high

        out.append(r)

    return pd.DataFrame(out)
