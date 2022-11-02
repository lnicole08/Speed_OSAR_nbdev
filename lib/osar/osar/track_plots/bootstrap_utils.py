#!/usr/bin/python
# -*-coding: utf-8 -*-
# Author: Joses Ho
# Email : joseshowh@gmail.com



def _create_jackknife_indexes(data):
    """
    Given an array-like, creates a jackknife bootstrap.

    For a given set of data Y, the jackknife bootstrap sample J[i]
    is defined as the data set Y with the ith data point deleted.

    Keywords
    --------
    data: array-like

    Returns
    -------
    Generator that yields all jackknife bootstrap samples.
    """
    from numpy import arange, delete

    index_range = arange(0, len(data))
    return (delete(index_range, i) for i in index_range)



def _calc_accel(jack_dist):
    from numpy import mean as npmean
    from numpy import sum as npsum
    from numpy import errstate

    jack_mean = npmean(jack_dist)

    numer = npsum((jack_mean - jack_dist)**3)
    denom = 6.0 * (npsum((jack_mean - jack_dist)**2) ** 1.5)

    with errstate(invalid='ignore'):
        # does not raise warning if invalid division encountered.
        return numer / denom



def _compute_1group_jackknife(x, func, *args, **kwargs):
    """
    Returns the jackknife bootstraps for func(x).
    """
    jackknives = [i for i in _create_jackknife_indexes(x)]
    out = [func(x[j], *args, **kwargs) for j in jackknives]
    del jackknives # memory management.
    return out
    
    
    
def _create_bootstrap_indexes(array, resamples=5000, random_seed=12345):
    """Given an array-like, returns a generator of bootstrap indexes
    to be used for resampling.
    """
    import numpy as np
    
    np.random.seed(random_seed)
    
    indexes = range(0, len(array))
    
    # reset seed
    np.random.seed()

    return (np.random.choice(indexes, len(indexes), replace=True)
            for i in range(0, resamples))




def compute_1group_bootstraps(x, func, resamples=5000, random_seed=12345,
                             *args, **kwargs):
    """Bootstraps func(x), with the number of specified resamples."""

    # Create bootstrap indexes.
    bs_index_kwargs = {'resamples': resamples, 'random_seed': random_seed}
    boot_indexes = _create_bootstrap_indexes(x, **bs_index_kwargs)

    out = [func(x[b], *args, **kwargs) for b in boot_indexes]
    del boot_indexes
    return out



def compute_1group_bias_correction(x, bootstraps, func, *args, **kwargs):
    from scipy.stats import norm
    metric = func(x, *args, **kwargs)
    prop_boots_less_than_metric = sum(bootstraps < metric) / len(bootstraps)

    return norm.ppf(prop_boots_less_than_metric)



def _compute_alpha_from_ci(ci):
    if ci < 0 or ci > 100:
        raise ValueError("`ci` must be a number between 0 and 100.")

    return (100. - ci) / 100.



def _compute_quantile(z, bias, acceleration):
    numer = bias + z
    denom = 1 - (acceleration * numer)

    return bias + (numer / denom)
    
    
def compute_interval_limits(bias, acceleration, n_boots, ci=95):
    """
    Returns the indexes of the interval limits for a given bootstrap.

    Supply the bias, acceleration factor, and number of bootstraps.
    """
    from scipy.stats import norm
    from numpy import isnan, nan

    alpha = _compute_alpha_from_ci(ci)

    alpha_low = alpha / 2
    alpha_high = 1 - (alpha / 2)

    z_low = norm.ppf(alpha_low)
    z_high = norm.ppf(alpha_high)

    kws = {'bias': bias, 'acceleration': acceleration}
    low = _compute_quantile(z_low, **kws)
    high = _compute_quantile(z_high, **kws)

    if isnan(low) or isnan(high):
        return low, high

    else:
        low = int(norm.cdf(low) * n_boots)
        high = int(norm.cdf(high) * n_boots)
        return low, high
        
        
        
def summary_ci_1group(x, func, resamples=5000, alpha=0.05, random_seed=12345,
                      sort_bootstraps=True, *args, **kwargs):
    """
    Given an array-like x, returns func(x), and a bootstrap confidence
    interval of func(x).

    Keywords
    --------
    x: array-like
        An numerical iterable.

    func: function
        The function to be applied to x.

    resamples: int, default 5000
        The number of bootstrap resamples to be taken of func(x).

    alpha: float, default 0.05
        Denotes the likelihood that the confidence interval produced
        _does not_ include the true summary statistic. When alpha = 0.05,
        a 95% confidence interval is produced.

    random_seed: int, default 12345
        `random_seed` is used to seed the random number generator during
        bootstrap resampling. This ensures that the confidence intervals
        reported are replicable.

    sort_bootstraps: boolean, default True



    Returns
    -------
    A dictionary with the following five keys:
        'summary': float.
            The outcome of func(x).

        'func': function.
            The function applied to x.

        'bca_ci_low': float
        'bca_ci_high': float.
            The bias-corrected and accelerated confidence interval, for the
            given alpha.

        'bootstraps': array.
            The bootstraps used to generate the confidence interval.
            These will be sorted in ascending order if `sort_bootstraps`
            was True.

    """
    from numpy import sort as npsort

    boots = compute_1group_bootstraps(x, func)
    bias = compute_1group_bias_correction(x, boots, func)

    jk = _compute_1group_jackknife(x, func)
    accel = _calc_accel(jk)
    del jk

    ci_idx = compute_interval_limits(bias, accel, resamples, alpha)

    boots_sorted = npsort(boots)

    low = boots_sorted[ci_idx[0]]
    high = boots_sorted[ci_idx[1]]

    if sort_bootstraps:
        B = boots_sorted
    else:
        B = boots
    del boots
    del boots_sorted

    out = {'summary': func(x), 'func': func,
            'bca_ci_low': low, 'bca_ci_high': high,
            'bootstraps': B}

    del B
    return out
