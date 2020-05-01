import pdb
import os
import math
from collections import deque
from bisect import insort, bisect_left
from itertools import islice

import numpy as np
from scipy import signal
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt

from . import z3dio
from . import timeio
from . import z3d_directory
from . import util

def interp_extrap(x,xp,fp):
    '''
    Wrapper for np.interp, that extrapolates beyond limits of xp
    Also ensures that xp and fp are sorted so that xp is increasing
    See np.info('interp') for more about x, xp, fp
    '''
    # first, sort xp and fp
    xi = np.argsort(xp)
    xps = xp[xi]
    fps = fp[xi]
    # find all x outside limits of xps
    lefti = x < xps[0]
    righti = x > xps[-1]
    # interpolate
    f = np.interp(x,xps,fps)
    # handle extrapolations
    # left
    m_left = (fps[1]-fps[0])/(xps[1]-xps[0])
    f[lefti] = fps[0] + (x[lefti]-xps[0])*m_left
    # right
    m_right = (fps[-1]-fps[-2])/(xps[-1]-xps[-2])
    f[righti] = fps[-2] + (x[righti]-xps[-2])*m_right
    return f


def interp_nans(x):
    '''
    interpolate over nans
    written in haste...
    '''
    x_nans = np.isnan(x)
    ttt = lambda z: z.nonzero()[0]
    x_interp = x.copy()
    x_interp[x_nans] = interp_extrap(ttt(x_nans),ttt(~x_nans),x[~x_nans])
    return x_interp


def trim_clips(txts,rxts,tx_clip=39.8,rx_clip=1.99):
    '''
    Crop out parts of txts and rxts where either clip
    By default, tx clips at +/- 40 and rx clips at +/- 2
    '''
    assert len(txts)==len(rxts), 'txts and rxts must be of the same length'



def unwrap_phase(phases,freqs=None,units='mrad'):
    '''
    phases, in units
    freqs, in Hz or log(Hz)
    Remove kinks in phase vs frequency plot due to -pi to pi wrap
    Assume that plot wraps no more than once between points
    '''
    scale = util.radian_conversion_factor(units)
    if freqs is None:
        freqs = np.arange(len(phases))
        idx_sort = np.arange(len(phases))
    else:
        # first, sort phases and freqs
        idx_sort = np.argsort(freqs)
    phases_sort = phases[idx_sort]*scale
    freqs_sort = freqs[idx_sort]
    # unwrap
    unwrapped_phases = np.unwrap(phases_sort)/scale
    # unsort
    idx_unsort = np.argsort(idx_sort)
    return unwrapped_phases[idx_unsort]


def unwrap_linear_phase(phases,freqs,units='mrad'):
    '''
    phases, in mrad
    freqs, in Hz
    Remove kinks in phase vs frequency plot due to -pi to pi
    Assume phase is due to time delay (linear wrt frequency: constant rad/Hz)
    Assume no aliasing, i.e., interpolate phase change linearly between first two frequencies
    Note: It seems that some phases are clipped to 0 to pi
    This function accounts for that, at the expense of being less robust
    '''
    scale = util.radian_conversion_factor(units)
    # first, sort phases and freqs
    idx_sort = np.argsort(freqs)
    phases_sort = phases[idx_sort]*scale
    freqs_sort = freqs[idx_sort]
    # Estimate delay from first point
    delay = phases_sort[0]/freqs_sort[0]
    # residual = (phases_sort - delay*freqs_sort)%(np.pi)
    residual = (phases_sort - delay*freqs_sort+np.pi/2.)%(np.pi)-np.pi/2.
    # Re-estimate delay from first quarter of points by least-squares
    head = int(len(phases_sort)/4)
    F = np.stack((freqs_sort[:head],np.ones(head)),axis=-1)
    [ddelay,b] = np.dot(np.linalg.inv(np.dot(F.T,F)),np.dot(F.T,residual[:head]))
    delay += ddelay
    # Unwrap
    residual = phases_sort - delay*freqs_sort - b
    phase_shift = (residual/(1.*np.pi)).round()*1.*np.pi
    unwrapped_phases = (phases_sort-phase_shift)/scale
    # unsort
    idx_unsort = np.argsort(idx_sort)
    return unwrapped_phases[idx_unsort]


def apply_board_cal(fcs,freqs,cal,board_number,sampling_rate):
    '''
    apply board calibration to Fourier coefficients fcs with frequencies freqs
    First axis of fcs and freqs should be same length
    '''
    assert len(fcs)==len(freqs), 'The first axis of fcs and freqs must be the same length!'
    if board_number[:2]=='0x':
        board_number = board_number[2:]
    cal_board = cal.query('CARDSNX==@board_number & ADFREQ==@sampling_rate')
    # cal_board.loc[:,'LOGFREQ'] = np.log2(cal_board.loc[:,'CALFREQ'])
    # interpolate phase in linear frequency - phase space
    cal_freqs = cal_board.CALFREQ.values
    unwrapped_cal_phase = unwrap_linear_phase(cal_board.PHASE.values,cal_freqs)
    unwrapped_interp_phase = interp_extrap(freqs,cal_freqs,unwrapped_cal_phase)/1000.
    interp_phase = unwrapped_interp_phase%(np.pi*2.) #unnecessary?
    # most linear in log frequency - log (a-magnitude) space
    # But, for stability, interpolate in linear frequency - magnitude space
    interp_mag = interp_extrap(freqs,cal_freqs,cal_board.MAG.values)
    # apply calibrations
    calibrations = interp_mag*np.exp(1j*interp_phase)
    if fcs.ndim <= 1:
        cal_fcs = fcs/calibrations
    elif fcs.ndim == 2:
        cal_fcs = fcs/calibrations[:,np.newaxis]
    else:
        newshape = [1]*fcs.ndim
        newshape[0] = fcs.shape[0]
        cal_fcs = fcs/calibrations.reshape(newshape)
    return cal_fcs


def apply_antenna_cal(fcs,freqs,antcal,coil_number):
    '''
    Apply antenna calibration to Fourier coefficients fcs with frequencies freqs
    First axis of fcs and freqs should be same length
    '''
    assert len(fcs)==len(freqs)
    antenna_cal = antcal.query('antenna_sn==@coil_number')
    # interpolate phase in log frequency space
    # antenna_cal.loc[:,'log_frequency'] = np.log2(antenna_cal.loc[:,'base_frequency'])
    log_antcal_freqs = np.log2(antenna_cal.base_frequency.values)
    log_freqs = np.log2(freqs)
    unwrapped_antcal_phase = unwrap_phase(antenna_cal.phz.values,log_antcal_freqs)
    unwrapped_interp_phase = interp_extrap(log_freqs,log_antcal_freqs,unwrapped_antcal_phase)/1000.
    interp_phase = unwrapped_interp_phase%(np.pi*2.) #unnecessary?
    # appears piecewise linear in log frequency - log magnitude space
    log_antcal_mag = np.log2(antenna_cal.mag.values)
    log_interp_mag = interp_extrap(log_freqs,log_antcal_freqs,log_antcal_mag)
    interp_mag = 2**log_interp_mag
    # apply calibrations
    calibrations = interp_mag*np.exp(1j*interp_phase)
    if fcs.ndim <= 1:
        cal_fcs = fcs/calibrations
    elif fcs.ndim == 2:
        cal_fcs = fcs/calibrations[:,np.newaxis]
    else:
        newshape = [1]*fcs.ndim
        newshape[0] = fcs.shape[0]
        cal_fcs = fcs/calibrations.reshape(newshape)
    return cal_fcs


def prewhiten(ts,keep_length=False):
    '''
    Whiten, after Myers et al 2011
    '''
    if keep_length:
        return np.insert(ts[1:]-ts[:-1],0,0)
    else:
        return ts[1:]-ts[:-1]


def postdarken(fc,freq,sampling_rate=4096):
    '''
    Darken, after Myers et al 2011
    '''
    return fc/(np.exp(2*np.pi*1j*freq/sampling_rate)-1)


def reshape_ts(ts,samples_per_cycle,i_shift=None,keep_length=False):
    '''
    Reshape time series for stacking, etc
    If keep_length, then pad with zeros and mask first and last rows as needed
    Otherwise, truncate samples to make ts/samples_per_cycle be an integer
    '''
    num_samples = len(ts)
    if i_shift is None:
        # choose i_shift to truncate samples at beginning of ts
        num_cycles = num_samples//samples_per_cycle
        i_shift = num_samples - num_cycles*samples_per_cycle
    else:
        i_shift = i_shift % samples_per_cycle
        num_cycles = (num_samples-i_shift)//samples_per_cycle
    if keep_length:
        front = (samples_per_cycle - i_shift) % samples_per_cycle
        back = (samples_per_cycle-(num_samples-num_cycles*samples_per_cycle-i_shift)) % samples_per_cycle
        padded = np.ma.append(np.ma.concatenate(([0]*front,ts)),[0]*back)
        num_padded_samples = len(padded)
        masked = np.ma.masked_array(padded)
        masked[:front] = np.ma.masked
        masked[num_padded_samples-back:] = np.ma.masked
        num_padded_cycles = num_padded_samples//samples_per_cycle
        if num_padded_samples != num_padded_cycles*samples_per_cycle:
            print(front,back,num_samples,num_padded_samples,num_cycles,num_padded_cycles)
            raise ValueError("number of samples isn't working out right...")
        cycles = masked.reshape(num_padded_cycles,samples_per_cycle)
    else:
        cycles = ts[i_shift:(num_cycles)*samples_per_cycle+i_shift].reshape(num_cycles,samples_per_cycle)
    return cycles


def trim_outliers(a,proportion=0.05):
    '''
    Trim the values at both tails by setting them equal to trimmed max and min
    CAUTION: this function overwrites the input array a
    TODO: allow other axis besides 0
    '''
    nobs = a.shape[0]
    lowercut = int(proportion * nobs)
    uppercut = nobs - lowercut
    if (lowercut > uppercut):
        raise ValueError("Proportion too big.")
    ai = np.argpartition(a, (lowercut, uppercut - 1), axis=0)
    minvals = np.take_along_axis(a,np.expand_dims(ai[lowercut,:],axis=0),axis=0)
    maxvals = np.take_along_axis(a,np.expand_dims(ai[uppercut-1,:],axis=0),axis=0)
    np.put_along_axis(a,ai[:lowercut,:],minvals,axis=0)
    np.put_along_axis(a,ai[uppercut:,:],maxvals,axis=0)


def trim_outliers_1D(a,proportion=0.05):
    '''
    Trim values at both tails
    Works for a masked 1D signal
    '''
    nobs = a.count()
    lowercut = int(proportion * nobs)
    uppercut = nobs - lowercut
    if (lowercut > uppercut):
        raise ValueError("Proportion too big.")
    [minval,maxval] = np.partition(a.compressed(), (lowercut, uppercut - 1))[[lowercut,uppercut-1]]
    a_out = np.ma.where(a<minval,minval,a)
    a_out = np.ma.where(a_out>maxval,maxval,a_out)
    return a_out


def mask_outliers(a,proportion=0.05):
    '''
    Mask values at both tails
    Works for a 1D signal
    '''
    nobs = a.count()
    lowercut = int(proportion * nobs)
    uppercut = nobs - lowercut
    if (lowercut > uppercut):
        raise ValueError("Proportion too big.")
    [minval,maxval] = np.partition(a.compressed(), (lowercut, uppercut - 1))[[lowercut,uppercut-1]]
    return np.ma.masked_outside(a,minval,maxval)


def running_mean(x, N, keep_length=False, pad='same'):
    '''
    compute running mean
    If keep_length, returns an array of the same size as x
    For pad='edge', for even N, bias mean backwards
    For pad='same', fill with first val
    '''

    cumsum = np.cumsum(np.insert(x, 0, 0, axis=0), axis=0) 
    if keep_length:
        xlen = len(x)
        front = N//2
        back = (N-1)//2
        rm = np.empty(np.shape(x))
        if pad=='edge':
            # first front samples
            front_counts = np.arange(N-front,N,dtype=float)
            back_counts = np.arange(N-1,N-1-back,-1,dtype=float)
            for dimension in range(cumsum.ndim-1):
                front_counts = front_counts[...,None]
                back_counts = back_counts[...,None]
            rm[:front] = cumsum[N-front:N] / counts
            # middle xlen - (N-1) samples
            rm[front:xlen-back] = (cumsum[N:] - cumsum[:-N]) / float(N)
            # last back samples
            rm[xlen-back:] = (cumsum[-1] - cumsum[-N:-N+back]) / back_counts
        if pad=='same':
            # middle xlen - (N-1) samples
            rm[front:xlen-back] = (cumsum[N:] - cumsum[:-N]) / float(N)
            # first front samples
            rm[:front] = rm[front]
            # last back samples
            rm[xlen-back:] = rm[xlen-back-1]
        return rm
    else:
        return (cumsum[N:] - cumsum[:-N]) / float(N)


def running_mean_fractional(x,n,keep_length=False):
    '''
    Compute running mean of x using a boxcar window of length n.
    Take into account fractional part of n.

    Parameters
    ----------
    x : array like
        The signal
    n : float
        length of window over which to compute running mean
    keep_length : boolean
        Whether to return an array that is the same length as x

    Returns
    -------
    rm : array
        The running mean over a boxcar window of length n
        If keep_length, then len(rm) = len(x)
        If keep_length=False:
            len(rm)=len(x)-2*ceil((n-1)/2)
        So, len(x)-len(rm) is always even when keep_length is False.

    '''

    cumsum = np.cumsum(np.insert(x, 0, 0, axis=0), axis=0) 
    c, k = math.modf((n-1)/2)
    k = int(k)
    if c == 0:
        # (n-1)/2 is an integer
        rm = running_mean(x,int(n),keep_length)
    else:
        if keep_length:
            nx = len(x)
            rm = np.empty(np.shape(x))
            front = int(np.ceil((n-1)/2))
            back = front
            # middle nx - front - back samples
            int_sum = cumsum[(2*k+2):-1]-cumsum[1:(-2*k-2)]
            frac_sum = c*(x[:(-2*k-2)]+x[(2*k+2):])
            rm[front:-back] = (int_sum + frac_sum)/float(n)
            # fill in first front samples and last back samples
            rm[:front] = rm[front]
            rm[-back:] = rm[-back-1]
        else:
            int_sum = cumsum[(2*k+2):-1]-cumsum[1:(-2*k-2)]
            frac_sum = c*(x[:(-2*k-2)]+x[(2*k+2):])
            rm = (int_sum + frac_sum)/float(n)
    return rm


def running_mean_windowed(x,window,mask=None):
    '''
    Compute running mean of x using a window
    Allow masking of x as well

    Parameters
    ----------
    x : array like
        The signal
    window : array like
        length of window over which to compute running mean
        must have same number of dimensions as x
    keep_length : boolean
        Whether to return an array that is the same length as x
    mask : array of floats
        mask to apply to x, where 1 is unmasked and 0 is masked
        other values are allowed
        must be same shape as x

    Returns
    -------
    rm : array
        The running mean over the given window
        Automatically normalizes the window
        If keep_length, then len(rm) = len(x)
        If keep_length=False:
            len(rm)=len(x)-2*ceil((n-1)/2)
        So, len(x)-len(rm) is always even when keep_length is False.

    '''

    num_dimensions = x.ndim
    assert num_dimensions == window.ndim, 'x and window must have the same number of dimensions'
    if mask is None:
        mask = np.ones(x.shape)

    assert x.shape==mask.shape, 'x and mask must have the same shape'

    if num_dimensions == 2:
        x_conv_raw = signal.convolve2d(x*mask,window,mode='same')
        mask_conv_raw = signal.convolve2d(mask,window,mode='same')
        conv = x_conv_raw/mask_conv_raw
    else:
        # raise ValueError('#s of dimensions other than 2 not yet implemented')
        x_conv_raw = signal.convolve(x*mask,window,mode='same')
        mask_conv_raw = signal.convolve(mask,window,mode='same')
        conv = x_conv_raw/mask_conv_raw

    return conv


def robust_running_mean(x, N, proportion=0.05):
    xlen = len(x)
    # print('Reshape and copy...')
    # x_cycles = np.ma.copy(reshape_ts(x,N,i_shift=0,keep_length=True))
    x_cycles = np.copy(reshape_ts(x,N,i_shift=0,keep_length=True))
    # print('Trim outliers...')
    # np.ma.apply_along_axis(trim_outliers_1D,0,x_cycles,proportion=proportion)
    trim_outliers(x_cycles,proportion=proportion)
    #cumsum = np.cumsum(np.insert(x_cycles.flatten(), 0, 0))
    #return (cumsum[N:] - cumsum[:-N]) / float(N)
    # print('Running mean...')
    return running_mean(x_cycles.flatten()[:xlen],N,keep_length=True)


def running_median_insort(seq, window_size):
    """Contributed by Peter Otten"""
    seq = iter(seq)
    d = deque()
    s = []
    result = []
    for item in islice(seq, window_size):
        d.append(item)
        insort(s, item)
        result.append(s[len(d)//2])
    m = window_size // 2
    for item in seq:
        old = d.popleft()
        d.append(item)
        del s[bisect_left(s, old)]
        insort(s, item)
        result.append(s[m])
    return result


# def running_median_waveform(x,N):
#     x_cycles = np.copy(reshape_ts(x,N))

# def mean_waveform(x,
#                   sampling_frequency,
#                   waveform_frequency,
#                   waveform_samples_per_waveform,
#                   waveform_offset
#                  ):


def demean_fractional(x,n,keep_length=False):
    '''
    remove running mean, as computed by running_mean_fractional
    '''
    rm = running_mean_fractional(x,n,keep_length=keep_length)
    if keep_length:
        dm = x-rm
    else:
        rm_truncation = int(np.ceil((n-1)/2))
        dm = x[rm_truncation:-rm_truncation] - rm
    return dm


def estimate_signal_samples_per_waveform(x,estimate=68,search_buffer=2,exclude_correlation=3):
    '''
    Estimate the number of samples in one waveform,
    using autocorrelation
    equivalent to estimating dominant frequency

    Parameters
    ----------
    x : array like
        The full signal
    estimate : int
        approximate number of samples, starting guess
    search_buffer : int
        number of samples in each direction to search at each iteration
        e.g. estimate=68, search_buffer=2:
        iteration 1: search 66 through 70
            choose 69
        iteration 2: search (136 through 140)/2
        etc.
    exclude_correlation : int
        number of multiples of estimate to exclude at the end of the autocorrelation
        necessary because the autocorrelations peaks are irregular near the end

    Returns
    -------
    estimate/factor : float
        The final estimated number of samples per waveform

    '''
    auto_correlation = np.correlate(x,x,mode='full')
    ac = auto_correlation[len(auto_correlation)//2:-exclude_correlation*estimate]
    factor = 2
    estimate = estimate * 2
    while (estimate+search_buffer) < len(ac):
        estimate = np.argmax(ac[estimate-search_buffer:estimate+search_buffer+1])+estimate-search_buffer
        factor *= 2
        estimate *= 2
    return estimate/factor


def estimate_frequency(x,estimate=60,search_radius=2):
    '''
    Estimate the closest dominant frequency to estimate,
    using Fourier analysis

    Parameters
    ----------
    x : array like
        The full signal
    estimate : float
        approximate dominant frequency, starting guess
    search_radius : float
        frequency range to search

    Returns
    -------
    frequency : float
        The estimated dominant frequency nearest estimate

    '''


def resample_to_fit_period(x,signal_samples_per_period=4096/60,
                  period_samples_per_period=64,signal_offset=0):
    '''
    Resample signal so that each period has an integer number of samples 
    Returns reshaped signal: each row is one period

    Parameters
    ----------
    x : 1D array_like
        The full signal.
    signal_samples_per_period : float, positive
        Length of one period, in signal samples.
    period_samples_per_period : int, positive
        Number of equally-spaced times to use to represent one period.
    signal_offset : optional number
        Sample number for first sample in signal x. This is useful for
        maintaining the period's phase across multiple pieces of a long
        signal.

    Returns
    -------
    resample : 2D array of floats
        The resampled signal, shaped so that each row is one period

    '''
    # useful constants
    num_signal_samples = len(x)
    period_samples_per_signal_sample = period_samples_per_period/signal_samples_per_period

    # determine period samples
    first_period_sample = int(np.ceil(signal_offset*period_samples_per_signal_sample))
    last_period_sample = int(np.floor((signal_offset+num_signal_samples-1)*period_samples_per_signal_sample))
    first_period_sample_time = first_period_sample/period_samples_per_signal_sample
    last_period_sample_time = last_period_sample/period_samples_per_signal_sample
    num_period_samples = last_period_sample-first_period_sample+1

    # create independent variable arrays in units of signal samples
    signal_sample_times = np.arange(num_signal_samples)+signal_offset
    period_sample_times = np.linspace(first_period_sample_time,last_period_sample_time,num=num_period_samples)
    period_interp = np.interp(period_sample_times,signal_sample_times,x)

    resample = reshape_ts(period_interp,
                          period_samples_per_period,
                          i_shift=first_period_sample,
                          keep_length=True)
    return resample


def mean_signal_per_period(x,signal_samples_per_period=4096/60,
                  period_samples_per_period=64,signal_offset=0):
    '''
    Estimate the average signal as it repeats over a certain period

    Parameters
    ----------
    x : array_like
        The full signal.
    signal_samples_per_period : float, positive
        Length of one period, in signal samples.
    period_samples_per_period : int, positive
        Number of equally-spaced times to use to represent one period.
    signal_offset : optional int
        Sample number for first sample in signal. This is useful for
        maintaining the period's phase across multiple pieces of a long
        signal.

    Returns
    -------
    mean_signal : ndarray of floats
        The estimated signal per period, with length of period_samples_per_period.

    '''

    resample = resample_to_fit_period(x,
                               signal_samples_per_period=signal_samples_per_period,
                               period_samples_per_period=period_samples_per_period,
                               signal_offset=signal_offset)
    return np.mean(resample, axis=0)


def running_mean_signal_per_period(x,signal_samples_per_period=4096/60,
                          period_samples_per_period=68,
                          periods_per_mean=15,signal_offset=0,mask=None):
    '''
    Estimate the average signal as it repeats over a certain period.
    Computes a running average via a boxcar window of a certain # of periods
    Allows for samples to be masked

    Parameters
    ----------
    x : array_like
        The full signal.
    signal_samples_per_period : float, positive
        Length of one period, in signal samples.
    period_samples_per_period : int, positive
        Number of equally-spaced times to use to represent one period.
    periods_per_mean: float, positive
        number of periods to average over in running mean
    signal_offset : optional int
        Sample number for first sample in signal. This is useful for
        maintaining the period's phase across multiple pieces of a long
        signal.
    mask: array of booleans of same length as signal
        Mask to be applied to signal if some samples should be discarded.
        True = masked, False = unmasked
        If None, no mask is applied.

    Returns
    -------
    running_mean_signal : ndarray of floats
        The estimated waveform, with length of period_samples_per_period.

    '''

    # useful constants
    num_signal_samples = len(x)
    period_samples_per_signal_sample = period_samples_per_period/signal_samples_per_period
    # convert mask to float array where 1 is unmasked and 0 is masked
    if mask is None:
        mask_float = np.ones(num_signal_samples)
    else:
        mask_float = (~np.array(mask)).astype(float)

    # resample to fit period and reshape
    resample = resample_to_fit_period(x, 
                                      signal_samples_per_period=signal_samples_per_period,
                                      period_samples_per_period=period_samples_per_period,
                                      signal_offset=signal_offset)
    # form resample mask based on mask: anywhere a masked sample was used for interpolation,
    # mask the resample
    resample_mask = resample_to_fit_period(mask_float, 
                                      signal_samples_per_period=signal_samples_per_period,
                                      period_samples_per_period=period_samples_per_period,
                                      signal_offset=signal_offset)
    # combine with np.ma mask from the reshaping
    resample_mask = np.floor(resample_mask)*((~resample.mask).astype(float))

    # compute running mean over columns
    # running_mean_masked = running_mean_fractional((resample*resample_mask),
    #                                               periods_per_mean,keep_length=True)
    # running_mean_mask = running_mean_fractional(resample_mask,periods_per_mean,
    #                                             keep_length=True)
    # running_mean_period = running_mean_masked/running_mean_mask

    hann_window = signal.hann(periods_per_mean+2)[:,np.newaxis]
    running_mean_period = running_mean_windowed(resample,hann_window,resample_mask)

    # pad by one period in each direction
    rm_shape = running_mean_period.shape
    rmp_pad_period = np.empty((rm_shape[0]+2,rm_shape[1]))
    rmp_pad_period[0,:] = running_mean_period[0,:]
    rmp_pad_period[1:-1,:] = running_mean_period
    rmp_pad_period[-1,:] = running_mean_period[-1,:]

    # resample back to signal sampling rate
    # Now, one "period" is the whole signal
    # a "signal sample" is a period sample
    # a "period sample" is a signal sample
    # offset is in period samples, and must include padding
    period_samples_per_signal = period_samples_per_signal_sample * num_signal_samples
    first_period_sample = int(np.ceil(signal_offset*period_samples_per_signal_sample))
    period_offset = first_period_sample - period_samples_per_period
    # print(period_samples_per_signal,num_signal_samples,period_offset)
    running_mean_signal = resample_to_fit_period(
        rmp_pad_period.flatten(), 
        signal_samples_per_period=period_samples_per_signal,
        period_samples_per_period=num_signal_samples,
        signal_offset=period_offset)

    # trim to original signal sample positions
    first_signal_sample = int(np.ceil(period_offset/period_samples_per_signal_sample))
    last_signal_sample = int(np.floor((period_offset+len(rmp_pad_period.flatten())-1)/period_samples_per_signal_sample))
    signal_start = signal_offset-first_signal_sample
    signal_end = signal_start+num_signal_samples
    trimmed_running_mean_signal = running_mean_signal.flatten()[signal_start:signal_end]
    return trimmed_running_mean_signal.data
    # return running_mean_signal.flatten()[signal_start:signal_end]
    # return (mask,resample.mask,resample_mask)


def unsample(resample,num_signal_samples,signal_samples_per_period=4096/60,
             period_samples_per_period=68,
             signal_offset=0,mask=None):
    # useful constants
    period_samples_per_signal_sample = period_samples_per_period/signal_samples_per_period
    # pad by one period in each direction
    rm_shape = resample.shape
    rmp_pad_period = np.empty((rm_shape[0]+2,rm_shape[1]))
    rmp_pad_period[0,:] = resample[0,:]
    rmp_pad_period[1:-1,:] = resample
    rmp_pad_period[-1,:] = resample[-1,:]

    # resample back to signal sampling rate
    # Now, one "period" is the whole signal
    # a "signal sample" is a period sample
    # a "period sample" is a signal sample
    # offset is in period samples, and must include padding
    period_samples_per_signal = period_samples_per_signal_sample * num_signal_samples
    first_period_sample = int(np.ceil(signal_offset*period_samples_per_signal_sample))
    period_offset = first_period_sample - period_samples_per_period
    unsampled_padded = resample_to_fit_period(
        rmp_pad_period.flatten(), 
        signal_samples_per_period=period_samples_per_signal,
        period_samples_per_period=num_signal_samples,
        signal_offset=period_offset)

    # trim to original signal sample positions
    first_signal_sample = int(np.ceil(period_offset/period_samples_per_signal_sample))
    last_signal_sample = int(np.floor((period_offset+len(rmp_pad_period.flatten())-1)/period_samples_per_signal_sample))
    signal_start = signal_offset-first_signal_sample
    signal_end = signal_start+num_signal_samples
    unsampled = unsampled_padded.flatten()[signal_start:signal_end]
    return unsampled.data


def waveform_average_filter(sig,sampling_rate=4096,filter_freq=60,
                            num_averaged_periods=15,num_period_samples=100,
                            mask=None):
    '''
    filter out a repeating waveform of known frequency by finding a running mean waveform
    '''
    # useful constants
    num_samples = len(sig)
    samples_per_noise_cycle = sampling_rate/filter_freq

    # de-mean at filter frequency
    demean = demean_fractional(sig,samples_per_noise_cycle,keep_length=True)
    # running mean waveform
    rm_sig = running_mean_signal_per_period(
        demean,
        signal_samples_per_period=samples_per_noise_cycle,
        period_samples_per_period=num_period_samples,
        periods_per_mean=num_averaged_periods,
        signal_offset=0,
        mask=mask)
    return sig-rm_sig.data


def subset_data(z3d,start,end,start_offset=0,end_offset=0,z3d_key='data'):
    '''
    Extract data from z3d during time interval start to end,
    offset by start_offset and end_offset seconds
    z3d_key specifies which dictionary key to access in the z3d.
    Note: this leaves out the record with timestamp end.
    '''
    if type(start)==pd.Timestamp:
        start += pd.Timedelta('1s')*start_offset
    else:
        raise ValueError('type of start expected to be pandas.Timestamp')
    if type(end)==pd.Timestamp:
        end += pd.Timedelta('1s')*end_offset
    else:
        raise ValueError('type of start expected to be pandas.Timestamp')
    sampling_rate = int(round(z3d['metadata']['A/D Rate']))
    mountain_times = timeio.gps_week_seconds_to_mountain(z3d['metadata']['GpsWeek'],z3d['gps_times'])
    # start_index = (mountain_times.index(start))*sampling_rate
    # end_index = (mountain_times.index(end))*sampling_rate
    start_index = np.argmax(np.array(mountain_times)>=start)
    end_index = np.argmax(np.array(mountain_times)>=end)
    # check for any True
    if start_index == 0:
        if mountain_times[0]<start:
            start_index = len(mountain_times)
    if end_index == 0:
        if mountain_times[0]<end:
            end_index = len(mountain_times)
    # data index depends on sampling rate, gps_times does not.
    if z3d_key == 'data':
        start_index *= sampling_rate
        end_index *= sampling_rate
    interval_data = z3d[z3d_key][start_index:end_index]
    return interval_data


def overlap_data(tx,rx,overlaps,i_pair):
    '''
    return tx and rx data that overlaps in time
    throw out first 2 secs and last sec, to avoid headaches
    '''
    tx_filename = tx.fullpath[overlaps.tx_ind.iloc[i_pair]]
    rx_filename = rx.fullpath[overlaps.rx_ind.iloc[i_pair]]

    tx_z3d=z3dio.read_z3d(tx_filename)
    rx_z3d=z3dio.read_z3d(rx_filename)

    '''
    sampling_rate = round(rx_z3d['metadata']['A/D Rate'])
    tx_mountain_times = timeio.gps_week_seconds_to_mountain(tx_z3d['metadata']['GpsWeek'],tx_z3d['gps_times'])
    rx_mountain_times = timeio.gps_week_seconds_to_mountain(rx_z3d['metadata']['GpsWeek'],rx_z3d['gps_times'])
    [tx_start,tx_end] = timeio.get_start_and_end_times_mountain(tx_z3d)
    [rx_start,rx_end] = timeio.get_start_and_end_times_mountain(rx_z3d)

    # rx_overlap_start_second = rx_mountain_times.index(overlaps.start[i_pair])
    rx_overlap_start_second = rx_mountain_times.index(overlaps.start[i_pair])+2
    rx_start_index = rx_overlap_start_second*sampling_rate
    # rx_overlap_end_second = rx_mountain_times.index(overlaps.end[i_pair])+1
    rx_overlap_end_second = rx_mountain_times.index(overlaps.end[i_pair])
    rx_end_index = rx_overlap_end_second*sampling_rate

    # tx_overlap_start_second = tx_mountain_times.index(overlaps.start[i_pair])
    tx_overlap_start_second = tx_mountain_times.index(overlaps.start[i_pair])+2
    tx_start_index = tx_overlap_start_second*sampling_rate
    # tx_overlap_end_second = tx_mountain_times.index(overlaps.end[i_pair])+1
    tx_overlap_end_second = tx_mountain_times.index(overlaps.end[i_pair])
    tx_end_index = tx_overlap_end_second*sampling_rate

    tx_overlap_data = tx_z3d['data'][tx_start_index:tx_end_index]
    rx_overlap_data = rx_z3d['data'][rx_start_index:rx_end_index]
    '''

    tx_overlap_data = subset_data(tx_z3d,overlaps.start.iloc[i_pair],overlaps.end.iloc[i_pair],start_offset=2)
    rx_overlap_data = subset_data(rx_z3d,overlaps.start.iloc[i_pair],overlaps.end.iloc[i_pair],start_offset=2)

    assert len(tx_overlap_data) == len(rx_overlap_data), 'tx and rx overlap data are of different lengths'
    return (tx_overlap_data,rx_overlap_data)


class overlap_reader:
    '''
    Class to temporarily store time series data during processing,
    to avoid reading the same file over and over
    '''
    def __init__(self,tx,rx,overlaps,read_all_tx = True,start_offset=2,end_offset=0):
        self.tx = tx
        self.rx = rx
        self.overlaps = overlaps
        self.tx_z3d = []
        self.rx_z3d = []
        self.tx_z3d_index = np.nan
        self.rx_z3d_index = np.nan
        self.read_all_tx = read_all_tx
        self.start_offset = start_offset
        self.end_offset = end_offset
        if read_all_tx:
            # tx_filenames = tx.fullpath
            print('Reading tx z3ds...')
            self.tx_z3d = [z3dio.read_z3d(tfn) for tfn in tx.fullpath]
            self.tx_z3d_index = tx.index.values
            print('Done!')

    def next_overlap_z3ds(self,i_pair):
        '''
        Return the z3d data and metadata for the next pair
        Used by next_overlap_data and next_overlap_gps_times
        '''
        tx_ind = self.overlaps.tx_ind.iloc[i_pair]
        tx_irow = self.tx.index.get_loc(tx_ind)
        rx_ind = self.overlaps.rx_ind.iloc[i_pair]
        # rx_irow = self.rx.index.get_loc(rx_ind)
        if self.read_all_tx:
            tx_z3d = self.tx_z3d[tx_irow]
        else:
            if self.tx_z3d_index == tx_ind:
                tx_z3d = self.tx_z3d
            else:
                print('Reading z3d data...')
                tx_filename = self.tx.fullpath[tx_ind]
                tx_z3d = z3dio.read_z3d(tx_filename)
                self.tx_z3d = tx_z3d
                self.tx_z3d_index = tx_ind
                print('Done!')
        if self.rx_z3d_index == rx_ind:
            rx_z3d = self.rx_z3d
        else:
            print('Reading z3d data...')
            rx_filename = self.rx.fullpath[rx_ind]
            rx_z3d = z3dio.read_z3d(rx_filename)
            self.rx_z3d = rx_z3d
            self.rx_z3d_index = rx_ind
            print('Done!')
        return(tx_z3d,rx_z3d)

    # def next_overlap_gps_times(self,i_pair):
    #     '''
    #     Return gps timestamps for tx and rx data, where they overlap
    #     '''
    #     tx_z3d,rx_z3d = self.next_overlap_z3ds(i_pair)
    #     tx_gps_times = subset_data(tx_z3d,
    #                                self.overlaps.start.iloc[i_pair],
    #                                self.overlaps.end.iloc[i_pair],
    #                                start_offset=self.start_offset,
    #                                end_offset=self.end_offset,
    #                                z3d_key='gps_times')
    #     rx_gps_times = subset_data(rx_z3d,
    #                                self.overlaps.start.iloc[i_pair],
    #                                self.overlaps.end.iloc[i_pair],
    #                                start_offset=self.start_offset,
    #                                end_offset=self.end_offset,
    #                                z3d_key='gps_times')
    #     assert len(tx_gps_times)==len(rx_gps_times), 'tx and rx gps times are of different lengths'
    #     return (tx_gps_times,rx_gps_times)

    def next_overlap_data(self,i_pair,z3d_key='data'):
        '''
        return tx and rx data that overlaps in time
        only read new z3ds if needed
        throw out first 2 secs and last sec, to avoid headaches
        i_pair is to access the ith row of overlaps, not the index 
        '''
        tx_z3d,rx_z3d = self.next_overlap_z3ds(i_pair)
        tx_overlap_data = subset_data(tx_z3d,
                                      self.overlaps.start.iloc[i_pair],
                                      self.overlaps.end.iloc[i_pair],
                                      start_offset=self.start_offset,
                                      end_offset=self.end_offset,
                                      z3d_key=z3d_key)
        rx_overlap_data = subset_data(rx_z3d,
                                      self.overlaps.start.iloc[i_pair],
                                      self.overlaps.end.iloc[i_pair],
                                      start_offset=self.start_offset,
                                      end_offset=self.end_offset,
                                      z3d_key=z3d_key)
        assert len(tx_overlap_data) == len(rx_overlap_data), 'tx and rx overlap data are of different lengths'
        # if len(tx_overlap_data) != len(rx_overlap_data):
            # print('Warning: tx and rx overlap data are of different lengths')
        return (tx_overlap_data,rx_overlap_data)


def on_time_offset_signed(sig,freq,sampling_rate=1):
    '''
    find positive on time, not just nearest polarity change
    '''
    num_samples = len(sig)
    samples_per_cycle = int(round(sampling_rate/freq))
    num_full_cycles = num_samples//samples_per_cycle
    ds = sig[1:]-sig[:-1]
    spike_train = np.zeros(samples_per_cycle)
    spike_train[0]=1.0
    spike_train[samples_per_cycle//2]=-1.0
    spike_train = np.tile(spike_train,num_full_cycles-1)
    conv = signal.convolve(ds,spike_train,mode='valid')
    #print(ds.shape,spike_train.shape,conv.shape)
    offset = np.argmax(conv)
    offset = offset % (samples_per_cycle)
    return offset


def on_time_offset(sig,freq,sampling_rate=1):
    num_samples = len(sig)
    samples_per_cycle = int(round(sampling_rate/freq))
    num_full_cycles = num_samples//samples_per_cycle
    ds = np.abs(sig[1:]-sig[:-1])
    spike_train = np.zeros(samples_per_cycle//2)
    spike_train[0]=1.0
    spike_train = np.tile(spike_train,num_full_cycles*2-1)
    conv = signal.convolve(ds,spike_train,mode='valid')
    #print(ds.shape,spike_train.shape,conv.shape)
    offset = np.argmax(conv)
    offset = offset % (samples_per_cycle//2)
    return offset


def mask_near_on_time(tx_signal,freq,sampling_rate=4096,pad_before=34,pad_after=100):
    '''
    create a mask for samples that occur near signal on-times
    '''
    num_samples = len(tx_signal)
    samples_per_tx_cycle = int(round(sampling_rate/freq))
    num_cycles = num_samples//samples_per_tx_cycle
    i_offset = on_time_offset(tx_signal,freq,sampling_rate)
    cycle_mask = np.array([True]*pad_after+[False]*(samples_per_tx_cycle//2-pad_before-pad_after)+[True]*pad_before)
    signal_mask = np.tile(cycle_mask,num_cycles*2+4)
    signal_mask = signal_mask[samples_per_tx_cycle-i_offset:num_samples+samples_per_tx_cycle-i_offset]
    return signal_mask



def DFC(sig,freq,sampling_rate=4096):
    '''
    Compute a fourier coefficient the old fashioned way
    '''
    num_samples = len(sig)
    return (sig*np.exp(-1j*2*np.pi*freq*np.arange(num_samples)/sampling_rate)).sum()


def IDFC(fc,freq,n_samples,n_scaling=None,sampling_rate=4096):
    '''
    Return the sinusoid corresponding to a lone fourier coefficient
    '''
    if n_scaling is None:
        n_scaling = n_samples
    return 2*fc*np.exp(1j*2*np.pi*freq*np.arange(n_samples)/sampling_rate)/n_scaling


def hat(N,mode='center'):
    if N % 2 == 0:
        hat = (np.arange(N)+0.5)/float(N)*2
    else:
        hat = (np.arange(N)+1)/float(N+1)*2
    if mode=='left':
        hat = hat[::-1]
        hat[:(N+1)//2] = 1.
    if mode=='right':
        hat[N//2:] = 1.
    if mode=='center':
        hat[N//2:]=hat[(N-1)//2::-1]
    return hat


def lockin(sig,txfreq,sampling_rate,pad_before=40,pad_after=100):
    '''
    Remove (60 Hz) noise in square waveform based on late times
    A notch filter with short bandwidth works okay, except for 4 Hz
    TODO: pass in i_offset; don't determine it from an rx signal!
    '''
    num_samples = len(sig)
    samples_per_cycle = int(round(sampling_rate/txfreq))
    i_offset = on_time_offset(sig,txfreq,sampling_rate=sampling_rate)
    num_cycles = num_samples//samples_per_cycle
    # cycles = reshape_ts(sig,samples_per_cycle,i_offset)
    #num_offset_cycles = (num_samples-i_offset)//samples_per_cycle
    #cycles = sig[i_offset:(num_offset_cycles)*samples_per_cycle+i_offset].reshape(num_offset_cycles,samples_per_cycle)
    cycle_mask = np.array([True]*pad_after+[False]*(samples_per_cycle//2-pad_before-pad_after)+[True]*pad_before)
    signal_mask = np.tile(cycle_mask,num_cycles*2+4)
    signal_mask = signal_mask[samples_per_cycle-i_offset:num_samples+samples_per_cycle-i_offset]
    masked = np.ma.masked_array(sig,signal_mask)
    running_mean_size = int(round(sampling_rate/60.))
    rm = running_mean(sig,running_mean_size,keep_length=True)
    # subtract running mean
    # masked_slice = slice(running_mean_size//2,num_samples-running_mean_size//2+1)
    # demeaned_masked = masked[masked_slice] - rm
    demeaned_masked = masked - rm
    demeaned_masked = mask_outliers(demeaned_masked,proportion=0.002)
    # compute fourier coefficients for windows
    fc60,n_unmasked=sliding_window(demeaned_masked,sampling_rate,txfreq,
                                   frequencies=60.,cycles_per_window=1,
                                   two_strides_per_window=True,window='boxcar')
    sinusoids = np.array([IDFC(fc,60.,samples_per_cycle,sampling_rate=sampling_rate) for fc in fc60])
    # correct for masked samples
    sinusoids = sinusoids*samples_per_cycle/n_unmasked[:,np.newaxis]
    hat_left = hat(samples_per_cycle,'left')
    hat_right = hat(samples_per_cycle,'right')
    hat_center = hat(samples_per_cycle,'center')
    tri_sin = sinusoids*hat_center
    tri_sin[0] = sinusoids[0]*hat_left
    tri_sin[-1] = sinusoids[-1]*hat_right
    full_sin = np.zeros((sinusoids.shape[0]+1)*samples_per_cycle//2)
    # I bet this next step could be sped up
    for i in np.arange(sinusoids.shape[0]):
        full_sin[i*samples_per_cycle//2:i*samples_per_cycle//2+samples_per_cycle] += np.real(tri_sin[i])
    # filt_data = sig[masked_slice]
    filt_data = sig
    filt_data = filt_data[:len(full_sin)] - full_sin
    return filt_data


def mean_cov(tf):
    '''
    Compute mean and covariance of transfer functions
    tf: complex transfer functions
    '''
    num_freq = tf.shape[0]
    mu = np.mean(tf,axis=1)
    cov = np.array([np.cov([tf.real[i_freq,:],tf.imag[i_freq,:]]) 
                    for i_freq in np.arange(num_freq)])
    return (mu,cov)


def LS(fct,fcr):
    '''
    Compute estimate and variance of transfer function via least-squares
    fct: transmitter Fourier coefficients
    fcr: receiver Fourier coefficients
    '''
    N = fct.shape[1]
    V = 1/np.sum(fct.real**2+fct.imag**2,axis=1)
    tfh = V*np.sum(np.conjugate(fct)*fcr,axis=1)
    # variance is confusing: what to divide by? N? N-1? N-2? 2N? 2N-1? 2N-2? 2N-4?
    res = fcr - fct*tfh[:,np.newaxis]
    var = np.sum(res.real**2+res.imag**2,axis=1)*V/(2*N-2)
    return (tfh,var)


def LS_cov(fct,fcr):
    '''
    Compute estimate and covariance of transfer function via least-squares
    fct: transmitter Fourier coefficients
    fcr: receiver Fourier coefficients
    '''
    num_freq = fct.shape[0]
    N = fct.shape[1]
    assert fct.shape == fcr.shape, 'fct and fcr must have the same shape'
    V = 1/np.sum(fct.real**2+fct.imag**2,axis=1)
    tfh = V*np.sum(np.conjugate(fct)*fcr,axis=1)
    # covariance is confusing: what to divide by? N? N-1? N-2? 2N? 2N-1? 2N-2? 2N-4?
    res = fcr - fct*tfh[:,np.newaxis]
    # covariance is confusing: no idea how to express in mat-vec form
    cov = []
    for i_freq in np.arange(num_freq):
        # cov_rx = np.cov([np.real(res[i_freq,:]),np.imag(res[i_freq,:])
        #                  for i_freq in np.arange(num_freq)])
        cov_rx = np.cov(np.real(res[i_freq,:]),np.imag(res[i_freq,:]),ddof=1)
        t = fct[i_freq,:]
        tT_cov_t_00 = np.sum(cov_rx[0,0]*t.real**2 + 
                             2*cov_rx[0,1]*t.real*t.imag + 
                             cov_rx[1,1]*t.imag**2)
        tT_cov_t_01 = np.sum(cov_rx[0,1]*(t.real**2-t.imag**2) + 
                             (cov_rx[1,1] - cov_rx[0,0])*t.real*t.imag)
        tT_cov_t_11 = np.sum(cov_rx[0,0]*t.imag**2 -
                             2*cov_rx[0,1]*t.real*t.imag +
                             cov_rx[1,1]*t.real**2)
        tT_cov_t = np.array([[tT_cov_t_00, tT_cov_t_01],
                             [tT_cov_t_01, tT_cov_t_11]])
        cov.append(tT_cov_t*V*V)
    # var = np.sum(res.real**2+res.imag**2,axis=1)*V/(2*N-2)
    return (tfh,cov)


def huber_norm(r,r0=1.5):
    '''
    Huber norm of array r
    = r**2/2 for abs(r)<r0
    = r0*abs(r) - r0**2/2 for abs(r)>=r0
    '''
    rabs = abs(r)
    return np.sum([ri**2/2. if ri<r0 else r0*ri-r0**2/2. for ri in rabs])


def huber_w(r,r0=1.5):
    '''
    Function w(r) as defined in Egbert & Booker, 1986
    1st derivative of Huber norm divided by r
    w(r) = p'(r)/r
    '''
    rabs = abs(r)
    return np.array([1. if ri <= r0 else r0/ri for ri in rabs])


def regression_m(*args,**kwargs):
    '''
    Alias for current regression_m implementation
    A better regression-m implementation can be replaced here without breaking calls
    '''
    return regression_m_slow(*args,**kwargs)


def regression_m_slow(fct,fcr,tf0=None,sig0=None,r0=1.5,beta=0.7784,convergence=0.000000002):
    '''
    Compute Regression-m estimate and covariance of transfer function
    using regression-m (Egbert & Booker 1986)
    fct: transmitter Fourier coefficients
    fcr: receiver Fourier coefficients
    tf0: initial estimate of transfer function (use LS if None)
    sig0: initial estimate of std (use LS if None)
    r0: cutoff value in Huber norm
    beta: scaling factor for variance estimate
    convergence: threshold of change in tf that dictates when convergence is reached
    '''
    # number of repeated measurements
    num_freq = fct.shape[0]
    N = fct.shape[1]
    tf_hat = np.empty((num_freq,2))
    cov_tf = np.empty((num_freq,2,2))
    all_iter = np.empty(num_freq)
    all_good = np.empty(num_freq)
    # for each frequency
    i_freq = 0
    for frx,ftx in zip(fcr,fct):
        # recast complex system as real system
        # RX = TX*TF -> [real(RX)] = [real(TX), -imag(TX)] * [real(TF)]
        #               [imag(RX)]   [imag(TX),  real(TX)]   [imag(TF)]
        frxr = np.empty(2*N)
        frxr[:N] = frx.real
        frxr[N:] = frx.imag
        ftxr = np.empty((2*N,2))
        ftxr[:N,0] = ftx.real
        ftxr[N:,0] = ftx.imag
        ftxr[:N,1] = -ftx.imag
        ftxr[N:,1] = ftx.real
        # initialize tf estimates and std estimates
        # tfls, varls = LS(ftxr,fcr)
        if tf0 is None:
            tfn = np.dot(np.linalg.inv(np.dot(ftxr.T,ftxr)),np.dot(ftxr.T,frxr))
        else: 
            tfn = tf0[i_freq]
        rx_pre = np.dot(ftxr,tfn)
        res = frxr - rx_pre
        if sig0 is None:
            varls = 1/(2*N-2)*np.dot(res.T,res)
            sign = np.sqrt(varls)
        else:
            sign = sig0[i_freq]
        # objective function starting value
        phin = huber_norm(res/sign)
        # start loop 
        n_iter = 0
        # dummy value to get loop started
        phinm1 = 2*phin/(1.-convergence)+1
        first_loop = True
        while first_loop or phin < phinm1*(1.-convergence):
            first_loop = False
            # modified observations
            rx_mod = rx_pre + huber_w(res/sign,r0)*res
            # new transfer function, LS variance
            tfn = np.dot(np.linalg.inv(np.dot(ftxr.T,ftxr)),np.dot(ftxr.T,rx_mod))
            # predicted observations
            rx_pre = np.dot(ftxr,tfn)
            # residual
            res = frxr - rx_pre
            # new variance
            varn = 1/(2*N-2)/beta*np.dot(res.T,res)
            # new standard deviation
            sign = np.sqrt(varn)
            # old objective function
            phinm1 = phin
            # new objective function
            phin = huber_norm(res/sign)
            n_iter += 1
        tf_hat[i_freq] = tfn
        # print('{} iterations'.format(n_iter))
        all_iter[i_freq] = n_iter
        # compute covariance at the end
        # rnTrn = np.dot(res.T,res)
        # rnTrnEB = np.dot(rx_mod.T,rx_mod)-np.dot(np.dot(ftxr.T,rx_mod).T,tfn)
        # wn = huber_w(res/sign,r0)
        # wrn = wn*res
        # wrnTwrn = np.dot(wrn.T,wrn)
        rnTrnTh = np.dot(rx_mod.T,rx_mod)-np.dot(rx_pre.T,rx_pre)
        # rnTrnneg = 2*np.dot(rx_pre.T,wrn)
        # compare them!
        # print(rnTrn)
        # print('EB: {}'.format(rnTrnEB))
        # print('th: {}'.format(rnTrnTh))
        # print('wr: {}'.format(wrnTwrn))
        # print('Zr: {}'.format(rnTrnneg))
        # print('{}% difference between residuals'.format((rnTrn-rnTrnEB)/rnTrn*100))
        # fraction of "good" observations
        fg = np.sum(res/sign<r0)/2/N
        # print('{}% good obs'.format(fg*100))
        all_good[i_freq] = fg
        cov_tf[i_freq] = rnTrnTh*np.linalg.inv(np.dot(ftxr.T,ftxr))/(2*N-2)/fg/fg
        i_freq += 1
    return tf_hat,cov_tf,all_iter,all_good



def regression_m_buggy(fct,fcr,tf0=None,sig0=None,r0=1.5,convergence=0.02):
    '''
    Compute Regression-m estimate and covariance of transfer function
    using regression-m (Egbert & Booker 1986)
    fct: transmitter Fourier coefficients
    fcr: receiver Fourier coefficients
    tf0: initial estimate of transfer function (use LS if None)
    sig0: initial estimate of std (use LS if None)
    r0: cutoff value in Huber norm
    convergence: threshold of change in tf that dictates when convergence is reached
    '''
    # number of repeated measurements
    N = fct.shape[1]
    # recast complex system as real system
    # RX = TX*TF -> [real(RX)] = [real(TX), -imag(TX)] * [real(TF)]
    #               [imag(RX)]   [imag(TX),  real(TX)]   [imag(TF)]
    # initialize tf estimates and std estimates
    if tf0 is None or sig0 is None:
        tfls, varls = LS(fct,fcr)
        if tf0 is None:
            tf0 = tfls
        if sig0 is None:
            sig0 = np.sqrt(varls)
    # start loop 
    # TODO: work out convergence criteria
    tfn = tf0
    sign = sig0
    #predicted RX
    rx_pre = fct*tfn[:,np.newaxis]
    # residual
    res = fcr - rx_pre
    n_iter = 0
    while False:
        # modified observations
        rx_mod = rx_pre + huber_w(res/sign[:,np.newaxis],r0)*res
        # new transfer function, LS variance
        tfn, varn = LS(fct,rx_mod)
        # robust std
        sign = np.sqrt(varn/beta)
        # predicted observations
        rx_pre = fct*tfn[:,np.newaxis]
        # residual
        res = fcr - rx_pre
        # TODO: don't compute this twice (once in LS)
        n_iter += 1
    #TODO: compute covariance at the end
    # This is maybe wrong - it should be 2x2 (real/imag - see vector/matrix reordering in E&B, '86
    rnTrn = np.sum(np.conjugate(res)*res,axis=1)
    rnTrnEB = np.sum(np.conjugate(rx_mod)*rx_mod,axis=1)-np.sum(np.sum(np.conjugate(fct)*rx_mod)*tfn,axis=1)
    # compare the twain!
    ng = 1# of "good" observations
    cov = 1/(2*N-2)*rnTrn*np.linalg.inv(np.sum(np.conjugate(fct)*fct,axis=1))/ng/ng
    return tfn,cov


def get_calibrated_fc(overlap_data,transmitter_frequency,tx_mask,
                      summary,settings):
    '''
    return calibrated Fourier coefficients, useful for computing transfer functions

    'summary' is a dataframe (a row as returned by z3d_directory.get_z3d_directory_info)
    with these columns:
        sampling_rate
        component
        box_number
        card_number
        antenna_number

    'settings' is a dictionary with these keywords:
        notch_bandwidth
        cycles_per_window
        strides_per_cycle
        window_shape
        cals
        antcal
    settings holds variables that don't change from overlap to overlap
    '''

    samples_per_cycle = int(round(summary.sampling_rate/transmitter_frequency))

    # filter 60 Hz noise
    print('Filtering...')
    if transmitter_frequency > 16:
        # notch filter
        f0 = 60.0
        bw = settings['notch_bandwidth'] # -3 dB bandwidth
        Q = f0/bw
        w0 = f0*2./summary.sampling_rate
        numerator,denominator = signal.iirnotch(w0,Q)
        filt_data = signal.filtfilt(numerator,denominator,overlap_data)
    else:
        # use mean running period
        filt_data = waveform_average_filter(
            overlap_data,sampling_rate=summary.sampling_rate,mask=tx_mask)

    # compute robust running mean
    print('Applying drift correction...')
    taur = robust_running_mean(filt_data,samples_per_cycle)
    # apply drift correction
    drift = filt_data-taur

    # compute short time Fourier transform
    print('Computing Fourier coefficients...')
    samples_per_window = int(round(samples_per_cycle*settings['cycles_per_window']))
    stride=int(round(samples_per_cycle/settings['strides_per_cycle']))
    if settings['window_shape'] == 'tukey_ten_sample_taper':
        samples_per_window += 10
        ten_sample_taper = 20./samples_per_window
        window_shape = ('tukey',ten_sample_taper)
    else:
        window_shape = settings['window_shape']

    stft_str = '{} sampling_rate, {} nperseg, {} noverlap'
    print(stft_str.format(summary.sampling_rate,samples_per_window,samples_per_window-stride))

    f,t,fc = signal.stft(drift,summary.sampling_rate,window=window_shape,nperseg=samples_per_window,
                          noverlap=samples_per_window-stride,boundary=None,padded=False)
    num_windows = len(t)
    num_freq = len(f)

    # apply calibrations
    print('Applying calibrations...')
    zen_cal = settings['cals'].loc[summary.box_number]
    fccb = apply_board_cal(fc,f,zen_cal,summary.card_number,summary.sampling_rate)
    if summary.component[0] == 'H':
        fcc = apply_antenna_cal(fccb,f,settings['antcal'],summary.antenna_number)
    else:
        fcc=fccb
    return (f,t,fcc)


def sliding_window(x,sampling_rate,tx_freq,frequencies=None,cycles_per_window=1,two_strides_per_window=False,window='boxcar'):
    # cycle defined as one transmitter waveform

    if frequencies==None:
        frequencies = [tx_freq]
    try:
        iterator = iter(frequencies)
    except:
        frequencies = [frequencies]

    for freq in frequencies:
        num_samples = len(x)
        # assume sampling_rate = n*transmitter_frequency
        samples_per_cycle = int(round(sampling_rate/tx_freq))
        samples_per_window = samples_per_cycle*cycles_per_window
        # times = np.arange(0,samples_per_window)/sampling_rate
        sinusoid_arguments = 2*np.pi*freq*np.arange(0,samples_per_window)/sampling_rate
        # reshape data into matrix with row size = samples_per_window,
        # ignore any final partial window rather than zero padding
        num_windows = num_samples//samples_per_window
        if two_strides_per_window:
            samples_per_half_window = samples_per_window//2
            overlap_windows1 = x[:num_windows*samples_per_window].reshape(num_windows,samples_per_window)
            num_windows2 = (num_samples-samples_per_half_window)//samples_per_window
            overlap_windows2 = x[samples_per_half_window:samples_per_half_window+num_windows2*samples_per_window].reshape(num_windows2,samples_per_window)
            overlap_windows = np.ma.empty((num_windows+num_windows2,samples_per_window))
            overlap_windows[::2,:] = overlap_windows1
            overlap_windows[1::2,:] = overlap_windows2
        else:
            overlap_windows = x[:num_windows*samples_per_window].reshape(num_windows,samples_per_window)
        # convolve complex exponential with windows
        exp_wave = np.exp(-1j*sinusoid_arguments)
        signal_complex = (overlap_windows*exp_wave).sum(axis=1)
        num_unmasked = overlap_windows.count(axis=1)
        return (signal_complex,num_unmasked)


def tf(tx,rx,sampling_rate,tx_freq,frequencies=None,cycles_per_window=1,cycles_overlap=0,window='boxcar'):
    # transfer function computed for a sliding window
    tx_complex = sliding_window(tx,sampling_rate,tx_freq,frequencies,cycles_per_window,cycles_overlap,window)
    rx_complex = sliding_window(rx,sampling_rate,tx_freq,frequencies,cycles_per_window,cycles_overlap,window)
    normalized = rx_complex/tx_complex
    return normalized


