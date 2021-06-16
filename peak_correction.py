# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 14:47:25 2021

@author: cfai2
"""

from interferogram_functions import prep_interferogram,  FFT_intr, import_MAP
import matplotlib.pyplot as plt
import numpy as np
from numpy import savetxt
import os

def time_to_index(time_data, time):
    return np.where(time_data <= time)[0][-1]

def find_peaks(raw_data, time_data, start_time=-1, end_time=20, 
               search_radius=3, peak_thr=1.5,peak_radius=10):
    """
    Search for peaks by identifying locations in raw_data where the data
    suddenly decreases after increasing for (search_radius) consecutive points.
    
    These locations are known as peak indices.

    Parameters
    ----------
    raw_data : np.ndarray
        PL data array
    time_data : 1D np.ndarray
        Vector of time steps.
    start_time : float, optional
        Zero-Shifted time where peak search should begin. 
        To ensure that the initial peak is included, choose a small negative time.
        The default is -1.
    end_time : float, optional
        Zero-shifted time where peak search should stop. 
        Choose a value such that the noisy data is not included. 
        The default is 20.
    search_radius : int, optional
        Number of consecutive increasing points needed to register a peak. 
        Higher values make the search less sensitive.
        The default is 3.
    peak_thr : float, optional
        How many times larger the final value in a peak compared to the first value
        must be for the peak to be considered a peak. 
        Higher values make the search less sensitive.
        The default is 1.5.
    peak_radius : int, optional
        A peak is defined as each peak_index +/- this many indices. The default is 5.

    Returns
    -------
    peaks : list(list(int, int, float))
        List of info for each peak - each peak stores a starting time index,
        end time index, and peak value

    """
    
    start_i = time_to_index(time_data, start_time)
    end_i = time_to_index(time_data, end_time)
    peaks = []
    peak_base = 0
    count = 0
    sign = 0
    for i in range(start_i, end_i):
        delta = (raw_data[i+1] - raw_data[i])
        if delta < 0:
            # Streak broken
            streak_long_enough = (count >= search_radius)
            peak_large_enough = raw_data[i] > peak_thr * peak_base
            if streak_long_enough and peak_large_enough:
                peaks.append(i)
            count = 0
            peak_base = 0
            
        elif delta > 0:
            if count == 0: peak_base = raw_data[i]
            # Continue streak
            count += 1
    
    for p, peak in enumerate(peaks):
        peaks[p] = [peak - peak_radius, peak + 10*peak_radius, raw_data[peak]]
    return peaks

def subtract_peaks(raw_data, peaks, reduce=0.9):
    largest_peak = [-1,-1,-1]
    for peak_info in peaks:
        if peak_info[2] > largest_peak[2]:
            largest_peak = peak_info
            
    if (largest_peak != peaks[0]):
        print("WARNING: largest peak is not first; start_time for find_peak() may be too low")
            
    largest_peak_values = raw_data[largest_peak[0]:largest_peak[1]]
    for peak_info in peaks[1:]:
        raw_data[peak_info[0]:peak_info[1]] -= largest_peak_values * (reduce * peak_info[2] / largest_peak[2])

    return


def find_peak(raw_data, peak_inds):
    peak = raw_data[peak_inds[0]:peak_inds[1]]
    plt.figure(5)
    plt.plot(time_data[peak_inds[0]:peak_inds[1]], peak)
    return peak

def subtract_peak(raw_data, peak, loc, shift):
    plt.figure(5)
    plt.plot(time_data[loc:loc+len(peak)], peak*shift)
    raw_data[loc:loc+len(peak)] -= peak * shift
    
    
    return
path = r"C:\Users\cfai2\Documents\src\Interferogram_FFT\20210615\105836"

BKGsub = True               #Background Subtract - Generally True
bkg_limit = -3              #ns before the TRPL peak to average the background data up to - see plot
TRPLmin_OM = 1e-6           #How many orders of magnitude to plot down in y-scale for TRPL curve


pos_data, time_data, map_data = import_MAP(path)

#Background Subtract
t_max = time_data[np.array(np.where(np.mean(map_data,axis=0)==np.max(np.mean(map_data,axis=0)))[0],dtype="int")]
time_data=time_data-t_max
BKGrange = np.array([time_data[0],bkg_limit],dtype='float')  #ns
if BKGsub:
    index = [(np.abs(time_data-np.min(BKGrange))).argmin(),(np.abs(time_data-np.max(BKGrange))).argmin()]
    BKGval = np.mean(map_data[:,np.min(index):np.max(index)],axis=1)
    map_data = map_data - np.reshape(BKGval, (len(BKGval), 1))
    
integralTRPL = np.sum(map_data,axis=0)
    
#Plot Wavelength Averaged decay to determine time range for background subtraction
plt.figure(3, dpi=120)
plt.plot(time_data,np.mean(map_data,axis=0))
plt.axvspan(np.min(BKGrange),np.max(BKGrange),facecolor='r',alpha=0.2)
plt.xlim(right=2)
#plt.xlabel('Time / ns')
#plt.ylabel('Counts / a.u.')
plt.yscale('log')

#Plot Full TRPL
plt.figure(4, dpi=120)
plt.plot(time_data,integralTRPL)
plt.ylim(np.max(integralTRPL)*TRPLmin_OM,2*np.max(integralTRPL))
#plt.xlabel('Time / ns')
#plt.ylabel('Counts / a.u.')
#plt.title("Integral TRPL")
plt.yscale('log')

peaks = find_peaks(integralTRPL, time_data)
print(peaks)
plt.figure(4)
for peak in peaks:
    plt.plot(time_data[peak[0]:peak[1]], integralTRPL[peak[0]:peak[1]])
    
subtract_peaks(integralTRPL, peaks)

plt.figure(6)
plt.title("Corrected #2")
plt.yscale('log')
plt.ylim(np.max(integralTRPL)*TRPLmin_OM,2*np.max(integralTRPL))
plt.plot(time_data, integralTRPL)

import sys
sys.exit()

peak = find_peak(integralTRPL, [878,1028])
subtract_peak(integralTRPL, peak, 1294, 0.0008)

subtract_peak(integralTRPL, peak, 1368, 0.00315)
plt.figure(6)
plt.title("Corrected #2")
plt.yscale('log')
plt.ylim(np.max(integralTRPL)*TRPLmin_OM,2*np.max(integralTRPL))
plt.plot(time_data, integralTRPL)
