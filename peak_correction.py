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

#Plot Full TRPL
plt.figure(5, dpi=120)
plt.plot(time_data,integralTRPL)
plt.ylim(np.max(integralTRPL)*TRPLmin_OM,2*np.max(integralTRPL))
#plt.xlabel('Time / ns')
#plt.ylabel('Counts / a.u.')
#plt.title("Integral TRPL")
plt.yscale('log')

peak = find_peak(integralTRPL, [878,1028])
subtract_peak(integralTRPL, peak, 1294, 0.0008)

subtract_peak(integralTRPL, peak, 1368, 0.00315)
plt.figure(6)
plt.title("Corrected #2")
plt.yscale('log')
plt.ylim(np.max(integralTRPL)*TRPLmin_OM,2*np.max(integralTRPL))
plt.plot(time_data, integralTRPL)
