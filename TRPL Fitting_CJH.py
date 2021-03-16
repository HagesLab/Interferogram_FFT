# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 17:33:10 2021

@author: c.hages
"""
#import libraries
import numpy as np
import matplotlib.pylab as plt
from scipy.optimize import curve_fit
import h5py
import os


path = r"C:\Users\c.hages\Dropbox (UFL)\UF\TRPL Computer\Aaron\144620"

importfilename = path + "\\" + os.path.split(path)[-1] + '_TRES.h5'
hf = h5py.File(importfilename, 'r')
TRPL_data = np.array(hf.get('TRPL'))
time_data = np.array(hf.get('Time Data'))
hf.close()

def Fit_1exp(TRPL_data,time_data,fitrange):

    def Exp1(time,A,tau):
        return -time/tau + np.log(A)

    #trim-data
    low_index, high_index = (np.abs(time_data-np.min(fitrange))).argmin() , (np.abs(time_data-np.max(fitrange))).argmin()
    TRPL_fit = TRPL_data[low_index:high_index]
    time_fit = time_data[low_index:high_index]

    popt, pcov = curve_fit(Exp1,time_fit,np.log(np.abs(TRPL_fit)))
    perr = np.sqrt(np.diag(pcov))

    TRPL_out = np.exp(Exp1(time_fit,*popt))
    label =  r'$\tau:\ $' + np.array2string(popt[1], precision=2, separator=',', suppress_small=True) + ' ns'

    return TRPL_out, time_fit, label, popt, perr

# =============================================================================
# def Fit_2exp(TRPL_data,fitrange):
#
#     def Exp2(time,A1,A2,tau_1,tau_2):
#         return np.log(A1*np.exp(-time/tau_1) + A2*np.exp(-time/tau_2))
#
#     #trim-data
#     low_index, high_index = (np.abs(TRPL_data[:,0]-np.min(fitrange))).argmin() , (np.abs(TRPL_data[:,0]-np.max(fitrange))).argmin()
#     TRPL_fit = TRPL_data[low_index:high_index,:]
#
#     popt, pcov = curve_fit(Exp2,TRPL_fit[:,0],np.log(np.abs(TRPL_fit[:,1])))
#     perr = np.sqrt(np.diag(pcov))
#
#     y_values = np.exp(Exp2(TRPL_fit[:,0],*popt))
#     output = np.vstack((TRPL_fit[:,0],y_values)).T
#     label =  r'$\tau_1:\ $' + np.array2string(np.min(popt[2:4]), precision=2, separator=',', suppress_small=True) + ' ns; ' + r'$\tau_2:\ $' + np.array2string(np.max(popt[2:4]), precision=2, separator=',', suppress_small=True) + ' ns'
#
#     return output, label, popt, perr
# =============================================================================

fit_range = [2,15]
TRPL_fit, time_fit, fit_label, popt, perr = Fit_1exp(TRPL_data,time_data,fit_range)
plt.figure(0, dpi=120)
plt.yscale('log')
plt.xlabel(r'Delay / ns')
plt.ylabel(r'Normalized Counts / a.u.')
# =============================================================================
# plt.ylim(5e-5,2)
# plt.xlim(-1, 40)
# =============================================================================
plt.plot(time_data,TRPL_data)
plt.plot(time_fit,TRPL_fit,'k--',label = fit_label)
plt.legend()

