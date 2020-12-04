# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 14:53:37 2020

@author: Chuck
"""

#import libraries
import numpy as np
import csv
import matplotlib.pylab as plt
from scipy.optimize import curve_fit

Plot_Title = r'CdTe @ 750 nm'

path = r"\\ad.ufl.edu\che\Hages-Lab\TRPL Data\Hages\Ferekides\20_11_24\180627\180627_AverageTRPL_800_to_875nm.csv"
tag = r'9-8-1-D'

norm = True
zero_shift = True

def TRPL_Manipulate(path,tag,norm=False,zero_shift=False,t_shift=0):

    with open(path,'r') as i:     #open a file in directory of this script for reading 
        TRPL = np.array(list(csv.reader(i,delimiter=",")),dtype=float)   #make a list of data in file
    
    if norm:
       TRPL_max = np.max(TRPL[:,1]) 
       TRPL[:,1] = TRPL[:,1]/TRPL_max

    if zero_shift:
        TRPL_max = np.max(TRPL[:,1])
        max_index = (np.abs(TRPL[:,1]-TRPL_max)).argmin()
        time_peak = TRPL[max_index,0]
        TRPL[:,0] = TRPL[:,0] - time_peak
    
    if t_shift != 0:
        TRPL[:,0] = TRPL[:,0] + t_shift
       
    return TRPL, tag

def Fit_1exp(TRPL_data,fitrange):
   
    def Exp1(time,A,tau):
        return -time/tau + np.log(A)
    
    #trim-data
    low_index, high_index = (np.abs(TRPL_data[:,0]-np.min(fitrange))).argmin() , (np.abs(TRPL_data[:,0]-np.max(fitrange))).argmin()
    TRPL_fit = TRPL_data[low_index:high_index,:]
    
    popt, pcov = curve_fit(Exp1,TRPL_fit[:,0],np.log(np.abs(TRPL_fit[:,1])))
    perr = np.sqrt(np.diag(pcov))
    
    y_values = np.exp(Exp1(TRPL_fit[:,0],*popt))
    output = np.vstack((TRPL_fit[:,0],y_values)).T
    label =  r'$\tau:\ $' + np.array2string(popt[1], precision=2, separator=',', suppress_small=True) + ' ns'
    
    return output, label, popt, perr 

def Fit_2exp(TRPL_data,fitrange):
   
    def Exp2(time,A1,A2,tau_1,tau_2):
        return np.log(A1*np.exp(-time/tau_1) + A2*np.exp(-time/tau_2))
    
    #trim-data
    low_index, high_index = (np.abs(TRPL_data[:,0]-np.min(fitrange))).argmin() , (np.abs(TRPL_data[:,0]-np.max(fitrange))).argmin()
    TRPL_fit = TRPL_data[low_index:high_index,:]
    
    popt, pcov = curve_fit(Exp2,TRPL_fit[:,0],np.log(np.abs(TRPL_fit[:,1])))
    perr = np.sqrt(np.diag(pcov))
    
    y_values = np.exp(Exp2(TRPL_fit[:,0],*popt))
    output = np.vstack((TRPL_fit[:,0],y_values)).T
    label =  r'$\tau_1:\ $' + np.array2string(np.min(popt[2:4]), precision=2, separator=',', suppress_small=True) + ' ns; ' + r'$\tau_2:\ $' + np.array2string(np.max(popt[2:4]), precision=2, separator=',', suppress_small=True) + ' ns'
    
    return output, label, popt, perr

TRPL, tag = TRPL_Manipulate(path,tag,norm,zero_shift,t_shift=0)

fit_range = [0,25]
TRPL_fit, fit_label, popt, perr = Fit_2exp(TRPL,fit_range)


plt.figure(0, dpi=120)
plt.title(Plot_Title)
plt.yscale('log')
plt.xlabel(r'Delay / ns')
plt.ylabel(r'Normalized Counts / a.u.')
plt.ylim(5e-5,2)
plt.xlim(-1, 40)
plt.plot(TRPL[:,0],TRPL[:,1],label = tag)
plt.plot(TRPL_fit[:,0],TRPL_fit[:,1],'k--',label = fit_label)
plt.legend()

