# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 14:53:37 2020

@author: Chuck
"""

#import libraries
import numpy as np
import csv
import matplotlib.pylab as plt

Plot_Title = r'CdTe 9-8 @ 750 nm'

path1 = r"\\ad.ufl.edu\che\Hages-Lab\TRPL Data\Hages\Ferekides\20_12_1\161903\161903_AverageTRPL_800_to_875nm.csv"
tag1 = r'9-8-1-A'

path2 = r"\\ad.ufl.edu\che\Hages-Lab\TRPL Data\Hages\Ferekides\20_12_1\155433\155433_AverageTRPL_800_to_875nm.csv"
tag2 = r'9-8-1-B'

path3 = r"\\ad.ufl.edu\che\Hages-Lab\TRPL Data\Hages\Ferekides\20_11_24\152430\152430_AverageTRPL_800_to_875nm.csv"
tag3 = r'9-8-1-C'

path4 = r"\\ad.ufl.edu\che\Hages-Lab\TRPL Data\Hages\Ferekides\20_11_24\180627\180627_AverageTRPL_800_to_875nm.csv"
tag4 = r'9-8-1-D'


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

TRPL1, tag1 = TRPL_Manipulate(path1,tag1,norm,zero_shift,t_shift=0)
TRPL2, tag2 = TRPL_Manipulate(path2,tag2,norm,zero_shift,t_shift=0)
TRPL3, tag3 = TRPL_Manipulate(path3,tag3,norm,zero_shift,t_shift=0)
TRPL4, tag4 = TRPL_Manipulate(path4,tag4,norm,zero_shift,t_shift=0)
    
plt.figure(0, dpi=120)
plt.title(Plot_Title)
plt.yscale('log')
plt.xlabel(r'Delay / ns')
plt.ylabel(r'Normalized Counts / a.u.')
plt.ylim(1e-3,2)
plt.xlim(-2, 40)
plt.plot(TRPL1[:,0],TRPL1[:,1],label = tag1)
plt.plot(TRPL2[:,0],TRPL2[:,1],label = tag2)
plt.plot(TRPL3[:,0],TRPL3[:,1],label = tag3)
plt.plot(TRPL4[:,0],TRPL4[:,1],label = tag4)
plt.legend()

