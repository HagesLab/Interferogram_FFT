# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 17:31:29 2020

@author: c.hages
"""
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import ticker
import h5py

path = r"20_12_1\161903"
importfilename = path + "\\" + os.path.split(path)[-1] + '_TRES.h5'
export = True

hf = h5py.File(importfilename, 'r')
TRES = np.array(hf.get('TRES Data'))
time_data=np.array(hf.get('Time Data'))
wave=np.array(hf.get('Wavelength'))
hf.close()

timemesh, wavemesh = np.meshgrid(time_data,wave)
wave_range =  [700 , 950]
time_range = [0,40]
TRPLmin_OM = 1e-3

#TimeShit to Peak - use this for BKG sub

AverageTRPL = True
rangeval = [800,875]  #nm
if AverageTRPL:
    index = [(np.abs(wave-np.min(rangeval))).argmin(),(np.abs(wave-np.max(rangeval))).argmin()]
    TRPLdata = np.sum(TRES[np.min(index):np.max(index),:],axis=0)
else:
    TRPLdata = np.sum(TRES,axis=0)

BKG_TRPL = True
BKGrange = [0,4.6]  #ns
if BKG_TRPL:
    index = [(np.abs(time_data-np.min(BKGrange))).argmin(),(np.abs(time_data-np.max(BKGrange))).argmin()]
    TRPLdata = TRPLdata-np.mean(TRPLdata[np.min(index):np.max(index)])

AveragePL = False
t_rangeval = [np.max(np.mean(TRES,axis=0)),25]  #ns
if AveragePL:
    index = [(np.abs(time_data-np.min(t_rangeval))).argmin(),(np.abs(time_data-np.max(t_rangeval))).argmin()]
    PLdata = np.sum(TRES[:,np.min(index):np.max(index)],axis=1)
else:
    PLdata = np.sum(TRES,axis=1)

fig = plt.figure(0,dpi=120)
grid = plt.GridSpec(2, 3, height_ratios=[3, 1],width_ratios=[1.5,2,0.43],wspace=0.05,hspace=0.05)

main_ax = fig.add_subplot(grid[:-1, 1:])
main_ax.set_title(path)
main_ax.set_xlim(wave_range)
main_ax.set_ylim(time_range)
#cs = main_ax.contourf(wavemesh,timemesh,np.log10(TRES),60,vmin=1, cmap='plasma')  
#cs = main_ax.contourf(wavemesh,timemesh,TRES,locator=ticker.LogLocator(base=10, numticks=30),vmin=1e1,cmap='plasma')  

# Instructive Matplotlib moment here:
import matplotlib.colors as mc
# This tells the colorbar how to distribute its colors - by default it's a linear scale, which, if your data is log scale,
# will assign one color to the highest order of magnitude datapoint and a second color to everything else
# Uncomment the BoundaryNorm, which is a linear scale, and comment out the LogNorm to see what I mean
#norm = mc.BoundaryNorm(np.geomspace(1e1, 1e6, 41), ncolors=41)
norm = mc.LogNorm(1.5e1, 1e4)

# Levels control how many different color shades there are: more levels = smoother gradient
cs = main_ax.contourf(wavemesh,timemesh,TRES,levels=np.geomspace(1e1, 1e4, 41),norm=norm ,cmap='plasma', extend='min')  
main_ax.tick_params(bottom='off')
main_ax.label_outer()
# ticks control what axis values are physically written on the colorbar
cbar = plt.colorbar(cs, ax=main_ax, ticks=np.geomspace(1e1, 1e5, 5))

TRPL = fig.add_subplot(grid[:-1, 0], xticklabels=[],sharey=main_ax)
TRPL.invert_xaxis()
TRPL.set_xscale('log')
TRPL.set_xlim(2*np.max(TRPLdata),1.5e2)
TRPL.plot(TRPLdata,time_data)
TRPL.set(ylabel='Time / ns', xlabel="TRPL")

PL = fig.add_subplot(grid[-1, 1], yticklabels=[], sharex=main_ax)
PL.plot(wave,PLdata)
PL.set(xlabel='Wavelength / nm')

if export:
    fig.savefig(r"{}\{}py.png".format(path, os.path.split(path)[-1]))
    np.savetxt(r"{}\{}_AverageTRPL_{}_to_{}nm.csv".format(path, os.path.split(path)[-1], *rangeval), np.vstack((time_data,TRPLdata)).T, delimiter=',')
    np.savetxt(r"{}\{}_AveragePL_{}_to_{}ns.csv".format(path, os.path.split(path)[-1], *t_rangeval), np.vstack((wave,PLdata)).T, delimiter=',')
    