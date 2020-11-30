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

path = r"C:\Users\Chuck\Desktop\20_11_24\173735"
importfilename = path + "\\" + os.path.split(path)[-1] + '_TRES.h5'
export = False

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
BKGrange = [0,4.5]  #ns
if BKG_TRPL:
    index = [(np.abs(time_data-np.min(BKGrange))).argmin(),(np.abs(time_data-np.max(BKGrange))).argmin()]
    TRPLdata = TRPLdata-np.mean(TRPLdata[np.min(index):np.max(index)])

AveragePL = False
rangeval = [np.max(np.mean(TRES,axis=0)),25]  #ns
if AveragePL:
    index = [(np.abs(time_data-np.min(rangeval))).argmin(),(np.abs(time_data-np.max(rangeval))).argmin()]
    PLdata = np.sum(TRES[:,np.min(index):np.max(index)],axis=1)
else:
    PLdata = np.sum(TRES,axis=1)

fig = plt.figure(0,dpi=120)
grid = plt.GridSpec(2, 3, height_ratios=[3, 1],width_ratios=[1,2,0.44],wspace=0.05,hspace=0.05)

main_ax = fig.add_subplot(grid[:-1, 1:])
main_ax.set_xlim(wave_range)
main_ax.set_ylim(time_range)
#cs = main_ax.contourf(wavemesh,timemesh,np.log(TRES),100,vmin=0, cmap='plasma')  
cs = main_ax.contourf(wavemesh,timemesh,TRES,60,locator=ticker.LogLocator(),vmin=1e1,cmap='plasma')  
main_ax.tick_params(bottom='off')
main_ax.label_outer()
plt.colorbar(cs, ax=main_ax)

TRPL = fig.add_subplot(grid[:-1, 0], xticklabels=[],sharey=main_ax)
TRPL.invert_xaxis()
TRPL.set_xscale('log')
TRPL.set_xlim(2*np.max(TRPLdata),np.max(TRPLdata)*TRPLmin_OM)
TRPL.plot(TRPLdata,time_data)
TRPL.set(ylabel='Time / ns')

PL = fig.add_subplot(grid[-1, 1], yticklabels=[], sharex=main_ax)
PL.plot(wave,PLdata)
PL.set(xlabel='Wavelength / nm')

if export:
    if not os.path.isdir("Exports"):
        os.mkdir("Exports")
    np.savetxt(r"Exports\{}_AverageTRPL.csv".format(os.path.split(path)[-1]), TRPLdata)
    np.savetxt(r"Exports\{}_AveragePL.csv".format(os.path.split(path)[-1]), PLdata)
    