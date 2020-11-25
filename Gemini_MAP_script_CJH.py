# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 23:08:58 2020

@author: Chuck
"""
from interferogram_functions import FFT_map, import_MAP, prep_map
import matplotlib.pyplot as plt
import numpy as np
import time
import h5py
import os

startTime = time.time()
path = r"C:\Users\Chuck\Dropbox (UFL)\UF\TRPL Computer\152430"
save_data = True
outputfilename = path + "\\" + os.path.split(path)[-1] + '_TRES.h5'

pos_data, time_data, map_data = import_MAP(path)

#trim time-scale
rangeval = [0,100]  #ns
index = [(np.abs(time_data-np.min(rangeval))).argmin(),(np.abs(time_data-np.max(rangeval))).argmin()]
map_data = map_data[:,np.min(index):np.max(index)]
time_data=time_data[np.min(index):np.max(index)]

#Manual Shift of Peak
pos_data = pos_data - 0.15619633484742135    #Taken from shift_factor output in the _INTR analysis script for this data

#Background Subtract TRPL Curves
BKGsub = "True"
BKGrange = [0,5]  #ns

if BKGsub == "True":
    index = [(np.abs(time_data-np.min(BKGrange))).argmin(),(np.abs(time_data-np.max(BKGrange))).argmin()]
    BKGval = np.mean(map_data[:,np.min(index):np.max(index)],axis=1)
    for i in range(len(BKGval)):
        map_data[i,:] = map_data[i,:]-BKGval[i]

apodization_width=1.75
apod_type="BH"    # "None" "Gauss" "Triangle" "Boxcar" "BH"
resample="True"
resample_factor=4
shift="True"
pad_test="True"
padfactor=4
mean_sub = "True"

# =============================================================================
# #Plot Raw Data (Background Subtracted and Shifted)
# raw_timemesh, raw_posmesh = np.meshgrid(time_data,pos_data)
# plt.figure(1, dpi=120)
# plt.contourf(raw_posmesh,raw_timemesh,np.log(map_data))
# plt.ylim(5,25)
# plt.ylabel('Time / ns')
# plt.xlabel('Position / mm')
# =============================================================================

preFFT_pos, preFFT_map = prep_map(pos_data,map_data,apodization_width,apod_type=apod_type,resample=resample,resample_factor=resample_factor,shift="False",pad_test=pad_test,padfactor=padfactor,mean_sub=mean_sub,plots="True")
print("Took {} sec".format(time.time() - startTime))

#Perform FFT
wave, build_TRES = FFT_intr(preFFT_pos, preFFT_map,plots="False",scale="linear",correct="True")
build_TRES=np.array(build_TRES,dtype="float")
print("Took {} sec".format(time.time() - startTime))

# Plot the results
timemesh, wavemesh = np.meshgrid(time_data,wave)
start_wave = 700
end_wave = 1000

plt.figure(2, dpi=120)
plt.contourf(wavemesh,timemesh,np.log(build_TRES))
plt.xlim(start_wave,end_wave)
plt.ylabel('Time / ns')
plt.xlabel('Wavelength / nm')
plt.colorbar()

#Plot averged PL over given range
AveragePL = "False"
rangeval = [0,100]  #ns

plt.figure(3, dpi=120)
if AveragePL == "True":
    index = [(np.abs(time_data-np.min(rangeval))).argmin(),(np.abs(time_data-np.max(rangeval))).argmin()]
    plt.plot(wave,np.mean(build_TRES[:,np.min(index):np.max(index)],axis=1))
elif AveragePL == "False":
    plt.plot(wave,np.mean(build_TRES,axis=1))
plt.xlim(start_wave,end_wave)
plt.ylabel('Counts / a.u.')
plt.xlabel('Wavelength / nm')
plt.yscale('linear')

#Plot averged PL decay over given range
AverageTRPL = "False"
rangeval = [800,900]  #nm

plt.figure(4, dpi=120)
if AverageTRPL == "True":
    index = [(np.abs(wave-np.min(rangeval))).argmin(),(np.abs(wave-np.max(rangeval))).argmin()]
    plt.plot(time_data,np.mean(build_TRES[np.min(index):np.max(index),:],axis=0))
elif AverageTRPL == "False":
    plt.plot(time_data,np.mean(build_TRES,axis=0))
plt.xlabel('Time / ns')
plt.ylabel('Counts / a.u.')
plt.yscale('log')

#Write 2D Data Set
if save_data:
    hf = h5py.File(outputfilename,'w')
    hf.create_dataset('TRES Data', data=build_TRES)
    hf.create_dataset('Time Data', data=time_data)
    hf.create_dataset('Wavelength', data=wave)
    hf.create_dataset('Metadata', data=[[start_wave,end_wave]])
    hf.close()

# =============================================================================
# #Make composite Plot??
# fig = plt.figure(figsize=(6, 6))
# grid = plt.GridSpec(4, 4, hspace=0.2, wspace=0.2)
# main_ax = fig.add_subplot(grid[:-1, 1:])
# y_hist = fig.add_subplot(grid[:-1, 0], xticklabels=[], sharey=main_ax)
# x_hist = fig.add_subplot(grid[-1, 1:], yticklabels=[], sharex=main_ax)
# =============================================================================
