# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 23:08:58 2020

@author: Chuck
"""
from interferogram_functions import FFT_intr, import_MAP, prep_map
import matplotlib.pyplot as plt
import numpy as np

path = r"C:\Users\Chuck\Dropbox (UFL)\Hages Lab Files\TRPL Data\Ruiquan\20201104\153345"
pos_data, time_data, map_data = import_MAP(path)
pos_data = pos_data - 0.03478240614951511    #Taken from shift_factor output in the _INTR analysis script for this data

apodization_width=0.6
apod_type="BH"    # "None" "Gauss" "Triangle" "Boxcar" "BH"
resample="True"
resample_factor=4
pad_test="True"
padfactor=4
mean_sub = "True"         

#Plot Raw Data
raw_timemesh, raw_posmesh = np.meshgrid(time_data,pos_data)
plt.figure(1, dpi=120)
plt.contourf(raw_posmesh,raw_timemesh,np.log(map_data))
plt.ylim(5,25)
plt.ylabel('Time / ns')
plt.xlabel('Position / mm')

preFFT_pos, preFFT_map = prep_map(pos_data,map_data,apodization_width,apod_type=apod_type,resample=resample,resample_factor=resample_factor,shift="False",pad_test=pad_test,padfactor=padfactor,mean_sub=mean_sub,plots="True")

#trim time-scale data because it is dreadfully slow
rangeval = [5,25]  #ns
index = [(np.abs(time_data-np.min(rangeval))).argmin(),(np.abs(time_data-np.max(rangeval))).argmin()]
preFFT_map = preFFT_map[:,np.min(index):np.max(index)]
time_data=time_data[np.min(index):np.max(index)]

#Perform FFT
build_TRES=[]
for i in range(preFFT_map.shape[1]):
    wave, FFT_intr_trim = FFT_intr(preFFT_pos,preFFT_map[:,i],plots="False",scale="linear",correct="True")
    build_TRES.append(FFT_intr_trim)
build_TRES=np.array(build_TRES,dtype="float").T

#Background Subtract
BKGsub = "False"
BKGrange = [5,10]  #ns

if BKGsub == "True":
    index = [(np.abs(time_data-np.min(rangeval))).argmin(),(np.abs(time_data-np.max(rangeval))).argmin()]
    BKGval = np.mean(build_TRES[:,np.min(index):np.max(index)])
    build_TRES = build_TRES-BKGval

# Plot the results
timemesh, wavemesh = np.meshgrid(time_data,wave)
start_wave = 400
end_wave = 1000

plt.figure(2, dpi=120)
plt.contourf(wavemesh,timemesh,build_TRES,100,vmin=0)
plt.xlim(start_wave,end_wave)
plt.ylabel('Time / ns')
plt.xlabel('Wavelength / nm')
plt.colorbar()

#Plot averged PL over given range
AveragePL = "True"
rangeval = [7.5,25]  #ns

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
AverageTRPL = "True"
rangeval = [610,750]  #nm

plt.figure(4, dpi=120)
if AverageTRPL == "True":
    index = [(np.abs(wave-np.min(rangeval))).argmin(),(np.abs(wave-np.max(rangeval))).argmin()]
    plt.plot(time_data,np.mean(build_TRES[np.min(index):np.max(index),:],axis=0))
elif AverageTRPL == "False":
    plt.plot(time_data,np.mean(build_TRES,axis=0))
plt.xlabel('Time / ns')
plt.ylabel('Counts / a.u.')
plt.yscale('log')

#Make composite Plot??
