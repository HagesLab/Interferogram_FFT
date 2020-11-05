# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 23:08:58 2020

@author: Chuck
"""
from interferogram_functions import FFT_intr, import_MAP, prep_map
import matplotlib.pyplot as plt
import numpy as np

path = r"C:\Users\Chuck\Dropbox (UFL)\UF\TRPL Computer\Measurments\20201016\145521"
pos_data, time_data, map_data = import_MAP(path)
pos_data = pos_data - 0.03478240614951511    #Taken from shift_factor output in the _INTR analysis script for this data
apodization_width=0.5
apod_type="BH"    # "None" "Gauss" "Triangle" "Boxcar" "BH"
resample="True"
resample_factor=2

#ADD MANUAL SHIFT by shift factor - currently done above manually
#WHAT ABOUT INTENSITY CORRECTION? Ask NIREOS
shift="False"
pad_test="True"
padfactor=4
mean_sub = "True"       

#Plot Raw Data
raw_timemesh, raw_posmesh = np.meshgrid(time_data,pos_data)
plt.figure(1, dpi=120)
plt.contourf(raw_posmesh,raw_timemesh,map_data)
plt.ylabel('Time / ns')
plt.xlabel('Position / mm')

preFFT_pos, preFFT_map = prep_map(pos_data,map_data,apodization_width,apod_type=apod_type,resample=resample,resample_factor=resample_factor,shift=shift,pad_test=pad_test,padfactor=padfactor,mean_sub=mean_sub,plots="True")

#trim time-scale data because it is dreadfully slow
#Add trimming based on time values in time_data insead of position
preFFT_map = preFFT_map[:,220:320]
time_data=time_data[220:320]

#Perform FFT
build_TRES=[]
for i in range(preFFT_map.shape[1]):
    wave, FFT_intr_trim = FFT_intr(preFFT_pos,preFFT_map[:,i],plots="False",scale="linear",correct="True")
    build_TRES.append(FFT_intr_trim)
build_TRES=np.array(build_TRES,dtype="float").T


#NEED TO ADD BACKGROND SUBTRACTION - CHECK SINGLE HISTOGRAM
timemesh, wavemesh = np.meshgrid(time_data,wave)
# Plot the results
start_wave = 400
end_wave = 1000

plt.figure(2, dpi=120)
plt.contourf(wavemesh,timemesh,np.log(build_TRES),100,vmin=0)
plt.xlim(start_wave,end_wave)
plt.ylabel('Time / ns')
plt.xlabel('Wavelength / nm')
plt.colorbar()

#Can I add a custom WL range to average over?
plt.figure(3, dpi=120)
plt.plot(wave,np.mean(build_TRES,axis=1))
plt.xlim(start_wave,end_wave)
plt.ylabel('Counts / a.u.')
plt.xlabel('Wavelength / nm')
plt.yscale('linear')

#Can I add a custom time range to average over?
plt.figure(4, dpi=120)
plt.plot(time_data,np.mean(build_TRES,axis=0))
plt.xlabel('Time / ns')
plt.ylabel('Counts / a.u.')
plt.yscale('log')

#Make composite Plot??