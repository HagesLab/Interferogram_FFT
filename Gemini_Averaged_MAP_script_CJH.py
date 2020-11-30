# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 11:37:33 2020

@author: Chuck
"""

from interferogram_functions import prep_interferogram, import_INTR, FFT_intr, import_MAP
import matplotlib.pyplot as plt
import os
import numpy as np

path = r"C:\Users\Chuck\Desktop\20_11_24\173735"
pos_data, time_data, map_data = import_MAP(path)

#Background Subtract TRPL Curves
BKGsub = True
BKGrange = [0,4.5]  #ns

if BKGsub:
    index = [(np.abs(time_data-np.min(BKGrange))).argmin(),(np.abs(time_data-np.max(BKGrange))).argmin()]
    BKGval = np.mean(map_data[:,np.min(index):np.max(index)],axis=1)

    map_data = map_data - np.reshape(BKGval, (len(BKGval), 1))
    
#Plot Wavelength Averaged decay to determine time range for background subtraction
plt.figure(0, dpi=120)
plt.plot(time_data,np.mean(map_data,axis=0))
plt.xlim(0,10)
plt.xlabel('Time / ns')
plt.ylabel('Counts / a.u.')
plt.yscale('log')
    

#New Time Averaged data from MAP
AVG_map_data = np.sum(map_data,axis=1)

apodization_width=1.75
apod_type="BH"    # "None" "Gauss" "Triangle" "Boxcar" "BH"
resample="True"
resample_factor=4
shift="True"
pad_test="True"
padfactor=4
mean_sub = "True"

preFFT_pos, preFFT_data, shiftfactor = prep_interferogram(pos_data,AVG_map_data,apodization_width,apod_type=apod_type,resample=resample,resample_factor=resample_factor,shift=shift,pad_test=pad_test,padfactor=padfactor,mean_sub=mean_sub,plots="True")
wave, FFT_intr_trim = FFT_intr(preFFT_pos,preFFT_data,plots="True",scale="linear",correct="True")

start_wave = 700
end_wave = 1000
plt.figure(3, dpi=120)
plt.plot(wave,FFT_intr_trim)
plt.xlim(start_wave,end_wave)
plt.yscale('linear')
plt.show()

# Best to turn this on only when you have found the desired params
save_params = True

if save_params:
    params = {"apod_width":apodization_width, "apod_type":apod_type, "do_resample":resample, "resample_factor":resample_factor,
              "do_shift":shift, "do_padzeros":pad_test, "pad_factor":padfactor, "do_mean_sub":mean_sub, "shift_factor":shiftfactor}
    with open(r"{}\Param_Import_metadata.txt".format(path), 'w+') as ofstream:
        ofstream.write("# Params used in Gemini_Averaged_MAP_script_CJH.py")
        for param, val in params.items():
            ofstream.write("\n{}:\t{}".format(param, val))
