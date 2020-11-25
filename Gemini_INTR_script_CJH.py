# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 12:20:31 2020

@author: Chuck
"""

from interferogram_functions import prep_interferogram, import_INTR, FFT_intr, import_MAP
import matplotlib.pyplot as plt
import os
import numpy as np

#Apply to "*_INTR.txt" data
path = r"\\ad.ufl.edu\che\Hages-Lab\TRPL Data\Hages\Ferekides\20_11_24\152430"
pos_data, intr_data = import_INTR(path)

apodization_width=1.75
apod_type="BH"    # "None" "Gauss" "Triangle" "Boxcar" "BH"
resample="True"
resample_factor=4
shift="True"
pad_test="True"
padfactor=4
mean_sub = "True"

preFFT_pos, preFFT_data, shiftfactor = prep_interferogram(pos_data,intr_data,apodization_width,apod_type=apod_type,resample=resample,resample_factor=resample_factor,shift=shift,pad_test=pad_test,padfactor=padfactor,mean_sub=mean_sub,plots="True")
wave, FFT_intr_trim = FFT_intr(preFFT_pos,preFFT_data,plots="True",scale="linear",correct="True")

start_wave = 700
end_wave = 1000
plt.figure(3, dpi=120)
plt.plot(wave,FFT_intr_trim)
plt.xlim(start_wave,end_wave)
plt.yscale('linear')
plt.show()



# =============================================================================
# #Apply to single line of MAP data
# pos_data, time_data, map_data = import_MAP(path)
# pos_data = pos_data - shiftfactor    #Taken from shift_factor output in the _INTR analysis script for this data
#
# #ADD MANUAL SHIFT by shift factor - currently done above manually
# #WHAT ABOUT INTENSITY CORRECTION? Ask NIREOS
# shift="False"
#
# #pick time to check
# rangeval = 12  #ns
# index = (np.abs(time_data-np.min(rangeval))).argmin()
#
# preFFT_pos, preFFT_data, shiftfactor = prep_interferogram(pos_data,map_data[:,index],apodization_width,apod_type=apod_type,resample=resample,resample_factor=resample_factor,shift=shift,pad_test=pad_test,padfactor=padfactor,mean_sub=mean_sub,plots="True")
# wave, FFT_intr_trim = FFT_intr(preFFT_pos,preFFT_data,plots="True",scale="linear",correct="True")
#
# plt.figure(4, dpi=120)
# plt.plot(wave,FFT_intr_trim)
# plt.xlim(start_wave,end_wave)
# plt.yscale('linear')
# plt.show()
# =============================================================================
