# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 12:20:31 2020

@author: Chuck
"""

from interferogram_functions import prep_interferogram, import_INTR, FFT_intr
import matplotlib.pyplot as plt

path = r"C:\Users\Chuck\Dropbox (UFL)\Hages Lab Files\TRPL Data\Ruiquan\20201104\153345"
pos_data, intr_data = import_INTR(path)
      
apodization_width=0.5
apod_type="BH"    # "None" "Gauss" "Triangle" "Boxcar" "BH"
resample="True"
resample_factor=2
shift="True"
pad_test="True"
padfactor=4
mean_sub = "True"       
   
preFFT_pos, preFFT_data, shiftfactor = prep_interferogram(pos_data,intr_data,apodization_width,apod_type=apod_type,resample=resample,resample_factor=resample_factor,shift=shift,pad_test=pad_test,padfactor=padfactor,mean_sub=mean_sub,plots="True")
wave, FFT_intr_trim = FFT_intr(preFFT_pos,preFFT_data,plots="True",scale="linear",correct="True")

start_wave = 400
end_wave = 1000
plt.figure(3, dpi=120)
plt.plot(wave,FFT_intr_trim)
plt.xlim(start_wave,end_wave)
plt.yscale('linear')
plt.show()