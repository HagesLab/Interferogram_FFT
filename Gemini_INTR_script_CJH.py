# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 12:20:31 2020

@author: Chuck
"""

from interferogram_functions import prep_interferogram, import_INTR, FFT_intr
import matplotlib.pyplot as plt

#Apply to "*_INTR.txt" data
path = r"C:\Users\Chuck\Desktop\20_11_24\180627"
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

# Best to turn this on only when you have found the desired params
save_params = False

if save_params:
    params = {"apod_width":apodization_width, "apod_type":apod_type, "do_resample":resample, "resample_factor":resample_factor,
              "do_shift":shift, "do_padzeros":pad_test, "pad_factor":padfactor, "do_mean_sub":mean_sub, "shift_factor":shiftfactor}
    with open(r"{}\INTR_metadata.txt".format(path), 'w+') as ofstream:
        ofstream.write("# Params used in Gemini_INTR_script_CJH.py")
        for param, val in params.items():
            ofstream.write("\n{}:\t{}".format(param, val))
            

