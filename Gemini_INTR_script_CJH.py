# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 11:37:33 2020

@author: Chuck
"""

from interferogram_functions import prep_interferogram,  FFT_intr, import_INTR
import matplotlib.pyplot as plt
import os
import numpy as np
from numpy import savetxt

path = r"C:/Users/cfai2/Documents/src/Interferogram_FFT/20210615/132214"
save_params = True          #Use this to create a txt file that can be imported into the "..._MAP_script" and export Plots
save_PL = True             # Save a .csv of wavelength/PL datasets - one PL per apodization

start_wave = 535            #For Plotting - keep in mind the LP filter value
end_wave = 1100             #For Plotting
pltzoomstate = False        #Zoom in around the zero position in interferogram to better observe oscillations
pltzoomrange = [-.25,.25]   #Range to zoom in on if pltzoomstate=True

apodization_width=[1, 0.7,0.5,0.3]     #Bounds (negative to positive) outside of which the data = 0, should be a list. Use many values in the list to compare Apod widths
apod_type="BH"              #Function to use for apodization: "None" "Gauss" "Triangle" "Boxcar" or "BH" (Default)
resample = True             #Enhance resolution by cubic interpolation
resample_factor=4           #Factor to increase data points by
shift= False                #Shift max value to be at 0 mm position - not sure it matters
pad_test = True             #Pad the data with zeros to enhance FFT resolution
padfactor = 14              #Total points will be filled with zeros untill 2**padfactor  -> 2**15 = 32k (default), 2**16=65k
baseline_sub_state = False  #Perform IModPoly baseline subtraction if not a linear baseline (poly = 1) - doesn't work with MAP FFT yet
mean_sub = True             #Shift the average value of the interferogram to be zero
plots = True               #Deactivate plots from FFT and prep - useful if using more than one apod width to compare.

pos_data, intr_data = import_INTR(path)

wave_list, FFT_intr_trim_list = [], []
for i in range(len(apodization_width)):
    preFFT_pos, preFFT_data, shiftfactor, baseline_fit = prep_interferogram(pos_data,intr_data,apodization_width[i],apod_type=apod_type,resample=resample,resample_factor=resample_factor,shift=shift,pad_test=pad_test,padfactor=padfactor,mean_sub=mean_sub,plots=plots,pltzoom=pltzoomstate,zoom_range=pltzoomrange,baseline_sub_state=baseline_sub_state)
    wave, FFT_intr_trim = FFT_intr(preFFT_pos,preFFT_data,plots=True,scale="linear",correct=True)
    wave_list.append(wave)
    FFT_intr_trim_list.append(FFT_intr_trim)

#Plot Full PL
plt.figure(5, dpi=120)
plt.ylabel("Counts / a.u.")
plt.xlabel("Wavelength / nm")
plt.title("Average PL")
for i in range(len(wave_list)):
    plt.plot(wave_list[i],FFT_intr_trim_list[i],label=apod_type+' Apod '+str(apodization_width[i])+' mm')
plt.xlim(start_wave,end_wave)
plt.yscale('linear')
if save_params:
    PLname = path + "\\" + os.path.split(path)[-1] + '_PLPlot.png'
    plt.savefig(PLname)
if len(apodization_width) > 1:
    plt.legend()

if save_PL:
    header = []
    outputfilename = path + "\\" + os.path.split(path)[-1] + '_PLdata.csv'
    data = np.empty((len(wave_list[0]), len(apodization_width) * 2))
    for i, apod in enumerate(apodization_width):
        data[:, 2*i] = wave_list[i]
        data[:, 2*i+1] = FFT_intr_trim_list[i].flatten()
        header.append("Wavelength [nm] apod={}".format(apod))
        header.append("PL [counts] apod={}".format(apod))
    np.savetxt(outputfilename, data, delimiter=',', header=",".join(header))

# Best to turn this on only when you have found the desired params
if save_params:
    outputfilename_meta = path + "\\" + os.path.split(path)[-1] + '_FFTmetadata.txt'
    params = {"apod_width":apodization_width, "apod_type":apod_type, "do_resample":resample, "resample_factor":resample_factor,
              "do_shift":shift, "do_padzeros":pad_test, "pad_factor":padfactor, "do_baseline_sub":baseline_sub_state, "do_mean_sub":mean_sub,"shift_factor":shiftfactor}
    with open(outputfilename_meta.format(path), 'w+') as ofstream:
        ofstream.write("# Params used in Gemini_INTR_script_CJH.py")
        for param, val in params.items():
            ofstream.write("\n{}:\t{}".format(param, val))

    if baseline_sub_state:
        outputfilename_baseline = path + "\\" + os.path.split(path)[-1] + '_BaselineFit.txt'
        savetxt(outputfilename_baseline,baseline_fit,delimiter='  ')