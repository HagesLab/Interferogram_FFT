# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 11:37:33 2020

@author: Chuck
"""

# New version: FFT fixed
from interferogram_functions import prep_interferogram,  FFT_intr, import_INTR
import matplotlib.pyplot as plt
import os
import numpy as np
from numpy import savetxt

path = r"E:\GEMENI DAQ\NIREOS Complete Example V12_MCS_TimeHarp_32bit Folder\Measurement\20220510/100111"

# path = r"E:\GEMENI DAQ\NIREOS Complete Example V12_MCS_TimeHarp_32bit Folder\Measurement\20211105/115108"
#path = r"F:\PL\Tao\20211012\125809"
save_params = True          #Use this to create a txt file that can be imported into the "..._MP_script" and export Plots
save_PL = True             # Save a .csv of wavelength/PL datasets - one PL per apodization

start_wave =500           #For Plotting - keep in mind the LP filter value
end_wave = 800             #For Plotting
pltzoomstate = False        #Zoom in around the zero position in interferogram to better observe oscillations
pltzoomrange = [-.25,.25]   #Range to zoom in on if pltzoomstate=True

apodization_width=[0.5]     #Bounds (negative to positive) outside of which the data = 0, should be a list. Use many values in the list to compare Apod widths
apod_type="BH"              #Function to use for apodization: "None" "Gauss" "Triangle" "Boxcar" or "BH" (Default)
resample = True             #Enhance resolution by cubic interpolation
resample_factor=4           #Factor to increase data points by
shift= False                #Shift max value to be at 0 mm position - not sure it matters
pad_test = True             #Pad the data with zeros to enhance FFT resolution
padfactor = 16              #Total points will be filled with zeros untill 2**padfactor  -> 2**15 = 32k (default), 2**16=65k
baseline_sub_state = False  #Perform IModPoly baseline subtraction if not a linear baseline (poly = 1) - doesn't work with MAP FFT yet
mean_sub = True             #Shift the average value of the interferogram to be zero
plots = True               #Deactivate plots from FFT and prep - useful if using more than one apod width to compare.

pos_data, intr_data = import_INTR(path)

wave_list, FFT_intr_trim_list = [], []
for i in range(len(apodization_width)):
    
    preFFT_pos, preFFT_data, shiftfactor, baseline_fit = prep_interferogram(pos_data,intr_data,apodization_width[i],apod_type=apod_type,resample=resample,resample_factor=resample_factor,shift=shift,pad_test=pad_test,padfactor=padfactor,mean_sub=mean_sub,plots=plots,pltzoom=pltzoomstate,zoom_range=pltzoomrange,baseline_sub_state=baseline_sub_state)
    wave, FFT_intr_trim = FFT_intr(preFFT_pos,preFFT_data,plots=True,scale="linear",correct=False)
    wave_list.append(wave)
    FFT_intr_trim_list.append(FFT_intr_trim)

#Plot Full PL
plt.figure(5, dpi=120)
plt.ylabel("Counts / a.u.")
plt.xlabel("Wavelength / nm")
plt.title("Average PL")
plt.grid()
#plt.xticks(np.arange(start_wave,end_wave+1,25))
for i in range(len(wave_list)):
    plt.plot(wave_list[i],FFT_intr_trim_list[i],label=apod_type+' Apod '+str(apodization_width[i])+' mm')
plt.xlim(start_wave,end_wave)
plt.yscale('linear')

##### The Raman Zone #####
# def wl_to_raman(wl_exc, wl):
#     return (1e7*(wl_exc**-1 - wl**-1))

# laser_exc = 532
# start_wave = 532            #For Plotting - keep in mind the LP filter value
# end_wave = 550             #For Plotting

# start_wave = wl_to_raman(laser_exc, start_wave)
# end_wave = wl_to_raman(laser_exc, end_wave)
# plt.figure(6, dpi=120)
# plt.ylabel("Counts / a.u.")
# plt.xlabel("Raman shift / cm^-1")
# plt.title("Raman?")
# plt.grid()
# #plt.xticks(np.arange(start_wave,end_wave+1,5))
# for i in range(len(wave_list)):
#     plt.plot(wl_to_raman(laser_exc, wave_list[i]),FFT_intr_trim_list[i],label=apod_type+' Apod '+str(apodization_width[i])+' mm')
# plt.xlim(start_wave,end_wave)
# plt.yscale('linear')
# plt.ylim(0, 5e5)
##########################

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