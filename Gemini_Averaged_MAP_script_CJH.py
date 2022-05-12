# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 11:37:33 2020
@author: Chuck
"""

from interferogram_functions import prep_interferogram, prep_map, FFT_intr, FFT_map
from interferogram_io import save_metadata, save_PL, import_MAP
from make_norm_spec import load_spectrum, interp
from scipy.integrate import simpson
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
from numpy import savetxt
import os

path = r"E:\GEMENI DAQ\NIREOS Complete Example V12_MCS_TimeHarp_32bit Folder\Measurement\20220512\122427"
save_params = 1          #Use this to create a txt file that can be imported into the "..._MAP_script"
export_PL = 1            # Save a .csv of wavelength/PL datasets - one PL per apodization
save_TRPL = 1            # Save a .csv of the avereraged Time/PL dataset

transfer_func = True     # Normalize by a transfer function specific to the optical path
BKGsub = True               #Background Subtract - Generally True
bkg_limit = -3              #ns before the TRPL peak to average the background data up to - see plot
start_wave = 550           #For Plotting
end_wave = 800             #For Plotting
pltzoomstate = False        #Zoom in around the zero position in interferogram to better observe oscillations
pltzoomrange = [-.25,.25]   #Range to zoom in on if pltzoomstate=True

apodization_width=[0.5]     #Bounds (negative to positive) outside of which the data = 0, should be a list. Use many values in the list to compare Apod widths. It will only save the first value in metadata for MAP_script!
apod_type="BH"              #Function to use for apodization: "None" "Gauss" "Triangle" "Boxcar" or "BH" (Default)
resample = True             #Enhance resolution by cubic interpolation
resample_factor=4           #Factor to increase data points by
shift= False                #Shift max value to be at 0 mm position - not sure it matters
pad_test = True             #Pad the data with zeros to enhance FFT resolution
padfactor = 15              #Total points will be filled with zeros untill 2**padfactor  -> 2**15 = 32k (default), 2**16=65k
baseline_sub_state = False   #Perform IModPoly baseline subtraction if not a linear baseline (poly = 1)
mean_sub = True             #Shift the average value of the interferogram to be zero
plots = True               #Deactivate plots from FFT and prep - useful if using more than one apod width to compare.

TRPLmin_OM = 1e-4           #How many orders of magnitude to plot down in y-scale for TRPL curve

# =============================================================================
# Program
# =============================================================================

exper_ID = os.path.split(path)[-1]
pos_data, time_data, map_data = import_MAP(path)

#Background Subtract
t_max = time_data[np.array(np.where(np.mean(map_data,axis=0)==np.max(np.mean(map_data,axis=0)))[0],dtype="int")]
time_data=time_data-t_max
BKGrange = np.array([time_data[0],bkg_limit],dtype='float')  #ns
if BKGsub:
    index = [(np.abs(time_data-np.min(BKGrange))).argmin(),(np.abs(time_data-np.max(BKGrange))).argmin()]
    BKGval = np.mean(map_data[:,np.min(index):np.max(index)],axis=1)
    map_data = map_data - np.reshape(BKGval, (len(BKGval), 1))

#Plot Wavelength Averaged decay to determine time range for background subtraction
plt.figure(3, dpi=120)
plt.plot(time_data,np.mean(map_data,axis=0))
plt.axvspan(np.min(BKGrange),np.max(BKGrange),facecolor='r',alpha=0.2)
plt.xlim(right=2)
#plt.xlabel('Time / ns')
#plt.ylabel('Counts / a.u.')
plt.yscale('log')

#New Time Averaged data from MAP
# index=(np.abs(time_data)).argmin()
# AVG_map_data = map_data[:,index]

AVG_map_data = np.sum(map_data,axis=1)

if transfer_func:
    norm_fname = "cuvet_norm_0.txt"
    norm_waves, norm = load_spectrum(norm_fname)
    
    # Truncate
    norm_waves, norm = interp(norm_waves, norm, start_wave, end_wave, 1)

wave_list, FFT_intr_trim_list = [], []
for i in range(len(apodization_width)):
    preFFT_pos, preFFT_data, shiftfactor, baseline_fit = prep_interferogram(pos_data,AVG_map_data,apodization_width[i],apod_type=apod_type,resample=resample,resample_factor=resample_factor,shift=shift,pad_test=pad_test,padfactor=padfactor,mean_sub=mean_sub,plots=plots,pltzoom=pltzoomstate,zoom_range=pltzoomrange,baseline_sub_state=baseline_sub_state)
    wave, FFT_intr_trim = FFT_intr(preFFT_pos,preFFT_data,plots=True,scale="linear",correct=False)
    
    if transfer_func:
        wave, FFT_intr_trim = interp(wave, FFT_intr_trim, start_wave, end_wave, 1)
        FFT_intr_trim /= norm
        
        wave_list.append(norm_waves)
    else:
        wave_list.append(wave)
        
    FFT_intr_trim_list.append(FFT_intr_trim)

if transfer_func:
    print("transfer")
    preFFT_pos, preFFT_map = prep_map(pos_data,map_data,apodization_width[0],apod_type,resample,resample_factor,shift,pad_test,padfactor,mean_sub)
    FFT_wave, FFT_map = FFT_map(preFFT_pos, preFFT_map)
    FFT_wave=FFT_wave[::-1]
    FFT_map = np.fliplr(np.array(FFT_map,dtype="float").T)
    
    FFT_wave, FFT_map = interp(FFT_wave, FFT_map.T, start_wave, end_wave, 1)

    
    FFT_map = (FFT_map.T / norm).T

    integralTRPL = simpson(FFT_map, x=FFT_wave, axis=0)
    if BKGsub:
        index = [(np.abs(time_data-np.min(BKGrange))).argmin(),(np.abs(time_data-np.max(BKGrange))).argmin()]
        BKGval = np.mean(integralTRPL[np.min(index):np.max(index)])
        integralTRPL = integralTRPL - BKGval
        
        
     
    
else:
    integralTRPL = np.sum(map_data,axis=0)

#Plot Full TRPL
plt.figure(4, dpi=120)
plt.plot(time_data,integralTRPL)
plt.ylim(np.max(integralTRPL)*TRPLmin_OM,2*np.max(integralTRPL))
plt.xlabel('Time / ns')
plt.ylabel('Counts / a.u.')
plt.title("Integral TRPL")
plt.yscale('log')
if save_params:
    PLname = path + "\\" + os.path.split(path)[-1] + '_TRPLPlot.png'
    plt.savefig(PLname)

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

if export_PL:
    PL_fname = os.path.join(path, '{}_PLdata.csv'.format(exper_ID))
    save_PL(PL_fname, wave_list, apodization_width, FFT_intr_trim_list)

if save_TRPL:
    header = []
    outputfilename = path + "\\" + os.path.split(path)[-1] + '_TRPLdata.csv'
    data = np.empty((len(time_data), 2))
    data[:, 0] = time_data
    data[:, 1] = integralTRPL.flatten()
    header.append("Time [ns]")
    header.append("PL [counts]")
    np.savetxt(outputfilename, data, delimiter=',', header=",".join(header))

# Best to turn this on only when you have found the desired params
if save_params:
    outputfilename_meta = os.path.join(path, '{}_FFTmetadata.txt'.format(exper_ID))
    params = {"apod_width":apodization_width[0], "apod_type":apod_type, "do_resample":resample, "resample_factor":resample_factor,
              "do_shift":shift, "do_padzeros":pad_test, "pad_factor":padfactor, "do_baseline_sub":baseline_sub_state, "do_mean_sub":mean_sub,"shift_factor":shiftfactor,"background_subtract":BKGsub,"background_range_low":np.min(BKGrange), "background_range_high":np.max(BKGrange)}
    
    save_metadata(outputfilename_meta, params, from_="Averaged MAP")
    if baseline_sub_state:
        outputfilename_baseline = os.path.join(path, '{}_BaselineFit.txt'.format(exper_ID))
        savetxt(outputfilename_baseline,baseline_fit,delimiter='  ')
