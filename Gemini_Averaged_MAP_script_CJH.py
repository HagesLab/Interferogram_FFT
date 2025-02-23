# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 11:37:33 2020
@author: Chuck
"""

from interferogram_functions import prep_interferogram, prep_map, FFT_intr, FFT_map, where_closest
from interferogram_io import save_metadata, save_PL, save_TRPL, import_MAP
from interferogram_vis import plot_PL_spectrum, plot_TRPL_decay
from make_norm_spec import load_spectrum, interp
from scipy.integrate import simpson
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
import math
from numpy import savetxt
import os

path = r"E:\GEMENI DAQ\NIREOS Complete Example V12_MCS_TimeHarp_32bit Folder\Measurement\20230110\194719"

pre2023 = 0 # Whether this data is taken before 2023. Uses an older calibration file if so.

save_params = 1          #Use this to create a txt file that can be imported into the "..._MAP_script"
export_PL = 1            # Save a .csv of wavelength/PL datasets - one PL per apodization
export_TRPL = 1            # Save a .csv of the avereraged Time/PL dataset

transfer_func = 0     # Normalize by a transfer function specific to the optical path
BKGsub = True               #Background Subtract - Generally True
bkg_limit = -3              #ns before the TRPL peak to average the background data up to - see plot
start_wave = 500           #For Plotting
end_wave = 600             #For Plotting
pltzoomstate = False        #Zoom in around the zero position in interferogram to better observe oscillations
pltzoomrange = [-.25,.25]   #Range to zoom in on if pltzoomstate=True

apodization_width=[1]     #Bounds (negative to positive) outside of which the data = 0, should be a list. Use many values in the list to compare Apod widths. It will only save the first value in metadata for MAP_script!
apod_type="BH"              #Function to use for apodization: "None" "Gauss" "Triangle" "Boxcar" or "BH" (Default)
resample = True             #Enhance resolution by cubic interpolation
resample_factor=4           #Factor to increase data points by
shift= False                #Shift max value to be at 0 mm position - not sure it matters
pad_test = True             #Pad the data with zeros to enhance FFT resolution
padfactor = 15              #Total points will be filled with zeros untill 2**padfactor  -> 2**15 = 32k (default), 2**16=65k
baseline_sub_state = False   #Perform IModPoly baseline subtraction if not a linear baseline (poly = 1)
mean_sub = True             #Shift the average value of the interferogram to be zero
plots = True               #Deactivate plots from FFT and prep - useful if using more than one apod width to compare.

TRPLmin_OM = 1e-6           #How many orders of magnitude to plot down in y-scale for TRPL curve
start_time = -10
end_time = 20
# =============================================================================
# Program
# =============================================================================

exper_ID = os.path.split(path)[-1]
pos_data, time_data, map_data = import_MAP(path)

# Shift t=0 onto max intensity
t_max = time_data[np.array(np.where(np.mean(map_data,axis=0)==np.max(np.mean(map_data,axis=0)))[0],dtype="int")]
time_data -= t_max

BKGrange = np.array([time_data[0],bkg_limit],dtype='float')  #ns

# Validation #
if BKGrange[0] >= BKGrange[1]:
    raise ValueError("Invalid BKGrange {} to {}".format(BKGrange[0], BKGrange[1]))
    
if np.any(np.diff(time_data) <= 0):
    raise ValueError("time_data is not monotonically ascending")
##############

if BKGsub:
    index = where_closest(time_data, BKGrange)
    BKGval = np.mean(map_data[:,index[0]:index[1]],axis=1)
    map_data = map_data - np.reshape(BKGval, (len(BKGval), 1))

#Plot Wavelength Averaged decay to determine time range for background subtraction
plt.figure(3, dpi=120)
plt.plot(time_data,np.mean(map_data,axis=0))
plt.axvspan(BKGrange[0],BKGrange[1],facecolor='r',alpha=0.2)
plt.xlim(right=2)
plt.xlabel('Time / ns')
plt.ylabel('Counts / a.u.')
plt.yscale('log')

#New Time Averaged data from MAP
# index=(np.abs(time_data)).argmin()
# AVG_map_data = map_data[:,index]

AVG_map_data = np.sum(map_data,axis=1)

if transfer_func:
    norm_fname = os.path.join(r"\\10.227.108.119\share\Data\PL\Transfer functions\Micro - Intf - Vis", 
                              "microPL_vis_intf.txt")
    norm_waves, norm = load_spectrum(norm_fname)
    norm_waves, norm = interp(norm_waves, norm, start_wave, end_wave, 1)

wave_list, FFT_intr_trim_list = [], []
for i in range(len(apodization_width)):
    preFFT_pos, preFFT_data, shiftfactor, baseline_fit = prep_interferogram(pos_data,AVG_map_data,apodization_width[i],apod_type=apod_type,resample=resample,resample_factor=resample_factor,shift=shift,
                                                                            pad_test=pad_test,padfactor=padfactor,mean_sub=mean_sub,plots=plots,pltzoom=pltzoomstate,zoom_range=pltzoomrange,
                                                                            baseline_sub_state=baseline_sub_state,
                                                                            pre2023=pre2023)
    wave, FFT_intr_trim = FFT_intr(preFFT_pos,preFFT_data,plots=True,scale="linear",correct=False)
    
    if transfer_func:
        wave, FFT_intr_trim = interp(wave, FFT_intr_trim, start_wave, end_wave, 1)
        FFT_intr_trim /= norm
        
        wave_list.append(norm_waves)
    else:
        wave, FFT_intr_trim = interp(wave, FFT_intr_trim, math.ceil(np.amin(wave)), math.floor(np.amax(wave)), 1)
        wave_list.append(wave)
        
    FFT_intr_trim_list.append(FFT_intr_trim)

# if transfer_func:
#     print("transfer")
#     preFFT_pos, preFFT_map = prep_map(pos_data,map_data,apodization_width[0],apod_type,resample,resample_factor,shift,pad_test,padfactor,mean_sub)
#     FFT_wave, FFT_map = FFT_map(preFFT_pos, preFFT_map)
#     FFT_wave=FFT_wave[::-1]
#     FFT_map = np.fliplr(np.array(FFT_map,dtype="float").T)
    
#     FFT_wave, FFT_map = interp(FFT_wave, FFT_map.T, start_wave, end_wave, 1)

    
#     FFT_map = (FFT_map.T / norm).T

#     integralTRPL = simpson(FFT_map, x=FFT_wave, axis=0)
#     if BKGsub:
#         index = where_closest(time_data, BKGrange)
#         BKGval = np.mean(integralTRPL[index[0]:index[1]])
#         integralTRPL = integralTRPL - BKGval
        
            
# else:
integralTRPL = np.sum(map_data,axis=0)

#Plot Full TRPL
PLname = os.path.join(path, '{}_TRPLPlot.png'.format(exper_ID))
plot_TRPL_decay(time_data, [integralTRPL], TRPLmin_OM, export=PLname, labels=[None],
                start_time = start_time, end_time = end_time)

#Plot Full PL
PLname = os.path.join(path, '{}_PLPlot.png'.format(exper_ID))
labels = ['{} Apod {} mm'.format(apod_type, apod) for apod in apodization_width]
plot_PL_spectrum(wave_list, FFT_intr_trim_list, labels, start_wave, end_wave, export=PLname)

if export_PL:
    PL_fname = os.path.join(path, '{}_PLdata.csv'.format(exper_ID))
    save_PL(PL_fname, wave_list, apodization_width, FFT_intr_trim_list)

if export_TRPL:
    trpl_fname = os.path.join(path, '{}_TRPLdata.csv'.format(exper_ID))
    save_TRPL(trpl_fname, time_data, integralTRPL)

# Best to turn this on only when you have found the desired params
if save_params:
    outputfilename_meta = os.path.join(path, '{}_FFTmetadata.txt'.format(exper_ID))
    params = {"apod_width":apodization_width[0], "apod_type":apod_type, "do_resample":resample, "resample_factor":resample_factor,
              "do_shift":shift, "do_padzeros":pad_test, "pad_factor":padfactor, "do_baseline_sub":baseline_sub_state, "do_mean_sub":mean_sub,"shift_factor":shiftfactor,"background_subtract":BKGsub,"background_range_low":np.min(BKGrange), "background_range_high":np.max(BKGrange)}
    
    save_metadata(outputfilename_meta, params, from_="Averaged MAP")
    if baseline_sub_state:
        outputfilename_baseline = os.path.join(path, '{}_BaselineFit.txt'.format(exper_ID))
        savetxt(outputfilename_baseline,baseline_fit,delimiter='  ')
