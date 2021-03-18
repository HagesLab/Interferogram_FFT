# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 11:37:33 2020

@author: Chuck
"""

from interferogram_functions import prep_interferogram, FFT_intr, import_MAP
import matplotlib.pyplot as plt
import numpy as np

path = r"C:\Users\Chuck\Dropbox (UFL)\UF\TRPL Computer\Calvin\142824"  
BKGsub = True               #Background Subtract - Generally True
bkg_limit = -5              #ns before the TRPL peak to average the background data up to - see plot
save_params = True          #Use this to create a txt file that can be imported into the "..._MAP_script"
start_wave = 1000           #For Plotting
end_wave = 1700             #For Plotting
pltzoomstate = False        #Zoom in around the zero position in interferogram to better observe oscillations

apodization_width=0.95      #Bounds (negative to positive) outside of which the data = 0
apod_type="BH"              #Function to use for apodization: "None" "Gauss" "Triangle" "Boxcar" or "BH" (Default)
resample="True"             #Enhance resolution by cubic interpolation
resample_factor=4           #Factor to increase data points by
shift="False"               #Shift max value to be at 0 mm position - not sure it matters
pad_test="True"             #Pad the data with zeros to enhance FFT resolution
padfactor = 15              #Total points will be filled with zeros untill 2**padfactor  -> 2**15 = 32k (default), 2**16=65k 
mean_sub = "True"           #Shift the average value of the interferogram to be zero

TRPLmin_OM = 1e-7           #How many orders of magnitude to plot down in y-scale for TRPL curve

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
plt.figure(0, dpi=120)
plt.plot(time_data,np.mean(map_data,axis=0))
plt.axvspan(np.min(BKGrange),np.max(BKGrange),facecolor='r',alpha=0.2)
plt.xlim(right=2)
#plt.xlabel('Time / ns')
#plt.ylabel('Counts / a.u.')
plt.yscale('log')

#New Time Averaged data from MAP
index=(np.abs(time_data)).argmin()
#AVG_map_data = np.mean(map_data,axis=1)
AVG_map_data = map_data[:,index]


preFFT_pos, preFFT_data, shiftfactor = prep_interferogram(pos_data,AVG_map_data,apodization_width,apod_type=apod_type,resample=resample,resample_factor=resample_factor,shift=shift,pad_test=pad_test,padfactor=padfactor,mean_sub=mean_sub,plots="True",pltzoom=pltzoomstate)
wave, FFT_intr_trim = FFT_intr(preFFT_pos,preFFT_data,plots="True",scale="linear",correct="True")

#Plot Full TRPL
plt.figure(4, dpi=120)
integralTRPL = np.sum(map_data,axis=0)
plt.figure(3, dpi=120)
plt.plot(time_data,integralTRPL)
plt.ylim(np.max(integralTRPL)*TRPLmin_OM,2*np.max(integralTRPL))
#plt.xlabel('Time / ns')
#plt.ylabel('Counts / a.u.')
#plt.title("Integral TRPL")
plt.yscale('log')

#Plot Full PL
plt.figure(4, dpi=120)
#plt.ylabel("Counts / a.u.")
#plt.xlabel("Wavelength / nm")
#plt.title("Average PL")
plt.plot(wave,FFT_intr_trim)
plt.xlim(start_wave,end_wave)
plt.yscale('linear')
plt.show()

# Best to turn this on only when you have found the desired params
if save_params:
    params = {"apod_width":apodization_width, "apod_type":apod_type, "do_resample":resample, "resample_factor":resample_factor,
              "do_shift":shift, "do_padzeros":pad_test, "pad_factor":padfactor, "do_mean_sub":mean_sub, "shift_factor":shiftfactor,"background_subtract":BKGsub,"background_range_low":BKGrange[0], "background_range_high":BKGrange[1]}
    with open(r"{}\Param_Import_metadata.txt".format(path), 'w+') as ofstream:
        ofstream.write("# Params used in Gemini_Averaged_MAP_script_CJH.py")
        for param, val in params.items():
            ofstream.write("\n{}:\t{}".format(param, val))
