# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 21:44:47 2020

@author: Chuck
"""
import numpy as np
import Interf_calibration1D
from scipy.interpolate import interp1d
import Interf_calibration2D
import matplotlib.pyplot as plt
import os
import pandas as pd
from numpy.polynomial import Polynomial
from scipy.optimize import curve_fit

try:
    from BaselineRemoval import BaselineRemoval
except ModuleNotFoundError:
    print("Warning: BaselineRemoval library not found")

def prep_interferogram(pos_data,intr_data,apodization_width,apod_type="BH",mean_sub=True,resample=True,resample_factor=4,shift=False,pad_test=True,padfactor=15,plots=True,pltzoom=False,zoom_range=[-.25,0.25],baseline_sub_state=False):
    intr_data, pos_data = Interf_calibration1D.interf_calibration1D(intr_data,pos_data)
    shiftfactor = 0
    raw_input_intr, raw_input_pos = intr_data, pos_data
    if baseline_sub_state:
        baseObj = BaselineRemoval(intr_data)
        intr_data = baseObj.IModPoly(1)
    if mean_sub:
        intr_data = intr_data-np.mean(intr_data)
    baseline_fit = raw_input_intr - intr_data
    if resample:
        fcubic = interp1d(pos_data, intr_data, kind='cubic')
        pos_data = np.linspace(pos_data[0],pos_data[-1],endpoint=True,num=len(pos_data)*resample_factor)
        intr_data = fcubic(pos_data)
    index_pos = np.argmin(abs(pos_data))
    if shift:
        index_intr = np.argmax(intr_data)
        shiftfactor=(pos_data[index_intr]-pos_data[index_pos])
        pos_data = pos_data-shiftfactor
    raw_intr = intr_data
    raw_pos = pos_data
    if apod_type != "None":
        if apod_type == "Gauss":
            intr_func = np.exp(-(2.24*(pos_data)/apodization_width)**2)
        if apod_type == "BH":
            intr_func = np.zeros(len(pos_data))
            for i in range(len(intr_func)):
                if abs(pos_data[i]) <= apodization_width:
                    A0, A1, A2, A3 = 0.42323, 0.49755, 0.07922, 0.01168
                    intr_func[i] = A0+A1*np.cos(np.pi*pos_data[i]/apodization_width)+A2*np.cos(2*np.pi*pos_data[i]/apodization_width)+A3*np.cos(3*np.pi*pos_data[i]/apodization_width)
                elif abs(pos_data[i]) > apodization_width:
                    intr_func[i] = 0
        if apod_type == "Triangle":
            intr_func = np.zeros(len(pos_data))
            for i in range(len(intr_func)):
                if abs(pos_data[i]) <= apodization_width:
                    intr_func[i] = 1-abs(pos_data[i])/apodization_width
                elif abs(pos_data[i]) > apodization_width:
                    intr_func[i] = 0
        if apod_type == "Boxcar":
            intr_func = np.zeros(len(pos_data))
            for i in range(len(intr_func)):
                if abs(pos_data[i]) <= apodization_width:
                    intr_func[i] = 1
                elif abs(pos_data[i]) > apodization_width:
                    intr_func[i] = 0
        intr_data = intr_data*intr_func
    posdiff=np.diff(pos_data)[0]
    if pad_test:
        padlen = 2**padfactor-len(intr_data)
        if padlen >0:
            intr_data = np.pad(intr_data,(0,int(padlen)),'constant', constant_values=(0))
            pos_data = np.append(pos_data,np.linspace(pos_data[-1]+posdiff,(pos_data[-1]+posdiff)+padlen*posdiff,num=padlen))

    left_axis_intr, right_axis_intr = intr_data[0:index_pos+1], intr_data[index_pos+1:]
    left_axis_pos, right_axis_pos = pos_data[0:index_pos+1], pos_data[index_pos+1:]
    preFFT_data, preFFT_pos = np.concatenate((right_axis_intr,left_axis_intr)), np.concatenate((right_axis_pos,left_axis_pos))

    if plots:
        #Plot shifted data
        plt.figure(0,dpi=120)
        plt.title("Raw Comparison")
        plt.plot(raw_input_pos, raw_input_intr,label="Raw")
        plt.plot(raw_pos, raw_intr,label="Corrected")
        plt.xlabel("Position / mm")
        plt.ylabel("Counts / a.u.")
        if pltzoom:
            plt.xlim(np.min(zoom_range),np.max(zoom_range))
        plt.legend()
        plt.show()

        #Plot interferogram data
        if apod_type != "None":
            plt.figure(1, dpi=120)
            plt.plot(raw_pos,raw_intr/np.max(raw_intr),label="Input Data")
            plt.plot(pos_data[0:len(raw_pos)],(intr_data/np.max(intr_data))[0:len(raw_pos)],label="Apod Data")
            plt.plot(raw_pos,intr_func,label="Apod Fxn")
            plt.title("Apodization")
            if pltzoom:
                plt.xlim(np.min(zoom_range),np.max(zoom_range))
            plt.xlabel('Position / mm')
            plt.ylabel('Counts / a.u.')
            plt.show()

        #Plot Pre_FFT Data set
        plt.figure(2, dpi=120)
        plt.title("Pre-FFT data")
        plt.ylabel("Counts / a.u.")
        plt.plot(preFFT_data)
        plt.show()

    return preFFT_pos, preFFT_data, shiftfactor, baseline_fit

def FFT_intr(preFFT_pos,preFFT_data, plots=False,correct=True,scale="linear"):
    # Treat single TS 1D case as 2D case
    if preFFT_data.ndim == 1:
        preFFT_data = preFFT_data.reshape((len(preFFT_data),1))

    #Do as much preparation outside of the loop over preFFT_data's timesteps as possible
    freq = np.fft.rfftfreq(preFFT_pos.shape[-1],np.diff(preFFT_pos)[0])

    #Import calibration for WL
    items = os.listdir(".")
    for names in items:
        if names.endswith("parameters_cal.txt"):
            filename = names
    ref = pd.read_csv(filename, sep="\t", header=None)
    first_row = (ref.iloc[0])
    second_row = (ref.iloc[1])
    wavelength = first_row.to_numpy(dtype='float64')
    reciprocal = second_row.to_numpy(dtype='float64')
    min_freq = np.min(reciprocal)
    max_freq = np.max(reciprocal)
    # Is freq[] always in ascending order? If so this can be made much simpler with array slicing
    select = np.intersect1d(np.where(freq >= min_freq), np.where(freq <= max_freq))

    fn = interp1d(reciprocal, wavelength, kind="linear")

    # Do FFT on the data
    FFT_intr_full = np.fft.rfft(preFFT_data, axis=0)

    if correct:
        #Phase Corrected
        FFT_real_full_raw = FFT_intr_full.real
        FFT_imag_full_raw = FFT_intr_full.imag
        phase_angle = np.arctan(FFT_imag_full_raw / FFT_real_full_raw)
        FFT_phase_corrected_full = FFT_intr_full*np.exp(-1j*phase_angle)
        FFT_final_full = FFT_phase_corrected_full.real
        
    else:
        FFT_real_full_raw = FFT_intr_full.real
        FFT_imag_full_raw = FFT_intr_full.imag
        FFT_final_full = np.sqrt(FFT_real_full_raw**2 + FFT_imag_full_raw**2)
        
    
    # Trim according to calibrations
    freq_trim = freq[select]
    FFT_intr_trim_full = FFT_final_full[select]

    #Compute WL
    wave = fn(freq_trim)

    #Plot
    if plots:
        if correct:
            title = "Phase Corrected FFT"
        else:
            title = "FFT"
        FFT_real_raw = FFT_real_full_raw[select]
        FFT_imag_raw = FFT_imag_full_raw[select]
        plt.figure(0, dpi=120)
        plt.plot(wave,FFT_intr_trim_full,"--",label="Full")
        plt.plot(wave,FFT_real_raw,label="Real")
        plt.plot(wave,FFT_imag_raw,label="Imag")
        plt.yscale(scale)
        plt.title(title)
        plt.xlabel("Wavelength (nm)")
        #plt.xlim((475,700))
        plt.grid()
        #plt.xticks(np.arange(475,701,25))
        plt.legend()

    return wave, FFT_intr_trim_full.reshape((1,len(FFT_intr_trim_full)))[0]

#Need to add baseline subtraction to MAP
def prep_map(pos_data,map_data,apodization_width,apod_type="BH",resample=True,resample_factor=2,shift=False,pad_test=True,padfactor=4,mean_sub=True):
    map_data, pos_data = Interf_calibration2D.interf_calibration2D(map_data,pos_data)
    prep_build=[]
    pos_data_raw = pos_data
    if mean_sub:
        map_data = map_data - np.mean(map_data, axis=0)

    if resample:
        fcubic = interp1d(pos_data_raw, map_data, kind='cubic', axis=0)
        pos_data_raw = np.linspace(pos_data_raw[0],pos_data_raw[-1],endpoint=True,num=len(pos_data_raw)*resample_factor)
        map_data = fcubic(pos_data_raw)

    for i in range(map_data.shape[1]): # For each timestep
        intr_data = map_data[:,i] # Take a column
        pos_data = pos_data_raw

        index_pos = np.argmin(abs(pos_data))
        if shift:
            index_intr = np.argmax(intr_data)
            #index_intr2 = np.argmax(intr_data, axis=1)
            shiftfactor=(pos_data[index_intr]-pos_data[index_pos])
            #shiftfactor2 = np.array([pos_data[i] - pos_data[index_pos] for i in index_intr2])
            #pos_data2
            pos_data = pos_data-shiftfactor
        if apod_type != "None":
            if apod_type == "Gauss":
                intr_func = np.exp(-(2.24*(pos_data)/apodization_width)**2)
            if apod_type == "BH":
                A0, A1, A2, A3 = 0.42323, 0.49755, 0.07922, 0.01168
                intr_func = np.where(abs(pos_data) <= apodization_width,
                                      A0+A1*np.cos(np.pi*pos_data/apodization_width)+A2*np.cos(2*np.pi*pos_data/apodization_width)+A3*np.cos(3*np.pi*pos_data/apodization_width),
                                      0)
            if apod_type == "Triangle":
                intr_func = np.where(abs(pos_data) <= apodization_width,
                                      1-abs(pos_data)/apodization_width,
                                      0)
            if apod_type == "Boxcar":
                intr_func = np.where(abs(pos_data) <= apodization_width, 1, 0)

            intr_data = intr_data*intr_func
        posdiff=np.diff(pos_data)[0]
        if pad_test:
            padlen = int(padfactor*(2**np.ceil(np.log(len(intr_data))/np.log(2))))-len(intr_data)
            intr_data = np.pad(intr_data,(0,int(padlen)),'constant', constant_values=(0))
            pos_data = np.append(pos_data,np.linspace(pos_data[-1]+posdiff,(pos_data[-1]+posdiff)+padlen*posdiff,num=padlen))

        index_pos = np.argmin(abs(pos_data))
        left_axis_intr, right_axis_intr = intr_data[0:index_pos+1], intr_data[index_pos+1:]
        left_axis_pos, right_axis_pos = pos_data[0:index_pos+1], pos_data[index_pos+1:]
        preFFT_data, preFFT_pos = np.concatenate((right_axis_intr,left_axis_intr)), np.concatenate((right_axis_pos,left_axis_pos))
        prep_build.append(preFFT_data)

    prep_build = np.array(prep_build, dtype='float').T
    return preFFT_pos, prep_build

def FFT_map(FFT_pos,FFT_data, plots=False,correct=False,scale="linear"):
    # Treat single TS 1D case as 2D case
    if FFT_data.ndim == 1:
        FFT_data = FFT_data.reshape((len(FFT_data),1))

    #Do as much preparation outside of the loop over preFFT_data's timesteps as possible
    freq = np.fft.rfftfreq(FFT_pos.shape[-1],np.diff(FFT_pos)[0])

    #Import calibration for WL
    items = os.listdir(".")
    for names in items:
        if names.endswith("parameters_cal.txt"):
            filename = names
    ref = pd.read_csv(filename, sep="\t", header=None)
    first_row = (ref.iloc[0])
    second_row = (ref.iloc[1])
    wavelength = first_row.to_numpy(dtype='float64')
    reciprocal = second_row.to_numpy(dtype='float64')
    min_freq = np.min(reciprocal)
    max_freq = np.max(reciprocal)
    # Is freq[] always in ascending order? If so this can be made much simpler with array slicing
    select = np.intersect1d(np.where(freq >= min_freq), np.where(freq <= max_freq))

    fn = interp1d(reciprocal, wavelength, kind="linear")

    # Do FFT on the data
    FFT_data = np.fft.rfft(FFT_data, axis=0)

    if correct:
        #Phase Corrected
        # FIXME: this should match FFT_intr
        # FFT_real_full_raw = FFT_data.real
        # FFT_imag_full_raw = FFT_data.imag
        # FFT_final_full_raw = np.sqrt(FFT_real_full_raw**2 + FFT_imag_full_raw**2)
        FFT_real_full = FFT_data.real*np.cos(np.angle(FFT_data))
        FFT_imag_full = FFT_data.imag*np.sin(np.angle(FFT_data))
        FFT_final_full = np.sqrt(FFT_real_full**2 + FFT_imag_full**2)
    else:
        #FFT_data = np.sqrt(FFT_data.real**2 + FFT_data.imag**2)
        FFT_data = np.abs(FFT_data) # equivalent to real^2+imag^2
        

    # Trim according to calibrations
    freq_trim = freq[select]
    #FFT_intr_trim_full = FFT_data[select]

    #Compute WL
    wave = fn(freq_trim)

    #Plot
    if plots:
        if correct:
            FFT_real_raw = FFT_real_full_raw[select]
            FFT_imag_raw = FFT_imag_full_raw[select]
            FFT_intr_raw = FFT_final_full_raw[select]
            plt.figure(0, dpi=120)
            plt.plot(wave,FFT_intr_raw,"--",label="Full")
            plt.plot(wave,FFT_real_raw,label="Real")
            plt.plot(wave,FFT_imag_raw,label="Imag")
            plt.yscale(scale)
            plt.title("Raw FFT")
            plt.xlabel("Wavelength (nm)")
            plt.legend()

        # FFT_real = FFT_data.real[select]
        # FFT_imag = FFT_data.imag[select]

        plt.figure(1, dpi=120)
        plt.plot(wave,FFT_data[select],"--",label="Full")
        plt.plot(wave,FFT_data.real[select],label="Real")
        plt.plot(wave,FFT_data.imag[select],label="Imag")
        plt.yscale(scale)
        if correct:
            plt.title("Phase Corrected FFT")
        else:
            plt.title("Raw FFT")
        plt.xlabel("Wavelength (nm)")
        plt.legend()

    return wave, FFT_data[select]



#Fit TRPL
def Fit_1exp(TRPL_data,time_data,fitrange):

    def Exp1(time,A,tau):
        return -time/tau + np.log(A)

    #trim-data
    low_index, high_index = (np.abs(time_data-np.min(fitrange))).argmin() , (np.abs(time_data-np.max(fitrange))).argmin()
    TRPL_fit = TRPL_data[low_index:high_index]
    time_fit = time_data[low_index:high_index]
    
    p = Polynomial.fit(time_fit, np.log(np.abs(TRPL_fit)), 1)
    popt = p.convert().coef
    popt[1] = -1 / popt[1]
    popt[0] = np.exp(popt[0])


    TRPL_out = np.exp(Exp1(time_fit,*popt))
    label =  r'$\tau:\ $' + np.array2string(popt[1], precision=2, separator=',', suppress_small=True) + ' ns'

    return TRPL_out, time_fit, label

def Fit_2exp(TRPL_data,time_data,fitrange):

    def Exp2(time,A,tau,B, tau2):
        # Log scale prevents overweighting of higher magnitude points
        return np.log(A * np.exp(-time/tau) + B * np.exp(-time/tau2))

    #trim-data
    low_index, high_index = (np.abs(time_data-np.min(fitrange))).argmin() , (np.abs(time_data-np.max(fitrange))).argmin()
    TRPL_fit = TRPL_data[low_index:high_index]
    time_fit = time_data[low_index:high_index]

    popt, pcov = curve_fit(Exp2,time_fit,np.log(np.abs(TRPL_fit)))
    perr = np.sqrt(np.diag(pcov))

    TRPL_out = np.exp(Exp2(time_fit,*popt))
    tau1 = popt[3]
    tau2 = popt[1]
    C1 = popt[2]
    C2 = popt[0]
    
    if tau1 > tau2:
        # Swap so tau1 always the smallest
        tau2, tau1 = tau1, tau2
        C2, C1 = C1, C2
    
    taustat = (C1*tau1+C2*tau2) / (C1+C2)
    
    label =  r'$\tau_1:\ $ {:.2f} ns, $\tau_2:\ $ {:.2f} ns, '.format(tau1, tau2)
    label += r"$\tau_{stat}:\ $" + "{:.2f} ns".format(taustat)

    return TRPL_out, time_fit, label, popt, perr

def RR(time,log_B, log_dN0):
    # See Ahrenkiel book pp. 64
    # Log scale prevents overweighting of higher magnitude points
    return log_B + 2*log_dN0 - 2*np.log(1+np.exp(log_B + log_dN0)*time) + np.log(2000e-7)

def Fit_RR(TRPL_data,time_data,fitrange):
    # time_data in [s] required!

    #trim-data
    low_index, high_index = (np.abs(time_data-np.min(fitrange))).argmin() , (np.abs(time_data-np.max(fitrange))).argmin()
    TRPL_fit = TRPL_data[low_index:high_index]
    time_fit = time_data[low_index:high_index]

    popt, pcov = curve_fit(RR,time_fit,np.log(np.abs(TRPL_fit)))
    perr = np.sqrt(np.diag(pcov))


    dN0 = np.exp(popt[1])
    B = np.exp(popt[0])
    TRPL_out = np.exp(RR(time_fit,*popt))
    
    
    label =  r'$B:\ $ {:.2e}'.format(B) + r' $\mathrm{cm^3 s^{-1}}$' + ' $\Delta N_0:\ $ {:.2e}'.format(dN0) + r' $\mathrm{cm^{-3}}$'
    
    return TRPL_out, time_fit, label, popt, perr

def where_closest(x : list, x0):
    
    if isinstance(x0, (int, float)):
        return (np.abs(x - x0)).argmin()
    
    elif isinstance(x0, (list, tuple, np.ndarray)):
        return [(np.abs(x - xi)).argmin() for xi in x0]