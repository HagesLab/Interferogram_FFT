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
import glob
import csv
import os
import pandas as pd

def import_INTR(path):
    allFiles = []
    b = glob.glob(path + "/*.txt")  # glob. lists the filename and path
    allFiles.extend(b)
    
    for i in range(len(allFiles)):
        namestest = allFiles[i].split('_')[-1]
        if namestest == r"POS.txt":
            with open(allFiles[i],'r') as j:
                pos_data = np.array(list(csv.reader(j,delimiter="\t")),dtype='float')[0,:]
                samplename = os.path.split(allFiles[i])[-1].split(".")[0]
                samplename = samplename.replace("_POS","")
        elif namestest == r"INTR.txt":
            with open(allFiles[i],'r') as j:
                intr_data = np.array(list(csv.reader(j,delimiter="\t")),dtype='float')[0,:]
    return pos_data, intr_data

    

def prep_interferogram(pos_data,intr_data,apodization_width,apod_type="BH",mean_sub="True",resample="True",resample_factor=2,shift="True",pad_test="True",padfactor=4,plots="True"):
    intr_data, pos_data = Interf_calibration1D.interf_calibration1D(intr_data,pos_data)
    shiftfactor = 0
    if mean_sub == "True":
        intr_data = intr_data-np.mean(intr_data)
    if resample == "True":
        fcubic = interp1d(pos_data, intr_data, kind='cubic')
        pos_data = np.linspace(pos_data[0],pos_data[-1],endpoint=True,num=len(pos_data)*resample_factor)
        intr_data = fcubic(pos_data)
    index_pos = np.argmin(abs(pos_data))
    if shift == "True":
        index_intr = np.argmax(intr_data)  
        shiftfactor=(pos_data[index_intr]-pos_data[index_pos])
        pos_data = pos_data-shiftfactor
    if apod_type != "None":
        raw_intr = intr_data
        raw_pos = pos_data
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
    if pad_test == "True":
        padlen = int(padfactor*(2**np.ceil(np.log(len(intr_data))/np.log(2))))-len(intr_data)
        intr_data = np.pad(intr_data,(0,int(padlen)),'constant', constant_values=(0))
        pos_data = np.append(pos_data,np.linspace(pos_data[-1]+posdiff,(pos_data[-1]+posdiff)+padlen*posdiff,num=padlen))
        
    left_axis_intr, right_axis_intr = intr_data[0:index_pos+1], intr_data[index_pos+1:]
    left_axis_pos, right_axis_pos = pos_data[0:index_pos+1], pos_data[index_pos+1:]
    preFFT_data, preFFT_pos = np.concatenate((right_axis_intr,left_axis_intr)), np.concatenate((right_axis_pos,left_axis_pos))   
    
    if plots == "True":
        #Plot interferogram data    
        plt.figure(1, dpi=120)
        if apod_type != "None":
            plt.plot(raw_pos,raw_intr/np.max(raw_intr))
            plt.plot(pos_data[0:len(raw_pos)],(intr_data/np.max(intr_data))[0:len(raw_pos)])
            plt.plot(raw_pos,intr_func)
        else:
            plt.plot(pos_data,intr_data)
        plt.xlabel('Position / mm')
        plt.ylabel('Counts / a.u.')
        plt.show()
        
        plt.figure(2, dpi=120)
        plt.plot(preFFT_data)
        plt.show()
        
    return preFFT_pos, preFFT_data, shiftfactor

def FFT_intr(preFFT_pos,preFFT_data,plots="False",correct="True",scale="linear"):
    #Perform FFT
    FFT_intr = np.fft.rfft(preFFT_data)
    if correct == "True":
        #Phase Corrected
        FFT_real = FFT_intr.real*np.cos(np.angle(FFT_intr))
        FFT_imag = FFT_intr.imag*np.sin(np.angle(FFT_intr))
        FFT_final = FFT_real + FFT_imag
    elif correct == "False":
        FFT_final = FFT_intr.real + FFT_intr.imag
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
    #Trim FFT data to calibrated WLs
    FFT_intr_trim=[]
    freq_trim=[]
    for i in range(len(FFT_intr)):
        if freq[i] >= np.min(reciprocal) and freq[i] <= np.max(reciprocal):
            FFT_intr_trim = np.append(FFT_intr_trim,FFT_final[i])
            freq_trim = np.append(freq_trim,freq[i])
    #Compute WL
    fn = interp1d(reciprocal, wavelength, kind="linear")
    wave = fn(freq_trim)
    
    if plots == "True":
        plt.figure(1, dpi=120)
        plt.plot(freq,FFT_intr.real,label="Real")
        plt.plot(freq,FFT_intr.imag,label="Imag")
        plt.yscale(scale)
        plt.title("Raw FFT")
        plt.xlabel("Freq.")
        plt.legend()
        plt.show()

        if correct == "True":
            plt.figure(2, dpi=120)
            plt.plot(freq,FFT_real,label="Real")
            plt.plot(freq,FFT_imag,label="Imag")
            #plt.plot(freq,np.angle(FFT_intr))
            plt.yscale(scale)
            plt.title("Corrected FFT")
            plt.xlabel("Freq.")
            plt.legend()
            plt.show()
    
    return wave, FFT_intr_trim

def import_MAP(path):
    allFiles = []
    b = glob.glob(path + "/*.txt")  # glob. lists the filename and path
    allFiles.extend(b)
    
    
    for i in range(len(allFiles)):
        namestest = allFiles[i].split('_')[-1]
        if namestest == r"POS.txt":
            with open(allFiles[i],'r') as j:
                pos_data = np.array(list(csv.reader(j,delimiter="\t")),dtype='float')[0,:]
                samplename = os.path.split(allFiles[i])[-1].split(".")[0]
                samplename = samplename.replace("_POS","")
        elif namestest == r"MAP.txt":
            with open(allFiles[i],'r') as j:
                map_data = np.array(list(csv.reader(j,delimiter="\t")),dtype='float')
        elif namestest == r"TIME.txt":
            with open(allFiles[i],'r') as j:
                time_data = np.array(list(csv.reader(j,delimiter="\t")),dtype='float')[0,:]
    
    # Discard the zero columns
    standard=np.std(map_data, axis=0)
    indexes=np.asarray(np.nonzero(standard))
    map_data=map_data[:,indexes[0][0]:indexes[0][-1]]
    time_data =time_data[indexes[0][0]:indexes[0][-1]]/1e3
    
    return pos_data, time_data, map_data
    

def prep_map(pos_data,map_data,apodization_width,apod_type="BH",mean_sub="True",resample="True",resample_factor=2,shift="True",pad_test="True",padfactor=4,plots="True"):
    map_data, pos_data = Interf_calibration2D.interf_calibration2D(map_data,pos_data)   
    prep_build=[]
    pos_data_raw = pos_data
    for i in range(map_data.shape[1]):
        intr_data = map_data[:,i]
        pos_data = pos_data_raw
        if mean_sub == "True":
            intr_data = intr_data-np.mean(intr_data)
        if resample == "True":
            fcubic = interp1d(pos_data, intr_data, kind='cubic')
            pos_data = np.linspace(pos_data[0],pos_data[-1],endpoint=True,num=len(pos_data)*resample_factor)
            intr_data = fcubic(pos_data)
        index_pos = np.argmin(abs(pos_data))
        if shift == "True":
            index_intr = np.argmax(intr_data)  
            shiftfactor=(pos_data[index_intr]-pos_data[index_pos])
            pos_data = pos_data-shiftfactor
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
        if pad_test == "True":
            padlen = int(padfactor*(2**np.ceil(np.log(len(intr_data))/np.log(2))))-len(intr_data)
            intr_data = np.pad(intr_data,(0,int(padlen)),'constant', constant_values=(0))
            pos_data = np.append(pos_data,np.linspace(pos_data[-1]+posdiff,(pos_data[-1]+posdiff)+padlen*posdiff,num=padlen))
        
        index_pos = np.argmin(abs(pos_data))    
        left_axis_intr, right_axis_intr = intr_data[0:index_pos+1], intr_data[index_pos+1:]
        left_axis_pos, right_axis_pos = pos_data[0:index_pos+1], pos_data[index_pos+1:]
        preFFT_data, preFFT_pos = np.concatenate((right_axis_intr,left_axis_intr)), np.concatenate((right_axis_pos,left_axis_pos))   
        prep_build.append(preFFT_data)
        
    return preFFT_pos, np.array(prep_build,dtype='float').T