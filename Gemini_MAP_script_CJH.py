# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 23:08:58 2020

@author: Chuck
"""
from interferogram_functions import FFT_map, import_MAP, prep_map, fetch_metadata, Fit_1exp
from make_norm_spec import interp, load_spectrum
import matplotlib.pyplot as plt
import numpy as np
import time
import h5py
import os
import matplotlib.colors
from matplotlib.ticker import LogLocator
from scipy import ndimage
import ast

path = r"E:\GEMENI DAQ\NIREOS Complete Example V12_MCS_TimeHarp_32bit Folder\Measurement\20220509\170820"
params_from_INTR_metadata = True        #Import metadata from "...Averaged_MAP..." script - if not using this there may be bugs.
save_data = True                        #Save all plots and TRES data
ImportTRES = False                       #Use this to prevent recalcualting the FFT - must have "..TRES.h5" already savded

# =============================================================================
# If not importing TRES data
# =============================================================================
#trim time-scale to have a smaller data set
rangeval = 100  #ns
#Plot Interferogram? (Background Subtracted and Shifted)
intfPlot = False
intrfxlims = "Full"   #if == "Full" no restriction, full data. Otherwise define range
#intrfxlims = [-0.25,0.25]

# =============================================================================
# If not importing params from metadata:
# =============================================================================

apodization_width=0.5
apod_type="BH"    # "None" "Gauss" "Triangle" "Boxcar" "BH"
resample=True
resample_factor=4
shift=False
pad_test=True
padfactor=16
mean_sub = True
baseline_sub_state = False
BKGsub = True
bkg_limit = -3  #ns  Before t_max
shift_factor = 0.005837467299965371       #Hard to compute shift on time-resolved data at each time, do it manually


# =============================================================================
# Plotting Metadata
# =============================================================================

#All Plots
timeRange = [-10,60]
PLRange = [550.,800.]

#TRES
min_value = 50
Gauss_Filter = True
sigmaval = 2   #For Gauss Filter

#PL Plot
AveragePL = False
rangevalPL = [[0,1],[4,20]]  #ns
NormPL = False

#TRPL Plot
Usemapdata= False      #To maximize TRPL decay
transfer_func = True # Only if not using mapdata
norm_fname = "cuvet_norm_0.txt"
AverageTRPL = False                     #Only if not using mapdata
rangevalTRPL = [[650,670]]  #nm    #Only if not using mapdata
NormTRPL = False
BKGTRPL = True
TRPLmin_OM = 1e-4
overrideTRPLrange = False    #For standalone TRPL plot
overidexrange = [-10,1000]    #Only if overriding range

#Fitting TRPL
FitTRPL = True
fit_range = [[5,12]]   #List length must match the number of TRPL curves (line 64) / if mapdata then length 1
fit_on_TRES = True

#Composite TRES
composite_legend = True


# =============================================================================
# Program
# =============================================================================
startTime = time.time()
outputfilename = path + "\\" + os.path.split(path)[-1] + '_TRES.h5'
pos_data, time_data, map_data = import_MAP(path)

# Auto read params from INTR
if params_from_INTR_metadata:
    params = fetch_metadata(path,"{}\\" + os.path.split(path)[-1] + '_FFTmetadata.txt')

if params_from_INTR_metadata:
    apodization_width = params['apod_width']
    apod_type = params['apod_type']
    resample = ast.literal_eval(params['do_resample'])
    resample_factor = params['resample_factor']
    shift = ast.literal_eval(params["do_shift"])
    pad_test = ast.literal_eval(params['do_padzeros'])
    padfactor = params['pad_factor']
    mean_sub = ast.literal_eval(params['do_mean_sub'])
    baseline_sub_state = ast.literal_eval(params['do_baseline_sub'])
    BKGsub = ast.literal_eval(params['background_subtract'])
    BKGrange = [params['background_range_low'],params['background_range_high']]
    shift_factor = params['shift_factor']

else:
    BKGrange = np.array([time_data[0],bkg_limit],dtype='float')  #ns

if save_data:
    allparam = {'rangeval':rangeval, 'intfPlot':intfPlot,'intrfxlims':intrfxlims, 'timeRange':timeRange, 'PLRange':PLRange, 'Gauss_Filter':Gauss_Filter, 'sigmaval':sigmaval, 'AveragePL':AveragePL, 'rangevalPL':rangevalPL,'NormPL':NormPL,'Usemapdata':Usemapdata,'AverageTRPL':AverageTRPL,'rangevalTRPL':rangevalTRPL,'NormTRPL':NormTRPL,'BKGTRPL':BKGTRPL,'TRPLmin_OM':TRPLmin_OM,'overrideTRPLrange':overrideTRPLrange , 'overidexrange':overidexrange, 'composite_legend':composite_legend,"apod_width":apodization_width, "apod_type":apod_type, "do_resample":resample, "resample_factor":resample_factor, "do_shift":shift, "do_padzeros":pad_test, "pad_factor":padfactor, "do_mean_sub":mean_sub, "shift_factor":shift_factor,"background_subtract":BKGsub,"background_range_low":BKGrange[0], "background_range_high":BKGrange[1],"baseline_sub_state":baseline_sub_state}
    with open((r"{}\\"+os.path.split(path)[-1]+"_MapParams.txt").format(path), 'w+') as ofstream:
        ofstream.write("# Params used to make plots in Gemini_MAP_script_CJH.py")
        for param, val in allparam.items():
            ofstream.write("\n{}:\t{}".format(param, val))


if ImportTRES:
    importfilename = path + "\\" + os.path.split(path)[-1] + '_TRES.h5'
    hf = h5py.File(importfilename, 'r')
    build_TRES = np.array(hf.get('TRES Data'))
    time_data=np.array(hf.get('Time Data'))
    t_max = time_data[np.array(np.where(np.mean(map_data,axis=0)==np.max(np.mean(map_data,axis=0)))[0],dtype="int")]
    index = [(np.abs(time_data-np.min(time_data))).argmin(),(np.abs(time_data-np.max(time_data))).argmin()+1]
    map_data = map_data[:,np.min(index):np.max(index)]
    wave=np.array(hf.get('Wavelength'))
    hf.close()


else:
    t_max = time_data[np.array(np.where(np.mean(map_data,axis=0)==np.max(np.mean(map_data,axis=0)))[0],dtype="int")]
    time_data = time_data - t_max

   #Manual Shift of Peak
    if shift:
        pos_data = pos_data - shift_factor    #Taken from shift_factor output in the _INTR analysis script for this data

   #Background Subtract TRPL Curves
    if BKGsub:
        index = [(np.abs(time_data-np.min(BKGrange))).argmin(),(np.abs(time_data-np.max(BKGrange))).argmin()]
        BKGval = np.mean(map_data[:,np.min(index):np.max(index)],axis=1)
        map_data = map_data - np.reshape(BKGval, (len(BKGval), 1))

   #if baseline_sub_state:


    #trim time-scale
    index = (np.abs(time_data-np.max(rangeval))).argmin()
    map_data = map_data[:,0:np.max(index)]
    time_data=time_data[0:np.max(index)]

    #Plot Raw Data (Background Subtracted and Shifted)
    if intfPlot:
        raw_timemesh, raw_posmesh = np.meshgrid(time_data,pos_data)
        plt.figure(1, dpi=120)
        plt.contourf(raw_posmesh,raw_timemesh,np.log(map_data))
        plt.ylabel('Time / ns')
        plt.xlabel('Position / mm')
        if intrfxlims != "Full":
            plt.xlim(intrfxlims.min(),intrfxlims.max())

    preFFT_pos, preFFT_map = prep_map(pos_data,map_data,apodization_width,apod_type,resample,resample_factor,shift,pad_test,padfactor,mean_sub)
    print("Took {} sec".format(time.time() - startTime))

    #Perform FFT
    wave, build_TRES = FFT_map(preFFT_pos, preFFT_map)
    wave=wave[::-1]
    build_TRES=np.fliplr(np.array(build_TRES,dtype="float").T)
    print("Took {} sec".format(time.time() - startTime))

    if transfer_func:
        norm_waves, norm_data = load_spectrum(norm_fname)
        norm_waves, norm_data = interp(norm_waves, norm_data, PLRange[0], PLRange[1], 1)
        wave, build_TRES = interp(wave, build_TRES.T, PLRange[0], PLRange[1], 1)
        build_TRES = (build_TRES.T / norm_data).T
        build_TRES = build_TRES.T
        


# Plot the results (TRES)
indexWL = [(np.abs(wave-np.min(PLRange))).argmin(),(np.abs(wave-np.max(PLRange))).argmin()]
WLPlot=wave[indexWL[0]:indexWL[1]]
indext = [(np.abs(time_data-np.min(timeRange))).argmin(),(np.abs(time_data-np.max(timeRange))).argmin()]
tplot=time_data[indext[0]:indext[1]]

timemesh, wavemesh = np.meshgrid(WLPlot,tplot)
TRESplot=build_TRES[indext[0]:indext[1],indexWL[0]:indexWL[1]]
min_value = np.amin(TRESplot)
TRESplot = np.where(TRESplot<min_value,min_value,TRESplot)

fig = plt.figure(2,dpi=120)
ax = fig.add_subplot()
if Gauss_Filter:
    TRESplot = ndimage.gaussian_filter(TRESplot, sigma=sigmaval)
norm= matplotlib.colors.LogNorm(vmin=min_value, vmax=TRESplot.max())
levels = np.logspace(np.log10(np.min(min_value)),np.log10(np.max(TRESplot)),num=50)
cs = ax.contourf(timemesh,wavemesh,TRESplot,levels=levels,norm=norm, cmap='plasma')
cbar = fig.colorbar(cs)
cbar.ax.yaxis.set_major_locator(LogLocator())
cbar.set_ticks(cbar.ax.yaxis.get_major_locator().tick_values(TRESplot.min(), TRESplot.max()))
ax.set_ylabel('Time / ns')
ax.set_xlabel('Wavelength / nm')

#Plot averged PL over given range
if AveragePL:
    plot_TRES = np.empty([len(rangevalPL),len(wave)])
    for i in range(len(rangevalPL)):
        index = [(np.abs(time_data-np.min(rangevalPL[i]))).argmin(),(np.abs(time_data-np.max(rangevalPL[i]))).argmin()]
        newarr= np.sum(build_TRES[np.min(index):np.max(index),:],axis=0)
        if NormPL:
            newarr = (newarr-newarr.min())/(newarr.max()-newarr.min())
        plot_TRES[i] = newarr
else:
    plot_TRES = np.sum(build_TRES,axis=0)
    if NormPL:
        plot_TRES = (plot_TRES-plot_TRES.min())/(plot_TRES.max()-plot_TRES.min())

plt.figure(3, dpi=120)
plt.xlim(min(PLRange),max(PLRange))
plt.title("Average PL")
plt.ylabel('Counts / a.u.')
plt.xlabel('Wavelength / nm')
plt.yscale('linear')
if AveragePL:
    for i in range(len(plot_TRES)):
        plt.plot(wave,plot_TRES[i],label=str(min(rangevalPL[i])) + " to " + str(max(rangevalPL[i])) + " ns")
    plt.legend()
else:
    plt.plot(wave,plot_TRES)
if save_data:
    PLname = path + "\\" + os.path.split(path)[-1] + '_PLPlot.png'
    plt.savefig(PLname)


#Plot integral TRPL decay over given range
if Usemapdata:
    AVGTRPL = np.fliplr(map_data.T)
    AverageTRPL=False
else:
    AVGTRPL = build_TRES
    


if AverageTRPL:
    integralTRPL = np.empty([len(rangevalTRPL),len(time_data)])
    for i in range(len(rangevalTRPL)):
        index = [(np.abs(wave-np.min(rangevalTRPL[i]))).argmin(),(np.abs(wave-np.max(rangevalTRPL[i]))).argmin()]
        TRPLarray = np.sum(AVGTRPL[:,np.min(index):np.max(index)],axis=1)
        if NormTRPL:
            TRPLarray = TRPLarray/TRPLarray.max()
        if BKGTRPL:
            index = [(np.abs(time_data-np.min(BKGrange))).argmin(),(np.abs(time_data-np.max(BKGrange))).argmin()]
            BKGval = np.mean(TRPLarray[np.min(index):np.max(index)])
            TRPLarray = TRPLarray - BKGval
        integralTRPL[i] = TRPLarray
else:
    integralTRPL = np.sum(AVGTRPL,axis=1)
    if NormTRPL:
        integralTRPL = integralTRPL/integralTRPL.max()
    if BKGTRPL:
        index = [(np.abs(time_data-np.min(BKGrange))).argmin(),(np.abs(time_data-np.max(BKGrange))).argmin()]
        BKGval = np.mean(integralTRPL[np.min(index):np.max(index)])
        integralTRPL = integralTRPL - BKGval




plt.figure(4, dpi=120)
plt.title("Integral TRPL")
plt.xlabel('Time / ns')
plt.ylabel('Counts / a.u.')
plt.yscale('log')
plt.ylim(np.max(integralTRPL)*TRPLmin_OM,2*np.max(integralTRPL))
plt.xlim(min(timeRange),max(timeRange))
if overrideTRPLrange:
    plt.xlim(min(overidexrange),max(overidexrange))
if AverageTRPL:
    for i in range(len(integralTRPL)):
        plt.plot(time_data,integralTRPL[i],label=str(min(np.array(rangevalTRPL[i],dtype='int16'))) + " to " + str(max(np.array(rangevalTRPL[i],dtype='int16'))) + " nm")
        plt.legend()
else:
    plt.plot(time_data,integralTRPL)
if save_data:
    TRPLname = path + "\\" + os.path.split(path)[-1] + '_TRPLPlot.png'
    plt.savefig(TRPLname)

if FitTRPL:
    if AverageTRPL:
        TRPL_fit_list, time_fit_list, fit_label_list = [],[],[]
        for i in range(len(integralTRPL)):
            TRPL_fit, time_fit, fit_label, popt, perr = Fit_1exp(integralTRPL[i],time_data,fit_range[i])
            TRPL_fit = list(TRPL_fit)
            time_fit = list(time_fit)
            fit_label = list(fit_label)
            TRPL_fit_list.append(TRPL_fit)
            time_fit_list.append(time_fit)
            fit_label_list.append(fit_label)
    else:
        TRPL_fit, time_fit, fit_label, popt, perr = Fit_1exp(integralTRPL,time_data,fit_range[0])

if FitTRPL:
    plt.figure(5, dpi=120)
    plt.title("Integral TRPL")
    plt.xlabel('Time / ns')
    plt.ylabel('Counts / a.u.')
    plt.yscale('log')
    plt.ylim(np.max(integralTRPL)*TRPLmin_OM,2*np.max(integralTRPL))
    plt.xlim(min(timeRange),max(timeRange))
    if overrideTRPLrange:
        plt.xlim(min(overidexrange),max(overidexrange))
    if AverageTRPL:
        for i in range(len(integralTRPL)):
            plt.plot(time_data,integralTRPL[i],label=str(min(np.array(rangevalTRPL[i],dtype='int16'))) + " to " + str(max(np.array(rangevalTRPL[i],dtype='int16'))) + " nm")
            plt.plot(np.array(time_fit_list[i]),np.array(TRPL_fit_list[i]),'k--',label = ''.join(fit_label_list[i]))
    else:
        plt.plot(time_data,integralTRPL)
        plt.plot(time_fit,TRPL_fit,'k--',label = fit_label)
    plt.legend()
    if save_data:
        TRPLFitname = path + "\\" + os.path.split(path)[-1] + '_TRPLFitPlot.png'
        plt.savefig(TRPLFitname)

#Composite TRES
fig = plt.figure(6,dpi=120)
grid = plt.GridSpec(2, 3, height_ratios=[1, 1/3],width_ratios=[0.8,2,0.44],wspace=0.05,hspace=0.05)

main_ax = fig.add_subplot(grid[:-1, 1:])
main_ax.set_xlim(PLRange)
main_ax.set_ylim(timeRange)
main_ax.label_outer()

cs = main_ax.contourf(timemesh,wavemesh,TRESplot,levels=levels,norm=norm, cmap='plasma')
cbar = fig.colorbar(cs)
cbar.ax.yaxis.set_major_locator(LogLocator())
cbar.set_ticks(cbar.ax.yaxis.get_major_locator().tick_values(TRESplot.min(), TRESplot.max()))
#cbar.ax.yaxis.set_major_locator(LogLocator())
#cbar.set_ticks(cbar.ax.yaxis.get_major_locator().tick_values(TRESplot.min(), TRESplot.max()))
#main_ax.tick_params(bottom='off')

TRPL = fig.add_subplot(grid[:-1, 0], xticklabels=[],sharey=main_ax)
TRPL.invert_xaxis()
TRPL.set_xscale('log')
TRPL.set_xlim(2*np.max(integralTRPL),np.max(integralTRPL)*TRPLmin_OM)
if fit_on_TRES and FitTRPL:
    if AverageTRPL:
        for i in range(len(integralTRPL)):
            TRPL.plot(integralTRPL[i],time_data,label=str(min(np.array(rangevalTRPL[i],dtype='int16'))) + " to " + str(max(np.array(rangevalTRPL[i],dtype='int16'))) + " nm")
            TRPL.plot(np.array(TRPL_fit_list[i]),np.array(time_fit_list[i]),'k--',label = ''.join(fit_label_list[i]))
    else:
        TRPL.plot(integralTRPL,time_data)
        TRPL.plot(TRPL_fit,time_fit,'k--',label = fit_label)
else:
    if AverageTRPL:
        for i in range(len(integralTRPL)):
            TRPL.plot(integralTRPL[i],time_data,label=str(min(rangevalTRPL[i])) + " to " + str(max(rangevalTRPL[i])) + " nm")
    else:
        TRPL.plot(integralTRPL,time_data)
TRPL.set(ylabel='Time / ns')

PL = fig.add_subplot(grid[-1, 1], yticklabels=[], sharex=main_ax)
if AveragePL:
    for i in range(len(plot_TRES)):
        PL.plot(wave,plot_TRES[i],label=str(min(rangevalPL[i])) + " to " + str(max(rangevalPL[i])) + " ns")
else:
    PL.plot(wave,plot_TRES)
PL.set(xlabel='Wavelength / nm')
if composite_legend:
    fig.legend(loc='lower left', bbox_to_anchor=(0, 0.1))
if save_data:
    TRESname = path + "\\" + os.path.split(path)[-1] + '_TRESPlot.png'
    plt.savefig(TRESname)

#Write 2D Data Set
if save_data:
    hf = h5py.File(outputfilename,'w')
    hf.create_dataset('TRES Data', data=build_TRES)
    hf.create_dataset('Time Data', data=time_data)
    hf.create_dataset('Wavelength', data=wave)
    hf.create_dataset('Metadata', data=[[min(PLRange),max(PLRange)]])
    hf.create_dataset('PL', data=plot_TRES)
    if AveragePL:
        hf.create_dataset('PL Metadata', data=rangevalPL)
    hf.create_dataset('TRPL', data=integralTRPL)
    if AverageTRPL:
        hf.create_dataset('TRPL Metadata', data=rangevalTRPL)
    hf.close()

