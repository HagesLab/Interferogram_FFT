# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 23:08:58 2020

@author: Chuck
"""
from interferogram_functions import FFT_map, prep_map, Fit_1exp, where_closest
from interferogram_io import fetch_metadata, import_MAP
from interferogram_vis import plot_PL_spectrum, plot_TRPL_decay, plot_TRES
from make_norm_spec import interp, load_spectrum
import matplotlib.pyplot as plt
import numpy as np
import time
import h5py
import os

from scipy.integrate import simpson
import ast

path = r"E:\GEMENI DAQ\NIREOS Complete Example V12_MCS_TimeHarp_32bit Folder\Measurement\20220519\154459"
params_from_INTR_metadata = True        #Import metadata from "...Averaged_MAP..." script - if not using this there may be bugs.
save_data = True                        #Save all plots and TRES data
ImportTRES = False                     #Use this to prevent recalcualting the FFT - must have "..TRES.h5" already savded

# =============================================================================
# If not importing TRES data
# =============================================================================
#trim time-scale to have a smaller data set
max_time_cutoff = 100  #ns
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
timeRange = [-10,80]
PLRange = [560.,800.]
transfer_func = True

#TRES
min_value = 0.001
Gauss_Filter = True
sigmaval = 2   #For Gauss Filter

#PL Plot
AveragePL = False
rangevalPL = [[0,2],[10,20]]  #ns
NormPL = 0

#TRPL Plot
Usemapdata= True      #To maximize TRPL decay - does the same as the Averaged_MAP script

norm_fname = "cuvet_norm_new.txt"

AverageTRPL = False                   #Only if not using mapdata
rangevalTRPL = [[520,540], [590,620]]  #nm    #Only if not using mapdata
NormTRPL = False
BKGTRPL = True
TRPLmin_OM = 1e-4
overrideTRPLrange = False    #For standalone TRPL plot
overidexrange = [-10,1000]    #Only if overriding range

#Fitting TRPL
FitTRPL = False
#fit_range = [[i-1,i+1] for i in range(2, 30, 2)]   #List length must match that of rangeValTRPL / if mapdata then length 1
fit_range = [[20,30]]
fit_on_TRES = True

#Composite TRES
composite_legend = True


# =============================================================================
# Program
# =============================================================================
startTime = time.time()

exper_ID = os.path.split(path)[-1]

pos_data, time_data, map_data = import_MAP(path)

# Auto read params from INTR
if params_from_INTR_metadata:
    metadata_path = os.path.join(path, '{}_FFTmetadata.txt'.format(exper_ID))
    params = fetch_metadata(metadata_path)

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
    
# Validation #
if BKGrange[0] >= BKGrange[1]:
    raise ValueError("Invalid BKGrange {} to {}".format(BKGrange[0], BKGrange[1]))
    
if PLRange[0] >= PLRange[1]:
    raise ValueError("Invalid PLrange {} to {}".format(PLRange[0], PLRange[1]))
    
if timeRange[0] >= timeRange[1]:
    raise ValueError("Invalid PLrange {} to {}".format(timeRange[0], timeRange[1]))
    
for r in rangevalPL:
    if r[0] >= r[1]:
        raise ValueError("Invalid rangevalPL {} to {}".format(r[0], r[1]))

for r in rangevalTRPL:
    if r[0] >= r[1]:
        raise ValueError("Invalid rangevalTRPL {} to {}".format(r[0], r[1]))
    
if np.any(np.diff(time_data) <= 0):
    raise ValueError("time_data is not monotonically ascending")
##############

if save_data:
    allparam = {'max_time_cutoff':max_time_cutoff, 'intfPlot':intfPlot,'intrfxlims':intrfxlims, 'timeRange':timeRange, 'PLRange':PLRange, 'Gauss_Filter':Gauss_Filter, 'sigmaval':sigmaval, 'AveragePL':AveragePL, 'rangevalPL':rangevalPL,'NormPL':NormPL,'Usemapdata':Usemapdata,'AverageTRPL':AverageTRPL,'rangevalTRPL':rangevalTRPL,'NormTRPL':NormTRPL,'BKGTRPL':BKGTRPL,'TRPLmin_OM':TRPLmin_OM,'overrideTRPLrange':overrideTRPLrange , 'overidexrange':overidexrange, 'composite_legend':composite_legend,"apod_width":apodization_width, "apod_type":apod_type, "do_resample":resample, "resample_factor":resample_factor, "do_shift":shift, "do_padzeros":pad_test, "pad_factor":padfactor, "do_mean_sub":mean_sub, "shift_factor":shift_factor,"background_subtract":BKGsub,"background_range_low":BKGrange[0], "background_range_high":BKGrange[1],"baseline_sub_state":baseline_sub_state}
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
        index = where_closest(time_data, BKGrange)
        BKGval = np.mean(map_data[:,index[0]:index[1]],axis=1)
        map_data = map_data - np.reshape(BKGval, (len(BKGval), 1))

    #if baseline_sub_state:


    #trim time-scale
    index = where_closest(time_data, max_time_cutoff)
    map_data = map_data[:,0:np.max(index)]
    time_data = time_data[0:np.max(index)]

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
        

# Averged PL over given range
if AveragePL:
    averaged_PL_spec = np.empty([len(rangevalPL),len(wave)])
    for i in range(len(rangevalPL)):
        index = where_closest(time_data, rangevalPL[i])
        newarr= np.sum(build_TRES[index[0]:index[1],:],axis=0)
        if NormPL:
            newarr = (newarr-newarr.min())/(newarr.max()-newarr.min())
        averaged_PL_spec[i] = newarr
else:
    averaged_PL_spec = np.sum(build_TRES,axis=0)
    if NormPL:
        averaged_PL_spec = (averaged_PL_spec-averaged_PL_spec.min())/(averaged_PL_spec.max()-averaged_PL_spec.min())
        
        
# integral TRPL decay over given range
if Usemapdata:
    AVGTRPL = np.fliplr(map_data.T)
    AverageTRPL=False

else:
    AVGTRPL = build_TRES

if AverageTRPL:
    pass
else:
    rangevalTRPL = [[np.min(wave), np.max(wave)]]
    
integralTRPL = np.empty([len(rangevalTRPL),len(time_data)])
for i in range(len(rangevalTRPL)):
    # Closest?
    index = where_closest(wave, rangevalTRPL[i])
    if Usemapdata:
        integralTRPL[i] = np.sum(AVGTRPL[:,index[0]:index[1]], axis=1)
    else:
        integralTRPL[i] = simpson(AVGTRPL[:,index[0]:index[1]], x=wave[index[0]:index[1]], axis=1)
        
    if NormTRPL:
        integralTRPL[i] = integralTRPL[i]/integralTRPL[i].max()
    if BKGTRPL:
        index = where_closest(time_data, BKGrange)
        BKGval = np.mean(integralTRPL[i][index[0]:index[1]])
        integralTRPL[i] -= BKGval
        
# Fit a monoexponential to integral TRPL
if FitTRPL:
    TRPL_fit_list, time_fit_list, fit_label_list = [],[],[]
    
    for i in range(len(integralTRPL)):
        TRPL_fit, time_fit, fit_label = Fit_1exp(integralTRPL[i],time_data,fit_range[i])
        TRPL_fit = list(TRPL_fit)
        time_fit = list(time_fit)
        fit_label = list(fit_label)
        TRPL_fit_list.append(TRPL_fit)
        time_fit_list.append(time_fit)
        fit_label_list.append(fit_label)
        
    fit = (time_fit_list, TRPL_fit_list, fit_label_list)
        
# Plot the results (TRES)
indexWL = where_closest(wave, PLRange)
WLPlot=wave[indexWL[0]:indexWL[1]]
indext = where_closest(time_data, timeRange)
tplot=time_data[indext[0]:indext[1]]

timemesh, wavemesh = np.meshgrid(WLPlot,tplot)
TRES=build_TRES[indext[0]:indext[1],indexWL[0]:indexWL[1]]
TRES = np.where(TRES<min_value,min_value,TRES)

plot_TRES(timemesh, wavemesh, TRES, Gauss_Filter=Gauss_Filter, sigma=sigmaval)

# Plot averaged PL
PLname = os.path.join(path, '{}_PLPlot.png'.format(exper_ID))
labels = ["{} to {} ns".format(PL_range[0], PL_range[1]) for PL_range in rangevalPL]
plot_PL_spectrum(wave, averaged_PL_spec, labels, PLRange[0], PLRange[1], export=PLname)

# Plot averaged TRPL
if AverageTRPL:
    labels = ["{} to {} nm".format(PL_range[0], PL_range[1]) for PL_range in rangevalTRPL]
else:
    labels = [None]
    
if overrideTRPLrange:
    start_time, end_time = overidexrange[0], overidexrange[1]
else:
    start_time, end_time = timeRange[0], timeRange[1]
    
PLname = os.path.join(path, '{}_TRPLPlot.png'.format(exper_ID))
plot_TRPL_decay(time_data, integralTRPL, TRPLmin_OM, labels=labels, 
                start_time=start_time, end_time=end_time, export=PLname)

if FitTRPL:
    TRPLFitname = os.path.join(path, '{}_TRPLFitPlot.png'.format(exper_ID))
    plot_TRPL_decay(time_data, integralTRPL, TRPLmin_OM, labels=labels, 
                    start_time=start_time, end_time=end_time, fit=fit, export=TRPLFitname)
    
#Composite TRES
fig = plt.figure(6,dpi=120)
grid = plt.GridSpec(2, 3, height_ratios=[1, 1/3],width_ratios=[0.8,2,0.44],wspace=0.05,hspace=0.05)

main_ax = fig.add_subplot(grid[:-1, 1:])
main_ax.set_xlim(PLRange)
main_ax.set_ylim(timeRange)
main_ax.label_outer()
plot_TRES(timemesh, wavemesh, TRES, ax=main_ax, fig=fig, Gauss_Filter=Gauss_Filter, sigma=sigmaval)

TRPL = fig.add_subplot(grid[:-1, 0], xticklabels=[],sharey=main_ax)
TRPL.invert_xaxis()
TRPL.set_xscale('log')
TRPL.set_xlim(2*np.max(integralTRPL),np.max(integralTRPL)*TRPLmin_OM)

for i in range(len(integralTRPL)):
    TRPL.plot(integralTRPL[i],time_data,label=labels[i])
    if fit_on_TRES and FitTRPL:
        TRPL.plot(np.array(TRPL_fit_list[i]),np.array(time_fit_list[i]),'k--',label = ''.join(fit_label_list[i]))

TRPL.set(ylabel='Time / ns')

PL = fig.add_subplot(grid[-1, 1], yticklabels=[], sharex=main_ax)
if AveragePL:
    for i in range(len(averaged_PL_spec)):
        PL.plot(wave,averaged_PL_spec[i],label="{} to {} ns".format(*rangevalPL[i]))
else:
    PL.plot(wave,averaged_PL_spec)
PL.set(xlabel='Wavelength / nm')
if composite_legend:
    fig.legend(loc='lower left', bbox_to_anchor=(0, 0.1))
if save_data:
    TRESname = os.path.join(path, '{}_TRESPlot.png'.format(exper_ID))
    plt.savefig(TRESname)

#Write 2D Data Set
if save_data:
    outputfilename = os.path.join(path, '{}_TRES.h5'.format(exper_ID))
    hf = h5py.File(outputfilename,'w')
    hf.create_dataset('TRES Data', data=build_TRES)
    hf.create_dataset('Time Data', data=time_data)
    hf.create_dataset('Wavelength', data=wave)
    hf.create_dataset('Metadata', data=[[min(PLRange),max(PLRange)]])
    hf.create_dataset('PL', data=averaged_PL_spec)
    if AveragePL:
        hf.create_dataset('PL Metadata', data=rangevalPL)
    hf.create_dataset('TRPL', data=integralTRPL)
    if AverageTRPL:
        hf.create_dataset('TRPL Metadata', data=rangevalTRPL)
    hf.close()

