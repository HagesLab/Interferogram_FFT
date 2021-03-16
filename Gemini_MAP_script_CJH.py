# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 23:08:58 2020

@author: Chuck
"""
from interferogram_functions import FFT_intr, import_MAP, prep_map, fetch_metadata
import matplotlib.pyplot as plt
import numpy as np
import time
import h5py
import os
import matplotlib.colors
from matplotlib.ticker import LogLocator
from scipy import ndimage
import ast


path = r"C:\Users\c.hages\Dropbox (UFL)\UF\TRPL Computer\Aaron\144620"
params_from_INTR_metadata = True
save_data = True
ImportTRES = True
TRES_param_override = False


# =============================================================================
# If not importing params from metadata:
# =============================================================================
manual_shift_factor = 0.005837467299965371
save_params = False

# =============================================================================
# If not importing TRES data
# =============================================================================
#trim time-scale
rangeval = 1500  #ns
#Plot Interferogram? (Background Subtracted and Shifted)
intfPlot = True
intrfxlims = "Full"
#intrfxlims = [-0.25,0.25]

# =============================================================================
# Plotting Metadata
# =============================================================================

#All Plots
timeRange = [-0.75,25]
PLRange = [600.,950.]

#TRES
min_value = 32
Gauss_Filter = True
sigmaval = 2   #For Gauss Filter

#PL Plot
AveragePL = False
rangevalPL = [[0,1],[5,30]]  #ns
NormPL = True

#TRPL Plot
Usemapdata=True
AverageTRPL = True                              #Only if not using mapdata
rangevalTRPL = [[600.,800.]]  #nm    #Only if not using mapdata
NormTRPL = False
BKGTRPL = True
TRPLmin_OM = 1e-4
overrideTRPLrange = False    #For standalone TRPL plot
overidexrange = [-10,1000]    #Only if overriding range

#Composite TRES
composite_legend = False


# =============================================================================
# Program
# =============================================================================
startTime = time.time()
outputfilename = path + "\\" + os.path.split(path)[-1] + '_TRES.h5'
pos_data, time_data, map_data = import_MAP(path)

#Auto read all params
if TRES_param_override:
    allparam = fetch_metadata(path,"{}\\Plot_Params.txt")
    manual_shift_factor = allparam['manual_shift_factor']
    save_params = allparam['save_params']
    rangeval = allparam['rangeval']
    intfPlot = ast.literal_eval(allparam['intfPlot'])
    intrfxlims = allparam["intrfxlims"]
    timeRange = ast.literal_eval(allparam['timeRange'])
    PLRange = ast.literal_eval(allparam['PLRange'])
    Gauss_Filter = ast.literal_eval(allparam['Gauss_Filter'])
    sigmaval = allparam['sigmaval']
    AveragePL = ast.literal_eval(allparam['AveragePL'])
    rangevalPL = ast.literal_eval(allparam['rangevalPL'])
    NormPL = ast.literal_eval(allparam['NormPL'])
    Usemapdata=ast.literal_eval(allparam['Usemapdata'])
    AverageTRPL = ast.literal_eval(allparam['AverageTRPL'])
    rangevalTRPL = ast.literal_eval(allparam['rangevalTRPL'])
    NormTRPL = ast.literal_eval(allparam['NormTRPL'])
    BKGTRPL = ast.literal_eval(allparam['BKGTRPL'])
    TRPLmin_OM = allparam['TRPLmin_OM']
    overrideTRPLrange = ast.literal_eval(allparam['overrideTRPLrange'])
    overidexrange = ast.literal_eval(allparam['overidexrange'])
    composite_legend = ast.literal_eval(allparam['composite_legend'])

# Auto read params from INTR
if params_from_INTR_metadata:
    params = fetch_metadata(path,"{}\\Param_Import_metadata.txt")

if params_from_INTR_metadata:
    apodization_width = params['apod_width']
    apod_type = params['apod_type']
    resample = params['do_resample']
    resample_factor = params['resample_factor']
    shift = params["do_shift"]
    pad_test = params['do_padzeros']
    padfactor = params['pad_factor']
    mean_sub = params['do_mean_sub']
    BKGsub = params['background_subtract']
    BKGrange = [params['background_range_low'],params['background_range_high']]
    shift_factor = params['shift_factor']
    save_params = False

else:
    if TRES_param_override:
        apodization_width = allparam['apod_width']
        apod_type = allparam['apod_type']
        resample = allparam['do_resample']
        resample_factor = allparam['resample_factor']
        shift = allparam["do_shift"]
        pad_test = allparam['do_padzeros']
        padfactor = allparam['pad_factor']
        mean_sub = allparam['do_mean_sub']
        BKGsub = allparam['background_subtract']
        BKGrange = [allparam['background_range_low'],allparam['background_range_high']]
        shift_factor = allparam['shift_factor']
        shift_factor = allparam['shift_factor']
    else:
        apodization_width=0.5
        apod_type="BH"    # "None" "Gauss" "Triangle" "Boxcar" "BH"
        resample="True"
        resample_factor=4
        shift="False"
        pad_test="True"
        padfactor=4
        mean_sub = "True"
        BKGsub = "True"
        BKGrange = [0,12]  #ns  Before t_max
        shift_factor = manual_shift_factor

if save_params:
    params = {"apod_width":apodization_width, "apod_type":apod_type, "do_resample":resample, "resample_factor":resample_factor,
              "do_shift":shift, "do_padzeros":pad_test, "pad_factor":padfactor, "do_mean_sub":mean_sub, "shift_factor":shift_factor,"background_subtract":BKGsub,"background_range_low":BKGrange[0], "background_range_high":BKGrange[1]}
    with open(r"{}\Param_Import_metadata.txt".format(path), 'w+') as ofstream:
        ofstream.write("# Params used in Gemini_MAP_script_CJH.py")
        for param, val in params.items():
            ofstream.write("\n{}:\t{}".format(param, val))

allparam = {'manual_shift_factor':manual_shift_factor,'save_params':save_params, 'rangeval':rangeval, 'intfPlot':intfPlot,'intrfxlims':intrfxlims, 'timeRange':timeRange, 'PLRange':PLRange, 'Gauss_Filter':Gauss_Filter, 'sigmaval':sigmaval, 'AveragePL':AveragePL, 'rangevalPL':rangevalPL,'NormPL':NormPL,'Usemapdata':Usemapdata,'AverageTRPL':AverageTRPL,'rangevalTRPL':rangevalTRPL,'NormTRPL':NormTRPL,'BKGTRPL':BKGTRPL,'TRPLmin_OM':TRPLmin_OM,'overrideTRPLrange':overrideTRPLrange , 'overidexrange':overidexrange ,'composite_legend':composite_legend,"apod_width":apodization_width, "apod_type":apod_type, "do_resample":resample, "resample_factor":resample_factor,
              "do_shift":shift, "do_padzeros":pad_test, "pad_factor":padfactor, "do_mean_sub":mean_sub, "shift_factor":shift_factor,"background_subtract":BKGsub,"background_range_low":BKGrange[0], "background_range_high":BKGrange[1]}
with open(r"{}\\Plot_Params.txt".format(path), 'w+') as ofstream:
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
    pos_data = pos_data - shift_factor    #Taken from shift_factor output in the _INTR analysis script for this data

   #Background Subtract TRPL Curves
    if BKGsub:
        index = [(np.abs(time_data-np.min(BKGrange))).argmin(),(np.abs(time_data-np.max(BKGrange))).argmin()]
        BKGval = np.mean(map_data[:,np.min(index):np.max(index)],axis=1)
        map_data = map_data - np.reshape(BKGval, (len(BKGval), 1))

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

    preFFT_pos, preFFT_map = prep_map(pos_data,map_data,apodization_width,apod_type=apod_type,resample=resample,resample_factor=resample_factor,shift="False",pad_test=pad_test,padfactor=padfactor,mean_sub=mean_sub,plots="True")
    print("Took {} sec".format(time.time() - startTime))

    #Perform FFT
    wave, build_TRES = FFT_intr(preFFT_pos, preFFT_map,plots="False",scale="linear",correct="True")
    wave=wave[::-1]
    build_TRES=np.fliplr(np.array(build_TRES,dtype="float").T)
    print("Took {} sec".format(time.time() - startTime))

# Plot the results (TRES)
indexWL = [(np.abs(wave-np.min(PLRange))).argmin(),(np.abs(wave-np.max(PLRange))).argmin()]
WLPlot=wave[indexWL[0]:indexWL[1]]
indext = [(np.abs(time_data-np.min(timeRange))).argmin(),(np.abs(time_data-np.max(timeRange))).argmin()]
tplot=time_data[indext[0]:indext[1]]

timemesh, wavemesh = np.meshgrid(WLPlot,tplot)
TRESplot=build_TRES[indext[0]:indext[1],indexWL[0]:indexWL[1]]
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
plt.show()

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
    plot_TRES = np.mean(build_TRES,axis=0)
    if NormPL:
        plot_TRES = (plot_TRES-plot_TRES.min())/(plot_TRES.max()-plot_TRES.min())

plt.figure(3, dpi=120)
plt.xlim(min(PLRange),max(PLRange))
plt.ylabel('Counts / a.u.')
plt.xlabel('Wavelength / nm')
plt.yscale('linear')
if AveragePL:
    for i in range(len(plot_TRES)):
        plt.plot(wave,plot_TRES[i],label=str(min(rangevalPL[i])) + " to " + str(max(rangevalPL[i])) + " ns")
    plt.legend()
else:
    plt.plot(wave,plot_TRES)


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
        plt.plot(time_data,integralTRPL[i],label=str(min(rangevalTRPL[i])) + " to " + str(max(rangevalTRPL[i])) + " nm")
    plt.legend()
else:
    plt.plot(time_data,integralTRPL)

#Composite TRES
fig = plt.figure(5,dpi=120)
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