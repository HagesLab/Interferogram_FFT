# -*- coding: utf-8 -*-
"""
Created on Thu May 12 13:56:13 2022

@author: cfai2304
"""
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors
from matplotlib.ticker import LogLocator
import numpy as np
from scipy import ndimage

def plot_PL_spectrum(waves, spectra, labels, start_wave, end_wave, export=None, interval=None):
    fig, ax = plt.subplots(dpi=120)
    ax.set_ylabel("Counts / a.u.")
    ax.set_xlabel("Wavelength / nm")
    ax.set_title("Average PL")
    ax.grid()
    if interval is not None: 
        ax.set_xticks(np.arange(start_wave,end_wave+interval / 2,interval))
        
    if isinstance(waves, list):
        for i in range(len(waves)):
            ax.plot(waves[i],spectra[i],label=labels[i])
    
    else:
        if isinstance(spectra, list) or (isinstance(spectra, np.ndarray) and spectra.ndim == 2):
            for i in range(len(spectra)):
                ax.plot(waves,spectra[i],label=labels[i])

        else:
            ax.plot(waves,spectra)
    ax.set_xlim(start_wave,end_wave)
    ax.set_yscale('linear')
    
    ax.legend()
    
    if export is not None:
        fig.savefig(export)

def plot_TRPL_decay(times, trpl, min_OM, labels=None, start_time=None, end_time=None, 
                    fit=None, export=None, interval=None):
    fig, ax = plt.subplots(dpi=120)
    ax.plot(times,trpl[0],label=None)
    for i in range(len(fit[0])):
        # ax.plot(times,trpl[i],label=labels[i])
        
        if fit is not None:
            ax.plot(fit[0][i], fit[1][i], 'k--', label=''.join(fit[2][i]))
         
    if start_time is not None and end_time is not None:
        ax.set_xlim((start_time, end_time))
         
    ax.set_ylim(np.max(trpl)*min_OM,2*np.max(trpl))
    ax.set_xlabel('Time / ns')
    ax.set_ylabel('Counts / a.u.')
    ax.set_title("Integral TRPL")
    ax.set_yscale('log')
    ax.legend()
    if export is not None:
        fig.savefig(export)

def plot_TRES(timemesh, wavemesh, TRESplot, ax=None, fig=None, Gauss_Filter=False, sigma=0):
    min_value = np.amin(TRESplot)
    max_value = np.amax(TRESplot)
    
    if fig is None:
        fig, ax = plt.subplots(dpi=120)
    if Gauss_Filter:
        TRESplot = ndimage.gaussian_filter(TRESplot, sigma=sigma)
    norm= matplotlib.colors.LogNorm(vmin=min_value, vmax=max_value)
    levels = np.logspace(np.log10(min_value),np.log10(max_value),num=50)
    cs = ax.contourf(timemesh,wavemesh,TRESplot,levels=levels,norm=norm, cmap='plasma')
    cbar = fig.colorbar(cs)
    cbar.ax.yaxis.set_major_locator(LogLocator())
    cbar.set_ticks(cbar.ax.yaxis.get_major_locator().tick_values(min_value, max_value))
    ax.set_ylabel('Time / ns')
    ax.set_xlabel('Wavelength / nm')
