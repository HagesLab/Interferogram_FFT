# -*- coding: utf-8 -*-
"""
Created on Thu May 12 13:56:13 2022

@author: cfai2304
"""
import matplotlib.pyplot as plt
import numpy as np

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
    
    for i in range(len(trpl)):
        ax.plot(times,trpl[i],label=labels[i])
        
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
