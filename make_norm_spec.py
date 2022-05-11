# -*- coding: utf-8 -*-
"""
Created on Tue May 10 10:14:43 2022

@author: cfai2304
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import interp1d

def load_spectrum(fname, normalize=False):
    supported_delimiters = [",", '\t']
    for delim in supported_delimiters:
        try:
            data = np.loadtxt(fname, delimiter=delim)
            break
        except ValueError:
            continue
        
    waves = data[:,0]
    data = data[:,1]
    if normalize:
        data /= np.amax(data)
        
    return waves, data

def truncate_waves(waves, data, min_w, max_w):
    keep = np.logical_and(waves >= min_w, waves <= max_w)
    waves = waves[keep]
    data = data[keep]
    return waves, data

def interp(x, y, min_x, max_x, new_spacing):
    new_grid = np.arange(min_x, max_x + new_spacing / 2, new_spacing, dtype=float)
    f = interp1d(x, y, axis=0)
    new_y = f(new_grid)
    return new_grid, new_y

def vis(waves, data, label):
    fig, ax = plt.subplots(dpi=120)
    ax.plot(waves, data, label=label)
    ax.legend()
    ax.set_ylabel("Intensity (a.u.)")
    ax.set_xlabel("Wavelength (nm)")
    ax.xaxis.set_minor_locator(MultipleLocator(100))
    ax.grid(axis='both', which='both')
    
def make_norm_spec(lamp_fname, resp_fname, norm_fname=None, show=True):
    
    lamp_waves, lamp_data = load_spectrum(lamp_fname)
    resp_waves, resp_data = load_spectrum(resp_fname)

    lamp_waves, lamp_data = interp(lamp_waves, lamp_data, min_w, max_w, 1)
    resp_waves, resp_data = interp(resp_waves, resp_data, min_w, max_w, 1)

    norm = resp_data / lamp_data
    if norm_fname is not None:
        np.savetxt(norm_fname, np.array([resp_waves, norm]).T, delimiter='\t')
        
    if show:
        vis(resp_waves, resp_data, "Observed intensity (cts)")
        vis(resp_waves, lamp_data, "Avantes known intensity (W/cm^2 nm)")
        vis(resp_waves, norm, "Transfer function")
        
    return resp_waves, norm, lamp_data, resp_data

if __name__ == "__main__":
    min_w = 550
    max_w = 800
    norm_fname = "cuvet_norm_0.txt"
    mode = "a" # "(m)ake" or "(a)pply"
    
    lamp_fname = os.path.join("F:", "1910019-NIR-CCVISNIR.txt")
    resp_fname = os.path.join(r"E:\GEMENI DAQ\NIREOS Complete Example V12_MCS_TimeHarp_32bit Folder\Measurement\20220510\100111", 
                              "100111_PLdata.csv")
    
    data_fname = os.path.join(r"E:\GEMENI DAQ\NIREOS Complete Example V12_MCS_TimeHarp_32bit Folder\Measurement\20220509\170820",
                  "170820_PLdata.csv")

    data_oname = os.path.join(r"E:\GEMENI DAQ\NIREOS Complete Example V12_MCS_TimeHarp_32bit Folder\Measurement\20220509\170820",
                      "170820_PLdata_norm.txt")

    
    if mode == "m":
        make_norm_spec(lamp_fname, resp_fname, norm_fname, show=True)
        
    elif mode == "a":
        norm_waves, norm = load_spectrum(norm_fname)
        
        
        
        data_waves, data = load_spectrum(data_fname)
        
        norm_waves, norm = interp(norm_waves, norm, min_w, max_w, 1)
        data_waves, data = interp(data_waves, data, min_w, max_w, 1)
        
        actual = data / norm
        
        vis(norm_waves, data, "Observed intensity (cts)")
        vis(norm_waves, actual, "Corrected intensity (W/cm^2 nm)")
        vis(norm_waves, norm, "Transfer function")
        
        np.savetxt(data_oname, np.array([norm_waves, actual]).T, delimiter='\t')
    
    else:
        raise ValueError("set mode to m (make a new norm_spec) or a (apply a previously created norm_spec)")