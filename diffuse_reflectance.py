# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 14:51:03 2023

@author: cfai2304
"""
import os
import numpy as np
import matplotlib.pyplot as plt

from interferogram_vis import plot_diffR
from interferogram_functions import kulbeka_munk

pl_path = r"Z:\Data\PL\Y349\7\Y349-7 (vis, intf)_PL_Corrected.csv"
pl = np.loadtxt(pl_path, delimiter=',', skiprows=1)
window = np.logical_and(pl[:,0] >= 500, pl[:,0] <= 900)
pl = pl[window]

std_path = r"Z:\Data\PL\Y349\diffR\std\averaged_PLdata.csv"
std = np.loadtxt(std_path, delimiter=",")

# CSV file with diffuse reflectance spectrum (Sample PL / diffR standard PL)
# two columns - first col wavelength [nm], second col R
path = r"Z:\Data\PL\Y349\diffR\7"
scan_ids = [154342, 155342, ]

for scan_id in scan_ids:
    diffR_path = os.path.join(path, f"{scan_id}", f"{scan_id}_PLdata.csv")
    diffR = np.loadtxt(diffR_path, delimiter=',')
    
    diffR[:, 1] = kulbeka_munk(diffR[:, 1], std[:, 1])
    
    plot_diffR(diffR, pl, title="Y349-7")
    
    
    np.savetxt(os.path.join(path, f"{scan_id}", f"{scan_id}_diffF.csv"),
               diffR.T, delimiter=",", header="Wavelength [nm],F [a.u.]")