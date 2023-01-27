# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 15:58:21 2023

This script subtracts a background interferogram from a base interferogram.

@author: cfai2304
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
from scipy.optimize import least_squares

# Load a base interferogram - an INTR.txt file and a POS.txt file
intr_dir = r"E:\GEMENI DAQ\NIREOS Complete Example V12_MCS_TimeHarp_32bit Folder\Measurement\20230120\bzs 532+550 lots washing"
intr_path = r"bzs_lots_washing_INTR.txt"
pos_path = r"bzs_lots_washing_POS.txt"

# Load a background interferogram - see "interf_vis_noise.txt" and 
# "interf_ir_noise.txt" for examples.
# These can be created by taking a dark measurement (i.e. no laser) on a detector,
# then applying a smoothing function like a savgol filter.
bkg = np.loadtxt("interf_vis_noise.txt")

pos = np.loadtxt(os.path.join(intr_dir, pos_path))
intr = np.loadtxt(os.path.join(intr_dir, intr_path))

# Interpolate the background to the base
bkg_int = griddata(bkg[:,0], bkg[:,1], pos, method='cubic')

fig, ax = plt.subplots(1, 1, dpi=200)
ax.plot(pos, intr, label="Base interferogram")

# Find best multipe of background to subtract
def resid(A, y, yhat):
    return np.sum((y - A*yhat)**2)

sol = least_squares(resid, 1, args=(intr, bkg_int))
best_noise_mult = sol.x[0]
print("Best bkg noise multiplier: {}".format(best_noise_mult))

intr -= best_noise_mult*bkg_int
ax.plot(pos,best_noise_mult*bkg_int, label="Matched background")
ax.set_title("Background matching attempt")

fig, ax = plt.subplots(1, 1, dpi=200)
ax.plot(pos, intr)
ax.set_title("Background subtracted")

# Export bkg subtracted
np.savetxt(os.path.join(intr_dir, "bzs_lots_washing_INTR_new.txt"), intr.reshape((1, len(intr))), delimiter='\t')