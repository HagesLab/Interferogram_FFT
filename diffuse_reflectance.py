# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 14:51:03 2023

@author: cfai2304
"""

import numpy as np
import matplotlib.pyplot as plt

# CSV file with diffuse reflectance spectrum (Sample PL / diffR standard PL)
# two columns - first col wavelength [nm], second col R
diffR_path = r"E:\GEMENI DAQ\NIREOS Complete Example V12_MCS_TimeHarp_32bit Folder\Measurement\20230120\bzs diffR\143556_diffR.csv"

diffR = np.loadtxt(diffR_path, delimiter=',')

# Convert spectrum from nm to eV
E = 1240 / diffR[:, 0]
R = diffR[:,1]

F = (1 - R)**2 / (2*R) # Kubelka Munk term
#F *= 1240 / E**2 # jacobian transform from wavelength to energy integral

pl_path = r"E:\GEMENI DAQ\NIREOS Complete Example V12_MCS_TimeHarp_32bit Folder\Measurement\20230120\bzs 532+550 lots washing\bzs 532+550 lots washing_PLdata.csv"
pl = np.loadtxt(pl_path, delimiter=',')

bkg_path = r"E:\GEMENI DAQ\NIREOS Complete Example V12_MCS_TimeHarp_32bit Folder\Measurement\20230120\mirror 532+550\mirror 532+550_PLdata.csv"
bkg = np.loadtxt(bkg_path, delimiter=',')

pl[:,1] -= 3.3*bkg[:,1]

window = np.logical_and(pl[:,0] >= 560, pl[:,0] <= 850)
pl = pl[window]

fig, ax = plt.subplots(1,1,figsize=(2.5,2.5), dpi=240)
ax2 = ax.twinx()
ax.plot(diffR[:,0], F, color='red')
ax2.plot(pl[:,0], pl[:,1])
#ax.set_xlim(1240/900, 1240/550)
ax.set_xlim(550,900)
ax.set_xticks([550] + [100*i for i in range(6, 10)])
ax.set_xlabel("wavelength [nm]")
# ax.set_ylabel("R")
#ax.set_xlabel("energy [eV]")
ax.set_ylabel("Kubelka Munk F", color='red')
ax2.set_ylabel("PL [cts]", color='steelblue')
ax.set_title("BaZrS3")
ax.set_ylim(0, 6)
