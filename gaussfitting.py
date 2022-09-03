# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 19:45:03 2022

@author: cfai2
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def gaussian(x, height, center, width, offset=0):
    return height*np.exp(-(x - center)**2/(2*width**2)) + offset

def n_gaussians(x, *args, **kwargs):
    offset = kwargs.get('offset', 0)
    yg = 0
    for i in range(len(args) // 3):
        yg += gaussian(x, args[3*i], args[3*i+1], args[3*i+2], offset)
        
    return yg

def fit_gaussian(x, y, guess):
    success = False
    try:
        optim, pcov = curve_fit(n_gaussians, x , y, p0=guess[:], bounds=(0,np.inf))
        success = True
    except RuntimeError:
        print("Warning: fitting failed")
        success = False
        return None
        
    if success:
        fitted_y  = n_gaussians(x, *optim)
        print(np.sum((fitted_y - y)**2))
        
    return optim, pcov, fitted_y
        
def plot_fit_gaussian(x, y, out, **kwargs):
    colors = kwargs.get("colors", None)
    title = kwargs.get("title", "")
    
    fig, ax = plt.subplots(1,1, figsize=(5,4), dpi=200)

    ax.plot(x, y,lw=5, c='red', label='measurement')
    if out is not None:
        optim, pcov, fitted_y = out
        ax.plot(x, fitted_y,
            lw=3, c='black', linestyle='dashed', label='3-gaussian fit')
    
        assert len(optim) % 3 == 0, "Incorrect num of fit params"
        
        gau = []
        for i in range(len(optim) // 3):
            gau.append(optim[3*i:3*i+3])
            
        gau.sort(key=lambda x:x[1])
        for i, g in enumerate(gau):
            yg = gaussian(x, *g)
            
            #label = "A: {:.2f}, M: {:.2f}, W: {:.2f}".format(*g)
            label = None
            ax.plot(x, yg, color=colors[i], label=label)
            ax.fill_between(x, 0, yg, facecolor=colors[i], alpha=0.5)
    
    ax.legend()
    ax.set_title(title)
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("counts [a.u.]")