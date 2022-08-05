# -*- coding: utf-8 -*-
"""
Created on Sun May 29 23:11:08 2022

@author: ruiqy
"""

import os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, leastsq
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator)

import matplotlib
matplotlib.rcParams.update({'font.size': 11})
matplotlib.rcParams.update({'font.family':'STIXGeneral'})
matplotlib.rcParams.update({'mathtext.fontset':'stix'})

def make_pl_tauc(ax, major_tick=0.5, minor_tick=0.1):
    """
    Input: a tab-limited txt file with two columns: first column wavelength and second column reflectance.
    Specify location of this file using file_location and filename.
    
    Use the trim range to narrow the data and the nominal_fitrange to select a range for linear fitting.



    Returns
    -------
    1. A plot of the normalized kubelka munk parameter (F*E)^2 with the linear fit superposed.
    2. The fitted bandgap (x-intercept) with uncertainty determined by bootstrapping of the nominal_fitrange.

    """
    file_location = r'Z:\Data\PL\Ruiquan BaZrS3 diagnostics\20220804\JN-20 1 hr\151346'
    filename = "diff_R.txt"
    def ylin(x,x_0,k):
        y = k*(x-x_0)
        return y
    
    data = np.loadtxt(os.path.join(file_location, filename), skiprows=1, delimiter='\t')
    wavelength = data[:,0]
    reflectance = data[:,1]

    color = next(ax._get_lines.prop_cycler)['color']
    
    alpha = (1 - reflectance) ** 2 / (2 * reflectance)
    energy = 1239.8 / wavelength

    xcoord = energy
    ycoord = (alpha * energy) ** 2
    
    xcoord = xcoord[::-1]
    ycoord = ycoord[::-1]
    
    trim = [1.4, 2.2]
    low_index, high_index = (np.abs(xcoord-np.min(trim))).argmin() , (np.abs(xcoord-np.max(trim))).argmin()
    ycoord = ycoord[low_index:high_index]
    xcoord = xcoord[low_index:high_index]
    ycoord /= np.amax(ycoord)
    
    ax.plot(
        xcoord, ycoord,
        # label=label,
        c=color
    )
    
    # Linear tauc fitting


    nominal_fitrange = [1.95,2.05]
    
    
        
    low_index, high_index = (np.abs(xcoord-np.min(nominal_fitrange))).argmin() , (np.abs(xcoord-np.max(nominal_fitrange))).argmin()
    ycoordfit = ycoord[low_index:high_index]
    xcoordfit = xcoord[low_index:high_index]
    
    popt, pcov = curve_fit(ylin,xcoordfit,ycoordfit)
    
    #popt, pcov = curve_fit(y,xcoordfit,ycoordfit)
    perr = np.sqrt(np.diag(pcov))
    ypred = ylin(xcoord, *popt)
    resid =  ylin(xcoordfit, *popt) - ycoordfit
    ax.plot(xcoord,ypred,'k--', linewidth=0.5)
    # ax.axvline(nominal_fitrange[0])
    # ax.axvline(nominal_fitrange[1])
    
    r_sq = 1 - (np.sum((ylin(xcoordfit, *popt)-ycoordfit)**2)) / (np.sum((ycoordfit - np.mean(ycoordfit))**2))
    print(filename)
    print(f"R2: {r_sq}")
    print(popt)
    print(perr)
    
    #pfit, perr = fit_bootstrap([2.3,10], xcoordfit, ycoordfit, ylin)

    ps = []
    np.random.seed(42)
    for i in range(1000):
        delta_shift = 0.05
        delta_contract = 0.05
        while True:
            delta_shift = np.random.uniform(-delta_shift, delta_shift)
            delta_contract = np.random.uniform(-delta_contract, delta_contract)
            fitrange = np.array([nominal_fitrange[0] + delta_shift - delta_contract, nominal_fitrange[1] + delta_shift + delta_contract])
            low_index, high_index = (np.abs(xcoord-np.min(fitrange))).argmin() , (np.abs(xcoord-np.max(fitrange))).argmin()
            ycoordfit = ycoord[low_index:high_index]
            xcoordfit = xcoord[low_index:high_index]
        

        
            if not i % 10: print(f"{i}: fitting with {fitrange} - shift {delta_shift}, expand {delta_contract}")
            randomDelta = np.random.normal(0., np.std(resid), len(ycoordfit))
        
            try:
                popt, pcov = curve_fit(ylin,xcoordfit,ycoordfit+randomDelta)
                ps.append(popt)
                break
            except Exception as e:
                print(e)
        
    mean_pfit = np.mean(ps, axis=0)

    # You can choose the confidence interval that you want for your
    # parameter estimates: 
    Nsigma = 1. # 1sigma gets approximately the same as methods above
                # 1sigma corresponds to 68.3% confidence interval
                # 2sigma corresponds to 95.44% confidence interval
    err_pfit = Nsigma * np.std(ps,axis=0) 

    pfit = mean_pfit
    perr = err_pfit
        
    print("\n# Fit parameters and parameter errors from bootstrap method :")
    print("pfit = ", pfit)
    print("perr = ", perr)
     
    ax.set_ylabel(r'$(F \cdot E)^{2}$ (a.u.)')
    ax.set_xlabel(r"Energy (eV)")
    # ax.set_xlim(1.5,4)
    ax.set_ylim(-0.05,1.05)
    
def fit_bootstrap(p0, datax, datay, function, yerr_systematic=0.0):

    errfunc = lambda p, x, y: function(x,*p) - y

    # Fit first time
    pfit, perr = leastsq(errfunc, p0, args=(datax, datay), full_output=0)


    # Get the stdev of the residuals
    residuals = errfunc(pfit, datax, datay)
    sigma_res = np.std(residuals)

    sigma_err_total = np.sqrt(sigma_res**2 + yerr_systematic**2)

    # 100 random data sets are generated and fitted
    ps = []
    for i in range(100):

        randomDelta = np.random.normal(0., sigma_err_total, len(datay))
        randomdataY = datay + randomDelta

        randomfit, randomcov = \
            leastsq(errfunc, p0, args=(datax, randomdataY),\
                             full_output=0)

        ps.append(randomfit) 

    ps = np.array(ps)
    mean_pfit = np.mean(ps, axis=0)

    # You can choose the confidence interval that you want for your
    # parameter estimates: 
    Nsigma = 1. # 1sigma gets approximately the same as methods above
                # 1sigma corresponds to 68.3% confidence interval
                # 2sigma corresponds to 95.44% confidence interval
    err_pfit = Nsigma * np.std(ps,0) 

    pfit_bootstrap = mean_pfit
    perr_bootstrap = err_pfit
    return pfit_bootstrap, perr_bootstrap 


if __name__ == "__main__":
    fig,ax = plt.subplots(figsize=(2, 2), dpi=300)
    make_pl_tauc(ax)


