# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 09:41:38 2022

@author: cfai2
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import curve_fit
from gaussfitting import fit_gaussian, plot_fit_gaussian
matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'font.family':'STIXGeneral'})
matplotlib.rcParams.update({'mathtext.fontset':'stix'})
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']


cases = {"620-RT-Pwr":{"path":r"Z:\Data\PL\JN05 BaZrs3 powders\20220720\620spot_power_at_293K.txt",
                       "title":"620nm spot - PScan at RT", "type":"power", "ub":1e7},
         "620-T":{"path":r"Z:\Data\PL\JN05 BaZrs3 powders\20220720\620spot_temperature_at_ND1.6.txt",
                  "title":"620nm spot - TScan", "type":"temp", "ub":1e6},
         "620-78T-Pwr":{"path":r"Z:\Data\PL\JN05 BaZrs3 powders\20220721\620spot_power_at_78K.txt",
                        "title":"620nm spot - PScan at 78K", "type":"power", "ub":1e7},
         "670-RT-Pwr":{"path":r"Z:\Data\PL\JN05 BaZrs3 powders\20220721\670spot_power_at_293K.txt",
                       "title":"670nm spot - PScan at RT", "type":"power", "ub":1e6},
         "670-T":{"path":r"Z:\Data\PL\JN05 BaZrs3 powders\20220721\670spot_temperature_at_ND1.6.txt",
                  "title":"670nm spot - TScan", "type":"temp", "ub":1e6},
         "670-78T-Pwr":{"path":r"Z:\Data\PL\JN05 BaZrs3 powders\20220721\670spot_power_at_78K.txt",
                        "title":"670nm spot - PScan at 78K", "type":"power", "ub":1e7},
         }

# Select one from above
casename = "620-T"
case = cases[casename]

## Read
all_data = np.loadtxt(case["path"], delimiter='\t')
datasets = {}
maxes = {}
for i in np.arange(0, len(all_data[0]), 2):
    waves = all_data[1:, i]
    pl = all_data[1:, i+1]
    exclevel = all_data[0, i+1]
    
    datasets[exclevel] = (waves, pl)
    
    where_max = np.argmax(pl)
    maxes[exclevel] = (waves[where_max], pl[where_max])

## Plot main
# fig, ax = plt.subplots(1,1, figsize=(5,4), dpi=200)

# for i, nd in enumerate(datasets):
    
#     if i % 2: continue
#     print(nd)
#     ax.plot(datasets[nd][0], datasets[nd][1])
#     ax.scatter(maxes[nd][0], maxes[nd][1], color='red', s=8)
    

    
# ax.set_ylabel("counts [a.u.]")
# ax.set_xlabel("wavelength [nm]")
# ax.set_title(case["title"])


## Fit multiGaussian to each
powers = []
peak_info = []
for nd in datasets:
    print(f"Fitting P={nd}")
    
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    
    scaley = max(datasets[nd][1])
    scalex = 900
    # amp, center, width
    guess = [#0.25, 575 / scalex, 10 / scalex, 
              0.25, 620 / scalex, 30 / scalex, 
              0.25, 670 / scalex, 50 / scalex,
              #25000, 700, 50,
              0.25, 820 / scalex, 30 / scalex, 
              ]  # I guess there are 3 peaks, 2 are clear, but between them there seems to be another one, based on the change in slope smoothness there
    if len(guess) % 3 != 0:
        raise ValueError("Incorrect num of fit params")
    
    out = fit_gaussian(datasets[nd][0] / scalex, datasets[nd][1] / scaley, guess)
    out[0][::3] *= scaley
    out[0][1::3] *= scalex
    out[0][2::3] *= scalex
    out[2][:] *= scaley
    title = f"P={nd} uW" if case['type'] == 'power' else f"T={nd} K"
    #title = ""
    plot_fit_gaussian(datasets[nd][0], datasets[nd][1], out, colors=colors, title=title)
    peaks = out[0]
    peaks = [peaks[3*i:3*i+3] for i in range(len(peaks) // 3)]
    peaks.sort(key=lambda x:x[1])
    powers.append(nd)
    peak_info.append(list(peaks))
    
peak_info = np.array(peak_info)

peak_info = [peak_info[:,i,:] for i in range(len(peak_info[0]))]
def lin(x, k, C):
    return C + k*x
kb = 8.617e-2
def EA2(x, E1, C1, E2, C2):
    """ For DAP-style transitions """
    return np.log(1 / (1 + C1 * x**1.5 * np.exp(-E1 / (kb*x)) + C2 * x **1.5 * np.exp(-E2 / (kb*x))))

def EA(x, E1, C1):
    """ For basic transitions """
    return np.log(1 / (1 + C1 * x**1.5 * np.exp(-E1 / (kb*x))))


def make_llplplot(x, y, **kwargs):
    fig = kwargs.get('fig', None)
    ax = kwargs.get('ax', None)
    color = kwargs.get('color', 'black')
    scantype = kwargs.get("scantype", None)
    offset = kwargs.get('offset', 1)
    ub = kwargs.get('ub', 1e6)
    
    if scantype == "power":
        popt, pcov = curve_fit(lin,np.log(x),np.log(y))
        perr = np.sqrt(np.diag(pcov))
        
        C = popt[1]
        k = popt[0]
        fitted = np.exp(C) * x ** k
        
        label = "{:.2f}".format(k) + r"$\ \pm\ $" + r"{:.2f}".format(perr[0])
        print("Idep", label)
        
    elif scantype == "temp":
        scaley = np.amax(y)
        popt, pcov = curve_fit(EA, x, np.log(y / scaley), p0=[0.6, 2.6], bounds=(0, np.inf))
        perr = np.sqrt(np.diag(pcov))
        
        E1 = popt[0]
        #E2 = popt[2]
        fitted = np.exp(EA(x, *popt)) * scaley
        
        label1 = "{:.2f}".format(E1) + r"$\ \pm\ $" + r"{:.2f}".format(perr[0])
        #label2 = "{:.2f}".format(E2) + r"$\ \pm\ $" + r"{:.2f}".format(perr[2])
        print("Idep", label1, "\n", "meV")
        print(popt)

    
    if fig is None or ax is None:
        fig, ax = plt.subplots(1,1,figsize=(5,5), dpi=200)
    
    if scantype == 'power':
        ax.scatter(x, y * offset, color=color)
        ax.plot(x, fitted * offset, color=color, linestyle='dashed', linewidth=1)
        
    elif scantype == "temp":
        ax.scatter(1000/x, y * offset, color=color)
        ax.plot(1000/x, fitted * offset, color=color, linestyle='dashed', linewidth=1)
        
    ax.set_ylabel("PL-counts [a.u.]")
    if scantype == 'power':
        ax.set_xlabel(r"excitation intensity [$\mathrm{\mu}$W]")
        ax.set_xscale('log')
    elif scantype == 'temp':
        ax.set_xlabel("1000 / T [1/K]")
    ax.set_yscale('log')
    ax.set_ylim(1e3, ub)
    ax.set_title(r"I_$\mathrm{PL}$ dependence")

    
def make_evplot(x,y, **kwargs):
    fig = kwargs.get('fig', None)
    ax = kwargs.get('ax', None)
    ax2 = kwargs.get('ax2', None)
    color = kwargs.get('color', 'black')
    scantype = kwargs.get("scantype", None)
    if fig is None or ax is None:
        fig, ax = plt.subplots(1,1,figsize=(5,5), dpi=200)
        
    if scantype == "temp":
        popt, pcov = curve_fit(lin,x,y)
        perr = np.sqrt(np.diag(pcov))
        
        C = popt[1]
        k = popt[0]
        fitted = x * k + C
        
        label = "{:.2f}".format(k*1e3) + r"$\ \pm\ $" + r"{:.2f}".format(perr[0]*1e3)
        print("evDep " , label, " meV / K" )
        ax.plot(x, fitted, color=color, linestyle='dashed', linewidth=1, label=label)
        
    elif scantype == "power":
        popt, pcov = curve_fit(lin,np.log10(x),y)
        perr = np.sqrt(np.diag(pcov))
        
        C = popt[1]
        k = popt[0]
        fitted = np.log10(x) * k + C
        
        label = "{:.2f}".format(k*1e3) + r"$\ \pm\ $" + r"{:.2f}".format(perr[0]*1e3)
        print("eVdep ", label, "meV / decade")
        ax.plot(x, fitted, color=color, linestyle='dashed', linewidth=1, label=label)
        
    ax.scatter(x,y, color=color)
    print("avg energy ", np.mean(y))
    
    ax.set_ylabel("peak energy [eV]")
    ax.set_ylim(1240/900, 1240/540)
    
    ax2.set_ylabel("peak wavelength [nm]")
    
    ax2.set_yticks([(1240/i) for i in np.arange(500,901,100)])
    ax2.set_yticklabels([i for i in np.arange(500,901,100)])
    ax2.set_ylim(ax.get_ylim())
    
    if scantype == 'power':
        ax.set_xlabel(r"excitation intensity [$\mathrm{\mu}$W]")
        ax.set_xscale('log')
    elif scantype == 'temp':
        ax.set_xlabel("temperature [K]")
    ax.set_title("Peak energy dependence")
    #ax.legend(loc='lower right', title=r"$\beta \times 10^4$", framealpha=0.5)
  
def make_fbplot(x,T, y, **kwargs):
    """ x = energy [eV]
    T = temperature [K]
    y = PL yield [cts] """
    fig = kwargs.get('fig', None)
    ax = kwargs.get('ax', None)
    ax2 = kwargs.get('ax2', None)
    color = kwargs.get('color', 'black')
    scantype = kwargs.get("scantype", None)
    ub = kwargs.get('ub', 1e6)
    
    def db(x, dE):
        return x**2 * (x-dE)**0.5 * np.exp((dE-x) / (kb*T))
    
    if fig is None or ax is None:
        fig, ax = plt.subplots(1,1,figsize=(5,5), dpi=200)
        

    popt, pcov = curve_fit(db,x*1e3,np.log(y),p0=[0.1])
    perr = np.sqrt(np.diag(pcov))
    
    dE = popt[0]
    fitted = np.exp(db(x, dE))
    
    label = "{:.2f}".format(dE) + r"$\ \pm\ $" + r"{:.2f}".format(perr[0])
    print("dbcheck " , label, " meV" )
    ax.plot(x, fitted, color=color, linestyle='dashed', linewidth=1, label=label)
        
    ax.scatter(x,y, color=color)
    
    ax.set_xlabel("peak energy [eV]")
    ax.set_xlim(1240/900, 1240/540)
    
    ax.set_ylabel("PL-counts [a.u.]")
    ax.set_yscale('log')
    ax.set_ylim(1e3, ub)
    ax.set_title(r"I_$\mathrm{PL}$ dependence")

## Unzip
x = np.array(powers)

figll, axll = plt.subplots(1,1,figsize=(4.5,4.5), dpi=200)
figev, axev = plt.subplots(1,1,figsize=(4.5,4.5), dpi=200)
#figfb, axfb = plt.subplots(1,1,figsize=(4.5,4.5), dpi=200)
ax2 = axev.twinx()

offsets = [10,1,0.3]

# Reassign peaks
x = [x,x,x]
if casename == "670-78T-Pwr":
    #peak_info[x, which_peak, A/mean/width]
    peak_info[1][4:, :] = peak_info[0][4:,:]
    peak_info[1][:2, :] = peak_info[0][:2, :]
    #peak_info = peak_info[:,1:,:]
    x = [x,x,x]
    start = 1
    
elif casename == "670-T":
   
    
    peak_info[2][5, :] = peak_info[1][5,:]
    peak_info[1][5,:] = peak_info[0][5,:]
    peak_info[1][6,:] = peak_info[0][6,:]
    peak_info[1][8,:] = peak_info[0][8,:]
    peak_info[1][9,:] = peak_info[0][9,:]
        
    start = 1
else:
    start = 0
    
for i in range(start, len(peak_info)):
    #figll, axll = plt.subplots(1,1,figsize=(4.5,4.5), dpi=200)
    y_amp = peak_info[i][:,0]
    y_wave = peak_info[i][:,1]
    make_llplplot(x[i], y_amp, fig=figll, ax=axll, color=colors[i], scantype=case['type'],
                  offset = offsets[i], ub=case['ub'])
    make_evplot(x[i], 1240 / y_wave, fig=figev, ax=axev, ax2=ax2, color=colors[i], scantype=case['type'])
    
    # if case['type'] == 'temp':
    #     T = x
    # elif "78T" in casename:
    #     T = 78
    # else:
    #     T = 298
    # make_fbplot(1240 / y_wave, T, y_amp, fig=figfb, ax=axfb, color=colors[i], scantype=case['type'],
    #               offset = offsets[i], ub=case['ub'])