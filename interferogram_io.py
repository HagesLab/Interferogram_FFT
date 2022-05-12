# -*- coding: utf-8 -*-
"""
Created on Thu May 12 11:29:26 2022

@author: cfai2304
"""
import numpy as np

def fetch_metadata(path):  #  "{}\\Param_Import_metadata.txt"  or "{}\\Plot_Params.txt"
    with open(path, "r") as ifstream:
        param_values_dict = {}
        for line in ifstream:
            if "#" in line: continue

            else:
                param = line[0:line.find(':')]
                new_value = line[line.find('\t') + 1:].strip('\n')

                try:
                    if "." in new_value:
                        param_values_dict[param] = float(new_value)
                    else:
                        param_values_dict[param] = int(new_value)
                except:
                    param_values_dict[param] = str(new_value)

    return param_values_dict

def save_metadata(path, params, from_="INTR"):
    with open(path, 'w+') as ofstream:
        ofstream.write("# Params used in {} script".format(from_))
        for param, val in params.items():
            ofstream.write("\n{}:\t{}".format(param, val))
            
    return

def save_PL(path, wave_list, apods, spectra_list):
    header = []
    data = np.empty((len(wave_list[0]), len(apods) * 2))
    for i, apod in enumerate(apods):
        data[:, 2*i] = wave_list[i]
        data[:, 2*i+1] = spectra_list[i].flatten()
        header.append("Wavelength [nm] apod={}".format(apod))
        header.append("PL [counts] apod={}".format(apod))
    np.savetxt(path, data, delimiter=',', header=",".join(header))
    return