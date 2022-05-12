# -*- coding: utf-8 -*-
"""
Created on Thu May 12 11:29:26 2022

@author: cfai2304
"""
import numpy as np
import glob
import csv
import os

def import_INTR(path):
    allFiles = []
    b = glob.glob(path + "/*.txt")  # glob. lists the filename and path
    allFiles.extend(b)

    for i in range(len(allFiles)):
        namestest = allFiles[i].split('_')[-1]
        if namestest == r"POS.txt":
            with open(allFiles[i],'r') as j:
                pos_data = np.array(list(csv.reader(j,delimiter="\t")),dtype='float')[0,:]
                samplename = os.path.split(allFiles[i])[-1].split(".")[0]
                samplename = samplename.replace("_POS","")
        elif namestest == r"INTR.txt":
            with open(allFiles[i],'r') as j:
                intr_data = np.array(list(csv.reader(j,delimiter="\t")),dtype='float')[0,:]
    return pos_data, intr_data

def import_MAP(path):
    allFiles = []
    b = glob.glob(path + "/*.txt")  # glob. lists the filename and path
    allFiles.extend(b)


    for i in range(len(allFiles)):
        namestest = allFiles[i].split('_')[-1]
        if namestest == r"MAP.txt":
            with open(allFiles[i],'r') as j:
                map_data = np.array(list(csv.reader(j,delimiter="\t")),dtype='float')
        elif namestest == r"POS.txt":
            with open(allFiles[i],'r') as j:
                pos_data = np.array(list(csv.reader(j,delimiter="\t")),dtype='float')[0,:]
                samplename = os.path.split(allFiles[i])[-1].split(".")[0]
                samplename = samplename.replace("_POS","")
        elif namestest == r"TIME.txt":
            with open(allFiles[i],'r') as j:
                time_data = np.array(list(csv.reader(j,delimiter="\t")),dtype='float')[0,:]

    # Discard the zero columns
    standard=np.std(map_data, axis=0)
    indexes=np.asarray(np.nonzero(standard))
    map_data=map_data[:,indexes[0][0]:indexes[0][-1]]
    time_data =time_data[indexes[0][0]:indexes[0][-1]]/1e3

    return pos_data, time_data, map_data

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