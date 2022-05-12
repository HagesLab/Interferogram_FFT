# -*- coding: utf-8 -*-
"""
Created on Thu May 12 11:29:26 2022

@author: cfai2304
"""

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
