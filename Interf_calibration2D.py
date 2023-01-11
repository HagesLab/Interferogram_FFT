import numpy as np
import pandas as pd
import Get_calibrated_position_axis
import os

'''interf_calibration2D calibrate the measurement map and the position axis to be ready to be the input of the Fourier transform
    INPUTS:
    - interf2d: interference map each column is an interferogram
    - pos_axis: position axis correspondent to the measured map
    OUTPUTS:
    - calib_interf2d: calibrated interference map
    - calib_pos_axis: calibrated position axis
'''


def interf_calibration2D(interf2d, pos_axis, pre2023=False):
    calib_pos_axis, left_index, right_index=Get_calibrated_position_axis.get_calibrated_position_axis(pos_axis, pre2023)
    calib_interf2d = interf2d[left_index-1:right_index,:]
    return calib_interf2d, calib_pos_axis
