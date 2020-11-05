import numpy as np
import pandas as pd
import os
import math
from scipy import interpolate
from scipy import signal
import Get_real_position_axis
from scipy import signal

'''get_calibrated_position_axis computes the calibrated position axis based on the encoder position axis.
    INPUTS:
    - position_axis: position axis correspondent to the interferogram
    OUTPUTS:
    - calibrated_position_axis: position axis calibrated on the encoder od the positioner
'''

def get_calibrated_position_axis(position_axis):

    # Load the calibrations file
    items = os.listdir(".")

    for names in items:
        if names.endswith("parameters_int.txt"):
            filename = names

    ref = pd.read_csv(filename, sep="\t", header=None)
    first_row = (ref.iloc[0])
    second_row = (ref.iloc[1])
    position_ref = first_row.to_numpy(dtype='float64')
    amplitude_ref = second_row.to_numpy(dtype='float64')

    position_axis = np.asarray(position_axis).squeeze()

    # Ratio between position axises
    factor = np.ceil((np.mean(np.diff(position_axis))) / (np.mean(np.diff(position_ref))))

    # Oversample of the position axis
    f = interpolate.interp1d(position_axis, position_axis, kind='cubic')
    oversmp_pos = f(np.linspace(position_axis[0], position_axis[-1], int(factor) * np.shape(position_axis)[0]))

    # Interpolated reference on the position axis
    f1 = interpolate.interp1d(position_ref, amplitude_ref, kind='cubic', fill_value="extrapolate")
    ref_interp = f1(oversmp_pos)

    # Find the peaks of the interpolated reference and it gets their indexes
    locs = signal.find_peaks(ref_interp)
    indexes=locs[0]
    left_index=indexes[0]
    right_index=indexes[-1]

    # Compute position axis calibrated on the encoder
    ref_interp_cut = ref_interp[left_index: right_index+1]

    calib_pos_axis_overs = Get_real_position_axis.get_real_position_axis(ref_interp_cut) * (
            position_axis[-1] - position_axis[0]) + position_axis[0]

    calibrated_position_axis=calib_pos_axis_overs[0:-1:int(factor)]

    # Compute the indexes factorized
    left_index = int(np.ceil(indexes[0]/ factor))
    right_index = int(left_index + len(calibrated_position_axis)-1)

    return calibrated_position_axis, left_index, right_index
