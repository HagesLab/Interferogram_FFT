import numpy as np

'''get_real_position_axis calculate the calibrated axis from a reference interferogram with high temporal coherence
    INPUTs:
    - reference: is the reference interferogram
    OUTPUTSs:
    - position_axis: is the calibrated position axis
'''


def normalize(a):
    b = (a-min(a))/(max(a)-min(a))
    return b


def get_real_position_axis(reference):
    fft_ref = np.fft.fft(reference)
    fft_ref[0: int(np.floor(len(reference)/2) - 1)] = np.zeros(int(np.floor(len(fft_ref) / 2) - 1))
    position_axis = normalize(np.unwrap(-np.angle(np.fft.ifft(fft_ref))))

    return position_axis
