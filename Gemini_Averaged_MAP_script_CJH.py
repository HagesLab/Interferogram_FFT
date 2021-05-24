# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 11:37:33 2020

@author: Chuck
"""


#New Time Averaged data from MAP
index=(np.abs(time_data)).argmin()
#AVG_map_data = np.mean(map_data,axis=1)
AVG_map_data = map_data[:,index]


preFFT_pos, preFFT_data, shiftfactor, baseline_fit = prep_interferogram(pos_data,AVG_map_data,apodization_width,apod_type=apod_type,resample=resample,resample_factor=resample_factor,shift=shift,pad_test=pad_test,padfactor=padfactor,mean_sub=mean_sub,plots=True,pltzoom=pltzoomstate,zoom_range=pltzoomrange,baseline_sub_state=baseline_sub_state)
wave, FFT_intr_trim = FFT_intr(preFFT_pos,preFFT_data,plots=True,scale="linear",correct=True)


#Plot Full TRPL
integralTRPL = np.sum(map_data,axis=0)
plt.figure(4, dpi=120)
plt.plot(time_data,integralTRPL)
plt.ylim(np.max(integralTRPL)*TRPLmin_OM,2*np.max(integralTRPL))
#plt.xlabel('Time / ns')
#plt.ylabel('Counts / a.u.')
#plt.title("Integral TRPL")
plt.yscale('log')

#Plot Full PL
plt.figure(5, dpi=120)
#plt.ylabel("Counts / a.u.")
#plt.xlabel("Wavelength / nm")
#plt.title("Average PL")
plt.plot(wave,FFT_intr_trim)
plt.xlim(start_wave,end_wave)
plt.yscale('linear')
if save_params:
    outputfilename_meta = path + "\\" + os.path.split(path)[-1] + '_FFTmetadata.txt'
    params = {"apod_width":apodization_width, "apod_type":apod_type, "do_resample":resample, "resample_factor":resample_factor,
              "do_shift":shift, "do_padzeros":pad_test, "pad_factor":padfactor, "do_baseline_sub":baseline_sub_state, "do_mean_sub":mean_sub,"shift_factor":shiftfactor,"background_subtract":BKGsub,"background_range_low":np.min(BKGrange), "background_range_high":np.max(BKGrange)}
    with open(outputfilename_meta.format(path), 'w+') as ofstream:
        ofstream.write("# Params used in Gemini_Averaged_MAP_script_CJH.py")
        for param, val in params.items():
            ofstream.write("\n{}:\t{}".format(param, val))

