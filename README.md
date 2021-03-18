# Interferogram_FFT
Analysis of spectral interferogram data using a fast-Fourier transform (FFT) on the NIREOS GEMENI interferometer

Data analysis follows the uploaded pdf [Understanding FTIR](https://github.com/HagesLab/Interferogram_FFT/blob/main/Understanding_FTIR.pdf) and the following link
[The FFT in FTIR](https://www.essentialftir.com/fftTutorial.html#:~:text=The%20Fast%20Fourier%20Transform%20(FFT)%20applied%20to%20FTIR%20Data&text=The%20starting%20point%20is%20the,point%2C%20or%20%27ZPD%27)

Recomended Procedure for TRPL MAP data:
1) Read the pdf and link referenced above.
3) Pick the appropriate wavelength range for your data - based on the detector used (Vis = 400-1000; NIR = 950-1700).
4) Run the "Gemini_Averaged_MAP_script_CJH.py" script.
      - This will sum the data over all times to maximize your signal for determining optimal FFT parameters. This is more effective than doing it at each time value since the FFT parameters should be the same for all times. Observe the plots that are generated.
          - Plot 1 -> Verify correct background subtraction region. Adjust if needed.
          - Plot 2 -> Verify apodization function. You can zoom in if needed.
          - Plot 3 -> Observe the data prep. Should be reflected around peak and padded with desired zeros in the middle
          - Plot 4 & 5 -> Will show the Raw and Phase-Corrected FFT data
          - Plot 6 -> TRPL computed by summing all data (integral TRPL). You can adjust range if needed.
          - Plot 7 -> Resulting integral PL from the FFT
 3) Observe your resulting PL data and adjust FFT paramters accordingly. 
      - This will generally involve the apodization width    
