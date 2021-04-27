# Interferogram_FFT
Analysis of spectral interferogram data using a fast-Fourier transform (FFT) on the NIREOS GEMENI interferometer

## Installing
Download scripts to your computer. Unzip the calibration files. 
The following libraries are used:
   -  matplotlib
   -  numpy
   -  time
   -  h5py (conda install h5py)
   -  os
   -  scipy
   -  ast
   -  BaselineRemoval (pip install BaselineRemoval)
   -  glob
   -  csv
   -  pandas

## Data Analysis
Data analysis follows the uploaded pdf [Understanding FTIR](https://github.com/HagesLab/Interferogram_FFT/blob/main/Understanding_FTIR.pdf) and the following link
[The FFT in FTIR](https://www.essentialftir.com/fftTutorial.html#:~:text=The%20Fast%20Fourier%20Transform%20(FFT)%20applied%20to%20FTIR%20Data&text=The%20starting%20point%20is%20the,point%2C%20or%20%27ZPD%27)

Recomended Procedure for TRPL MAP data:
1) Read the [pdf](https://github.com/HagesLab/Interferogram_FFT/Theory/blob/main/Understanding_FTIR.pdf) and [link](https://www.essentialftir.com/fftTutorial.html#:~:text=The%20Fast%20Fourier%20Transform%20(FFT)%20applied%20to%20FTIR%20Data&text=The%20starting%20point%20is%20the,point%2C%20or%20%27ZPD%27) referenced above.

2) Run the *Gemini_Averaged_MAP_script_CJH.py* script.
* This will sum the data over all times to maximize your signal for determining optimal FFT parameters. This is more effective than doing it at each time value since the FFT parameters should be the same for all times. Observe the plots that are generated.
* Insert a filepath to the data directory which contains the "..._MAP.txt" "..._POS.txt", and "..._TIME.txt" data
* Pick the appropriate wavelength range for your data, based on the detector used (Vis = 400-1000; NIR = 950-1700).
* Pick appropraite FFT paramters for tha data.
* Run the script and observe the output plots:
    * Plot 1 -> Verify correct background subtraction region. Adjust if needed. The region in red is used to detrmine background, ensure no signal here.
    
    <img src="https://github.com/HagesLab/Interferogram_FFT/blob/main/Readme%20Images/BKGSub.png" width="350">
    
    * Plot 2 -> Verify apodization function. You can zoom in if needed. The blue curve is the measured interferogram, the green curve is the apodization function, anb the orange curve is the new interferogram after applying the apodization function. 
    
    <img src="https://github.com/HagesLab/Interferogram_FFT/blob/main/Readme%20Images/Apod.png" width="350">
    
     * Plot 3 -> Observe the data prep. Should be reflected around peak and padded with desired zeros in the middle.    
     
    <img src="https://github.com/HagesLab/Interferogram_FFT/blob/main/Readme%20Images/Pad.png" width="350">
    
     * Plot 4 & 5 -> Will show the Raw and Phase-Corrected FFT data
     * Plot 6 -> TRPL computed by summing all data (integral TRPL). You can adjust range if needed.
     * Plot 7 -> Resulting integral PL from the FFT
     
    <img src="https://github.com/HagesLab/Interferogram_FFT/blob/main/Readme%20Images/PL%20good.png" width="350"> <img src="https://github.com/HagesLab/Interferogram_FFT/blob/main/Readme%20Images/Poor%20Apod.png" width="373">

3) Observe your resulting PL data and adjust FFT paramters accordingly. 
* This will generally involve the apodization width    

4) If you are satisfied with your FFT, ensure that "save_params = True". This will save the metadata into the "path" directory of our data for use in the MAP script.
5) Run the *Gemini_MAP_script_CJH.py* script.
* The basic way to run this script is to insert the directory for your data which also contains the FFT metadata from the preceding steps. Set "params_from_INTR_metadata = True" to import the metadata for the FFT at each time point.
* If it is the first time running the MAP data analysis, make sure to set *ImportTRES = False*.
* Input a "rangeval" which will trim the data to an upper time to limit the computational load.
* Check plotting ranges in the "Plotting Metadata" section of the code
* Run the script. Adjust "Plotting Metadata" to get the desired plots.
   * To use entire dataset for the TRPL data (for good resolution), set *Usemapdata=True*. This will use the rawdata prior to losing a lot of signal in the FFT, apodiztion, etc. If *Usemapdata=False* and *AverageTRPL = False*, you will plot the full integral TRPL, but it is computed after the FFT (not recommended).
   * You can average the TRPL or PL data over a given range or ranges, input as a list. You must toggle the legend on/off manually when multiple curves are shown.
   * For the standalone TRPL plot, it is posisble to pick a different plot range.
   * It is also possible to fit TRPL data (currently with a single exponential) over a given range. Toggle manually if you want the fit to be in the composite TRES plot. 
<!---      - Here are some examples plots
     
            -   TRES with TRPL fit:
           
          <img src="https://github.com/HagesLab/Interferogram_FFT/blob/main/Readme%20Images/Apod.png" width="350">
          
            -   TRES with averaging TRPL over multiple wavelengths:
            
          <img src="https://github.com/HagesLab/Interferogram_FFT/blob/main/Readme%20Images/Apod.png" width="350">
          
            -   TRES with averaging TRPL over multiple time regions:
           
          <img src="https://github.com/HagesLab/Interferogram_FFT/blob/main/Readme%20Images/Apod.png" width="350">
          
            -   TRES with TRPL fit:
            
          <img src="https://github.com/HagesLab/Interferogram_FFT/blob/main/Readme%20Images/Apod.png" width="350">   -->
