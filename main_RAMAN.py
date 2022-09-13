"""
-This program draws vertically stacked XRD plots with references at the bottom.

-All files (.xy files from VESTA export or from powder XRD) should be first placed in the folder 'data\to_be_plotted'.

-Reference files should be named starting with std_ref.

-Curves are normalized and shifted automatically. The plotting order is based on name.

-Guides will be drawn automatically for the strongest 3 peaks for each reference plot.

-After you see a prompt saying 'Input the figure name (Don't forget the suffix .jpg)', enter the name with a .jpg suffix
and hit Enter. The figure is saved in the program folder.

-The 2θ range can be adjusted at line 59.

-Figure size can be changed as needed at line 31.

-Figure dpi can be changed as needed at line 117.

-Additional guides can be drawn as needed at line 111.

version 1 by Ruiquan Yang
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

# 初始化一张空白图，赋予固定且通用的属性
fig, ax = plt.subplots(1,1, figsize=(8.264, 8.264), tight_layout=True)  # The default size of Origin figure (8.264, 11.69) in inch

plot_settings = {"xlabel":'1/cm',
                 "ylabel":'Intensity (a.u.)',
                 "xmin": 50,
                 "xmax": 500,
                 "show_legend":True,
                 "export":False} #input("Input the figure name")}

# Import data
# Place data to be plotted in the folder 'data\to_be_plotted' folder。Using os.walk to traverse this folder.

file_location = r'Z:\Data\Raman Data\Calvin\BaS3 powder\220912\plot_these'
file_ext = ".txt"

# df = pd.DataFrame()  # Initialize an empty DataFrame
pd.set_option('expand_frame_repr', False)
pd.set_option('display.max_rows', None)

# Iterate reading files, organizing data, and plotting
for root, dirs, files in os.walk(file_location):
    # print(files)
    for filename in files:
        # print(filename.endswith(".xy"))
        if filename.endswith(file_ext):
            # Read .xy files with pandas. The data was separated by indefinite separator.
            df_loop = pd.read_csv(
                filepath_or_buffer=os.path.join(file_location, filename),
                sep='\s+',
                names=['1/cm', 'intensity'],
                dtype=float,
                usecols=[0, 1])

            # Range of 2θ to be shown. [10,60] default. Index was reset after the data truncation.
            df_loop = df_loop[(df_loop['1/cm'] >= plot_settings["xmin"]) & (df_loop['1/cm'] <= plot_settings["xmax"])]
            df_loop.reset_index(inplace=True, drop=True)
            # print(df_loop)

            # Find the strongest peak positions within each 1-width span of 2θ.
            if filename.startswith('std_ref'):
                max_intensity = []  # Initialize an empty list for the strongest peak intensity in each 1-width window.
                max_2theta = []  # Initialize an empty list for the strongest 3 peak positions of each std_ref.

                # Split the 2θ range into 1-width windows.
                for pst in range(int(df_loop['1/cm'].min()), int(df_loop['1/cm'].max()), 1):

                    if df_loop['1/cm'].max() > (pst + 1):  # Verify the current window is within the 2θ range.
                        # The strongest intensity in each window.
                        max_intensity.append(df_loop.loc[(pst <= df_loop['1/cm']) & (df_loop['1/cm'] < (pst + 1)), \
                                                         'intensity'].max())
                    else:
                        continue


                # Add the 3 strongest peak positions into a list. Adjust here if you want to plot more ref lines.
                max_2theta = df_loop.loc[
                    df_loop['intensity'].isin(pd.Series(max_intensity).nlargest(20).tolist()), '1/cm']
                # print(max_intensity)

                # Plot the guides.
                

                for i in max_2theta:
                    fig1 = plt.axvline(x=i, ls='--', lw=0.5, alpha=0.3)

            # Normalize the intensity between [0,100]
            df_loop['intensity_norm'] = (df_loop['intensity'] - df_loop['intensity'].min()) / \
                                        (df_loop['intensity'].max() - df_loop['intensity'].min()) * 100

            # Shift the data so that the plots can be shift automatically.
            df_loop['intensity_norm'] += 100 * files.index(filename)
            del df_loop['intensity']

            # df[['2theta_' + filename.rstrip('.xy'), 'intensity_' + filename.rstrip('.xy')]] = df_loop[
            #     ['2theta', 'intensity_norm']]

            # All data from different files are filed into one DataFrame for further manipulation as needed.
            # df = pd.concat([df, df_loop], axis=1)

            # Plot the curve for current file.
            fig1 = plt.plot(
                df_loop['1/cm'], df_loop['intensity_norm'],
                label=filename.rstrip(file_ext),
                lw=0.5)


# other setting for the figure
ax.set_xlabel(plot_settings["xlabel"])
ax.set_ylabel(plot_settings["ylabel"])
if plot_settings["show_legend"]:
    ax.legend(loc='best')  # show legend using file names
#fig1 = plt.axvline(x=26.02, ls='--', lw=0.5, alpha=0.3)  # Add additional guides besides those from the reference.


# Show the figure
# plt.show()

# Save the figure
if plot_settings["export"]:
    plt.savefig(plot_settings["export"], dpi=300, format='jpg')
