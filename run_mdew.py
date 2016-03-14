__author__ = 'tdwilkinson'
"""
Runs Mdwarf_EqW.py on a directory full of images!
Change filepath to edit where it runs. Edit linelist to change output file
Change the call at the bottom to plot_per_image = True if wanted

output: finalEW.txt
"""

import Mdwarf_EqW as ME
import os
import numpy as np
import pandas as pd

# global variables = which lines to measure
filepath = '/home/tessa/astronomy/mdwarf/highres_data'

halpha = [6562.801]
tio1 = [6703, 6708, 6718, 6723]
tio2 = [7043, 7046, 7058, 7061]
tio3 = [7079, 7084, 7092, 7097]
tio4 = [7115, 7120, 7130, 7135]
tio5 = [7126, 7135, 7042, 7046]
cah2 = [6814, 6846, 7042, 7046]
cah3 = [6960, 6990, 7042, 7046]
caoh = [6345, 6354, 6230, 6240]
o2 = [6830, 6850, 6868, 6881]
n0 = [5895]
na1 = [5865, 5880]
na2 = [5910, 5925]

# hardcoded output names for now:
linelist = [halpha, tio1, tio2, tio3, tio4, tio5, cah2, cah3, caoh, o2, n0]
outcolumns = ['ha', 'tio1', 'tio2', 'tio3', 'tio4', 'tio5', 'cah2', 'cah3', 'caoh', 'o2', 'n0']
# chi = 10**([0, -3.93438, -4.01516, -4.13192, -4.19592, -4.56243,-4.75570, -5.28066, -5.21965, -5.41719])


ew_dict = {}

def get_wavelength_calibrated_fits_files(input_directory):
    """
        Give the directory containing multiple folders of observing dates and data. This definition will go through
        and look for the wavelength calibrated fits images and save their path locations to a list.
        output = list of paths to fits files
        """
    check = []
    fits_path = []
    for root, dirs, files in os.walk(input_directory, topdown=False):
        for name in files:
            if (os.path.join(root, name)) not in check and name.startswith('dkic') and name.endswith('.fits'):
                check.append(os.path.join(root, name))
                fits_path.append(os.path.join(root, name))
    return fits_path


def save_ew(fitsimage, line, to_save):
    """output: text file with tab separated values of each kic number followed by EqW measurements"""

    fi = fitsimage.find('kic')
    tmp = fitsimage[fi + 3:]
    dot = tmp.find('.')
    kic = int(tmp[:dot])  # gets just the kic number

    measurement = {}
    measurement[line] = to_save

    if kic not in ew_dict.keys():
        ew_dict[kic] = measurement
    else:
        stack = np.hstack((ew_dict[kic], measurement))
        ew_dict[kic] = stack


def ew_per_directory(parent_directory, plot_per_image=True):
    """
    Runs the class EqWidth though Mdwarf_EqW.py

    output: a text file in the output_directory for each .fits file starting with dkic
            found in the parent directory
    """
    path_list = get_wavelength_calibrated_fits_files(parent_directory)
    for image in path_list:
        for line in linelist:
            if len(line) == 1:
                try:
                    feature_width, peak = ME.EqWidth(image).find_feature_region(line[0])
                    base1, base2, description, prewave, postwave = ME.EqWidth(image).map_feature(feature_width, peak, line[0])
                    continuum, c1, c2 = ME.EqWidth(image).find_continuum(base1, base2, 5)
                    ew, err = ME.EqWidth(image).measure_ew(feature_width, peak, continuum, c1, c2, prewave, postwave, base1, base2, plot_region = False)
                    save_ew (image, line[0], [ew, err])
                except Exception:
                    pass
            else:
                num = ME.EqWidth(image).avg_flux_per_region(line[0], line[1])
                den = ME.EqWidth(image).avg_flux_per_region(line[2], line[3])
                if num != None != den:
                    ratio = num/den
                    save_ew(image, line[0], [ratio])

                if plot_per_image:
                    ME.EqWidth(image).plot_ew(line, peak, base1, base2, continuum, c1, c2)

ew_per_directory(filepath, plot_per_image=False)

# write output to file:
df = pd.DataFrame(index = sorted(ew_dict.keys()), columns = outcolumns)
for k, v in sorted(ew_dict.items()):
    for lines in v:
        for value, measurements in lines.items():
            for n in xrange(len(outcolumns)):
                if value == linelist[n][0]:
                    df.loc[k, outcolumns[n]] = measurements
df.to_csv('Mdwarf_out.txt', sep = '\t')

