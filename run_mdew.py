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

# global variables = which lines to measure
filepath = '/home/tessa/astronomy/mdwarf/highres_data'

halpha = 6562.801
tio5 = 0

linelist = [halpha]

# set up output file and header
ew_dict = {}
newfile = open('finalEW.txt', 'w')
newfile.write('kic_number' + '\t' + 'ha_ew' + '\t' + 'ha_ew_err' + '\t' + 'tio5' + '\t' + 'tio5_err''\n')


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


def save_ew(fitsimage, ew):
    """output: text file with tab separated values of each kic number followed by EqW measurements"""

    fi = fitsimage.find('kic')
    tmp = fitsimage[fi + 3:]
    dot = tmp.find('.')
    kic = int(tmp[:dot])  # gets just the kic number

    if kic not in ew_dict.keys():
        ew_dict[kic] = ew
    else:
        stack = np.hstack((kic, ew))
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
            try:
                feature_width, peak = ME.EqWidth(image).find_feature_region(line)
            except Exception:
                pass
            else:
                base1, base2, description = ME.EqWidth(image).map_feature(feature_width, peak, line)
                continuum, c1, c2 = ME.EqWidth(image).find_continuum(base1, base2, 5)
                ew = ME.EqWidth(image).measure_ew(feature_width, peak, continuum, c1, c2)

                print description
                save_ew(image, ew)

                if plot_per_image:
                    ME.EqWidth(image).plot_ew(line, peak, base1, base2, continuum, c1, c2)


# TODO: make sure multiple files for one star or stored or appended in list (see todo above)
# TODO: separate ones with magnetic activity
# TODO: add in looking at other bands
# TODO: plot output values (maybe in new file)

# change this line to search a directory that contains directories or files.
ew_per_directory(filepath, plot_per_image=False)

# write output to file:
for k, v in sorted(ew_dict.items()):
    newfile.write(str(k) + '\t' + str(v[0]) + '\t' + str(v[1]) + '\n')
newfile.close()
