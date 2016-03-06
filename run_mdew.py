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
tio1 = 6708
tio2 = 7046
tio3 = 7084
tio4 = 7120
tio5 = 7046
cah2 = 7042
caoh = 6830
o2 = 6830
n0 = 5895

linelist = [halpha, tio1, tio2, tio3, tio4, tio5, cah2, caoh, o2, n0]
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


def save_ew(fitsimage, line, ew):
    """output: text file with tab separated values of each kic number followed by EqW measurements"""

    fi = fitsimage.find('kic')
    tmp = fitsimage[fi + 3:]
    dot = tmp.find('.')
    kic = int(tmp[:dot])  # gets just the kic number

    measurement = {}
    measurement[line] = ew


    if kic not in ew_dict.keys():
        ew_dict[kic] = measurement
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
            try:
                feature_width, peak = ME.EqWidth(image).find_feature_region(line)
            except Exception:
                pass
            else:
                base1, base2, description, prewave, postwave = ME.EqWidth(image).map_feature(feature_width, peak, line)
                continuum, c1, c2 = ME.EqWidth(image).find_continuum(base1, base2, 5)
                ew = ME.EqWidth(image).measure_ew(feature_width, peak, continuum, c1, c2, prewave, postwave, plot_region = False)

                #print description
                save_ew(image, line, ew)

                if plot_per_image:
                    ME.EqWidth(image).plot_ew(line, peak, base1, base2, continuum, c1, c2)


# TODO: make sure multiple files for one star or stored or appended in list (see todo above)
# TODO: separate ones with magnetic activity
# TODO: add in looking at other bands
# TODO: plot output values (maybe in new file)
# change this line to search a directory that contains directories or files.
ew_per_directory(filepath, plot_per_image=False)

# write output to file:

# set up output file and header

newfile = open('finalEW.txt', 'w')
newfile.write('kic_number' + '\t' + 'ha' + '\t' + 'ha_err' + '\t' + 'tio1'+ '\t' + 'tio1_err'+ '\t' + \
              'tio2' + '\t' + 'tio2_err' + '\t' + 'tio3' + '\t' + 'tio3_err' + '\t' + ' tio4'  + '\t' + \
              'tio4_err' + '\t' + 'tio5' + '\t' + 'tio5_err' + '\t' + 'cah2' + '\t' + 'cah2_err' + '\t' +\
             'caoh' + '\t' + 'caoh_err' + '\t' + 'o2' + '\t' + 'o2_err' + '\t' + 'n0' + '\t' + 'no_err' +'\n')

# triple for loop! find a better way?

for k, v in sorted(ew_dict.items()):
    newfile.write(str(k) + '\t' )
    for l in linelist:
        for lines in v:
            if lines.keys()[0] == l:
                 print lines.keys(), l
                 print 'wooo'
                 newfile.write(str(lines.values()[0]) + '\t')
            else:
                print lines.keys()[0], l
                newfile.write(str(np.nan) + '\t')
                print 'nah'
    newfile.write('\n')
    #newfile.write(str(k) + '\t' +  str(v) + '\n')

newfile.close()
