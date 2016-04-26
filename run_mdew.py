__author__ = 'tdwilkinson'
"""
Runs Mdwarf_EqW.py recursively by searching for spectra in folders and measuring their properties.

# will have to change the file name parts to be more general

"""

import Mdwarf_EqW as ME
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# global variables = which lines to measure
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
outcolumns = ['ha', 'tio1', 'tio2', 'tio3', 'tio4', 'tio5', 'cah2', 'cah3', 'caoh', 'o2', 'n0', 'SpT']
ew_dict = {}

def get_wavelength_calibrated_fits_files(input_directory):
    """
        recursive search for fits files
        Give the directory containing multiple folders of observing dates and data. This definition will go through
        and look for the wavelength calibrated fits images and save their path locations to a list.
        output = list of paths to fits files
        """
    check = []
    fits_path = []
    for root, dirs, files in os.walk(input_directory, topdown=False):
        for name in files:
            if (os.path.join(root, name)) not in check and name.endswith('.fits'):
                if name.startswith('f.kic') or name.startswith('spec') or name.startswith('dkic') or name.startswith('f.'):
                        check.append(os.path.join(root, name))
                        fits_path.append(os.path.join(root, name))
    return fits_path


def save_ew(fitsimage, line, to_save):
    """output: text file with tab separated values of each kic number followed by EqW measurements"""

    try:
        fi = fitsimage.find('kic')
        tmp = fitsimage[fi + 3:]
        dot = tmp.find('.')
        kic = int(tmp[:dot])  # gets just the kic number
    except ValueError: #lamost data is different
        if 'spec' in str(fitsimage):
            kic = fitsimage[65:-5] # for lamost
        else:
            kic = fitsimage[57:-10]# for standard stars

    measurement = {}
    measurement[line] = to_save

    if kic not in ew_dict.keys():
        ew_dict[kic] = measurement
    else:
        stack = np.hstack((ew_dict[kic], measurement))
        ew_dict[kic] = stack


def ew_per_directory(path_list, print_values = False, plot_spectra = False, plot_per_line = False):
    """
    Runs the class EqWidth though Mdwarf_EqW.py

    output: a text file in the output_directory for each .fits file starting with dkic
            found in the parent directory

    plot_per_line = debugging tool for checking how well line values are measured
    """
    # change based on lowres or highres or calstar
    for image in path_list:
        if print_values:
            print image
        me = ME.EqWidth(image)

        for line in linelist:
            if len(line) == 1:
                try:
                    feature_width, peak = me.find_feature_region(line[0])
                    base1, base2, description, prewave, postwave = me.map_feature(feature_width, peak, line[0])
                    continuum, c1, c2 = me.find_continuum(base1, base2, 5)
                    ew, err = me.measure_ew(feature_width, peak, continuum, c1, c2, prewave, postwave, base1, base2, plot_region = False)
                    save_ew (image, line[0], [ew, err])
                    if print_values:
                        print 'ew', ew
                    if plot_per_line:
                        me.plot_ew(line,peak,  base1, base2, continuum, c1, c2)
                except Exception:
                    pass

            else:
                num = me.avg_flux_per_region(line[0], line[1])
                den = me.avg_flux_per_region(line[2], line[3])
                if num != None != den:
                    ratio = num/den
                    if print_values:
                        print 'ratio' + str(line[0]), ratio
                    save_ew(image, line[0], [ratio])

            try:
                if ME.EqWidth(image).subclass: # lamost data
                    save_ew(image, 'SpT', ME.EqWidth(image).subclass)
            except:
                pass


        if plot_spectra:
            plt.plot(me.x_wavelengths, me.y_data)
            plt.title(str(image[50:]))
            plt.show()

def write_to_outfile(outfile_name):
    """
    """

    assert type(outfile_name) == str
    # write output to file:
    df = pd.DataFrame(index = sorted(ew_dict.keys()), columns = outcolumns)
    for k, v in sorted(ew_dict.items()):
        for lines in v:
            try:
                for value, measurements in lines.items():
                    for n in xrange(len(outcolumns)):
                        try:
                            if value == linelist[n][0]:
                                df.loc[k, outcolumns[n]] = measurements
                        except IndexError:
                            if value == 'SpT':
                                df.loc[k, outcolumns[-1]] = measurements

            except AttributeError:
                pass

    # change as different out files are desired
    df.to_csv(outfile_name + '.csv', sep = '\t')


# The code in this function is executed when this file is run as a Python program
def main( Print = False, Plot = False):
    """

    """

    # file paths of various files
    #filepath = '/home/tessa/astronomy/mdwarf/highres_data'
    #filepath = '/home/tessa/astronomy/mdwarf/lowres_data'
    #filepath = '/home/tessa/astronomy/mdwarf/LAMOST'
    #filepath = '/home/tessa/astronomy/mdwarf/calibration_data'
    #filepath = '/home/tessa/astronomy/mdwarf/'
    current_directory = os.getcwd()

    path_list = get_wavelength_calibrated_fits_files(current_directory)
    ew_per_directory(path_list, print_values = Print, plot_spectra = Plot, plot_per_line = False)
    write_to_outfile('all_mdwarf_data')

if __name__ == "__main__":
    main()
