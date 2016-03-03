__author__ = 'Tessa Wilkinson'
"""
- run in python console
- change the path at the bottom of the line to search a directory of UT directories for images.
 outfile: finalEW.txt
"""
import os
import numpy as np
import pyfits
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import polyfit, polyval
from scipy.interpolate import UnivariateSpline as spline
import pandas as pd
from numpy import trapz  # trapezoidal sum of the area
from scipy.optimize import curve_fit
from scipy import asarray as ar, exp
import itertools
from numpy.polynomial.polynomial import polyfit, polyval

# global variables = which lines to measure
halpha = 6562.801
tio5 = 0

# set up output file and header
ew_dict = {}
newfile = open('finalEW.txt', 'w')
newfile.write('kicnumber' + '\t' + 'ha_ew' + '\t' + 'ha_ew_err' + '\n')


class EqWidth():
    """
    """
    def __init__(self, image_path):

        # define data:
        image, header = pyfits.getdata(image_path, header=True)
        self.y_data = image[0][0, :]
        minw = header['CRVAL1']
        dw = header['CD1_1']

        # define wavelengths
        flux = np.arange(0, len(self.y_data))
        self.x_wavelengths = flux * dw + minw



    def find_feature_region(self, line):
        """
        """

        initial_region = {}
        feature_width = {}
        # this defines how many pixels out from the line to look
        lower, upper = setrange(line, 7)
        # store the region near feature
        for i, j in zip(self.x_wavelengths, self.y_data):
            if i >= lower and i <= upper:
                initial_region[i] = j

        # find feature pixels near feature
        for i, j in zip(self.x_wavelengths, self.y_data):
            if j == min(initial_region.values()):
                line_min = [i, j]
            elif j == max(initial_region.values()):
                line_max = [i, j]

        # choose max or min feature based on depth
        if np.abs(np.mean(initial_region.values()) - line_min[1]) > np.abs(np.mean(initial_region.values()) - line_max[1]):
            feature_peak = line_min
        else:
            feature_peak = line_max

        # this defines how many pixels out from the feature to look
        lower, upper = setrange(feature_peak[0], 7)
        for i, j in zip(self.x_wavelengths, self.y_data):
            if i >= lower and i <= upper:
                feature_width[i] = j

        return feature_width, feature_peak

    def map_feature(self, feature_width, feature_peak, line):
        """
        maps out the line feature to find where the curve turns over
        feature_width is a dictionary
        feature_peak is a 2d array
        returns: start and end pixel arrays of feature

        """
        # store waves and flux values ordered and separated
        wave = []
        flux = []
        region_keys = sorted(feature_width)
        for i in region_keys:
            wave.append(i)
            flux.append(feature_width[i])

        # store areas before the wave and after the wave
        prewave = []
        postwave = []
        region_keys = sorted(feature_width)
        for i in region_keys:
            if i <= feature_peak[0]:
                prewave.append([i, feature_width[i]])
            elif i >= feature_peak[0]:
                postwave.append([i, feature_width[i]])

        discription, start_point, end_point = base_point(prewave, postwave, feature_peak)

        self.base1 = start_point
        self.base2 = end_point

        return start_point, end_point, discription


    def find_continuum(self, base1, base2, amount):
        """
        Define region just outside of halpha to be continuum
        base1 and base2 are 1D array of wavelength,flux values
        :return:
        """
        # save the data in before and after pockets
        cont_range_preline = {}
        cont_range_postline = {}

        for i, j in zip(self.x_wavelengths, self.y_data):
            if i >= (int(base1[0]) - amount) and i <= int(base1[0]):
                cont_range_preline[i] = j
            elif i >= int(base2[0]) and i <= (int(base2[0]) + amount):
                cont_range_postline[i] = j

        # pre, post = remove_cosmic_ray(cont_range_preline, cont_range_postline)

        # define continuum
        c1 = np.mean(cont_range_preline.values())
        c2 = np.mean(cont_range_postline.values())
        continuum = np.mean([c1, c2])

        return continuum, c1, c2

    def measure_ew(self, feature_width, peak, base1, base2, continuum, c1, c2):
        '''measure area under feature'''

        # take all data points in ha region and subtract continuum
        haflux = {}
        for wavelengths, fluxx in feature_width.iteritems():
            feature = fluxx - continuum
            haflux[wavelengths] = feature

        # sum all the feature data points
        area = trapz(haflux.values(), haflux.keys())
        fluxerr =  np.std([c1,c2])  # std in continuum noise

        x = np.ones(len(haflux.values()))
        linerr = np.sqrt(trapz(haflux.values(), (fluxerr * x) ** 2))
        ha_err = np.sqrt((linerr / (continuum)) ** 2 + (area / (continuum) ** 2 * np.std([c1, c2])) ** 2)
        ew = [area / continuum, ha_err]

        return ew, area

    def plot_ew(self, line, peak, base1, base2, continuum, c1, c2, ew):
        ''' plot to see what's being calculated   '''

        # plot halpha
        plt.axvline(x= line, color='pink')

        # plot ha feature based on min/max value
        plt.axvline(x= peak[0], color='b' )

        # plot boundries of halpha measurements
        plt.axvline(x=base1[0], color='g')
        plt.axvline(x=base2[0], color='g')
        plt.axvline(x=base1[0] - 5, linestyle='--', color = 'pink')
        plt.axvline(x=base2[0] + 5, linestyle='--', color = 'pink')

        # plot continuum
        plt.axhline(c1, color='orange')  # go to get_ew2.py
        plt.axhline(c2, color='orange')
        plt.axhline(continuum, color='red')

        plt.plot(self.x_wavelengths, self.y_data, color='black', linestyle='-', marker=',')
        plt.xlim(6530, 6600)
        plt.show()

    def final_plot():
        """

        :return: plot output values of ew per star
        """

        kicnumber, ha, ha_err = np.loadtxt('finalEW.txt', dtype=(int, float, float), skiprows=1, unpack=True)

        #print kicnumber, ha, ha_err



########################################################################################################################
########################################################################################################################

def base_point(prewave, postwave, featurepeak):
    """ Looks at the two halfs of a line feature and finds the turn point at the bases    """

    # find the slopes by comparing the first point in the wave to the peak in the feature
    preslope = slope(prewave[0], featurepeak)
    postslope = slope(featurepeak, postwave[-1])
    if preslope < 0 and postslope > 0:
        str_feature = 'absorption'
        dip = True
    elif preslope > 0 and postslope < 0:
        str_feature = 'emission!'
        dip = False
    else:
        print 'borked'

    # now check where the turning point is : where its no longer + or -
    if dip:
        for i, e in list(enumerate(postwave)):
            if i < len(postwave)-2 and postwave[i][1] > postwave[i+1][1]:
                greenline2 = postwave[i]
                break
        for i, e in reversed(list(enumerate(prewave))):
             if i > 0 and prewave[i][1] > prewave[i-1][1]:
                greenline1 = prewave[i]
                break
    else:
        for i, e in list(enumerate(postwave)):
            if i < len(postwave)-2 and postwave[i][1] < postwave[i+1][1]:
                greenline2 = postwave[i]
                break
        for i, e in reversed(list(enumerate(prewave))):
             if i > 0 and prewave[i][1] < prewave[i-1][1]:
                greenline1 = prewave[i]
                break

    return str_feature, greenline1, greenline2

def remove_cosmic_ray(region1, region2):
    """ bounds here for testing. outputs two dictionaries"""

    upperbound = 1.25
    lowerbound = 0.75

    preline_nocr = {}
    postline_nocr = {}

    # remove cosmic rays from the continuum regions by setting limits
    for n, (w, flux) in enumerate(dict_lines_region1.iteritems()):
        if dict_lines_region1.values()[n - 1] < flux * upperbound and dict_lines_region1.values()[
                    n - 1] > flux * lowerbound:
            preline_nocr[w] = flux
    for n, (w, flux) in enumerate(dict_lines_region2.iteritems()):
        if dict_lines_region2.values()[n - 1] < flux * upperbound and dict_lines_region2.values()[
                    n - 1] > flux * lowerbound:
            postline_nocr[w] = flux

    return preline_nocr, postline_nocr

def roll(array, roll_amount):
    parameters = [array]
    ar = np.roll(array, 1)
    return ar

def slope(list1, list2):
    '''  change in y / change in x , comparing start to mid and mid to end points'''
    w1, f1 =  list1
    w2, f2 = list2
    return (f2 - f1) / (w2 - w1)

def gaus(x, a, x0, sigma):
    return a * exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def func(x, a, b, c):
    return a * np.exp(-b * x) + c


def setrange(line, amount):
    # for measuring continuum
    lower, upper = [line - amount, line + amount]
    return lower, upper


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

    if kic not in ew_dict.keys():
        ew_dict[kic] = ew

    for k, v in sorted(ew_dict.items()):
        newfile.write(str(k) + '\t' + str(v[0]) + '\t' + str(v[1]) + '\n')


def ew_per_directory(parent_directory, plot_per_image = True):
    """
    output: a text file in the output_directory for each .fits file starting with dkic
            found in the parent directory
    """

    linelist = [halpha]
    path_list = get_wavelength_calibrated_fits_files(parent_directory)
    for image in path_list:
        for line in linelist:
            feature_width, peak = EqWidth(image).find_feature_region(line)
            base1, base2, discription = EqWidth(image).map_feature(feature_width, peak, line)
            continuum, c1, c2 = EqWidth(image).find_continuum(base1, base2, 5)
            ew = EqWidth(image).measure_ew(feature_width, peak, base1, base2, continuum, c1, c2)

            print discription
            save_ew(image, line,ew)


        if plot_per_image:
            EqWidth(image).plot_ew(line, peak, base1, base2, continuum, c1, c2, ew)


# TODO: make sure multiple files for one star or stored or appended in list (see todo above)
# TODO: separate ones with magnetic activity
# TODO: add in looking at other bands
# TODO: plot output values (maybe in new file)


# change this line to search a directory that contains directories or files.
ew_per_directory('/home/tessa/astronomy/mdwarf/highres_data', plot_per_image=True)
newfile.close()
