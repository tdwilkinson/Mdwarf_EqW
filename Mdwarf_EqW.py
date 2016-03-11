__author__ = 'Tessa Wilkinson'
"""
- python run_mdew.py in a directory with images or a directory of directories with images
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


class EqWidth:
    """ The whole sha-bang that finds feature equivalent widths """

    def __init__(self, image_path):

        # define data:
        image, header = pyfits.getdata(image_path, header=True)
        self.y_data = image[0][0, :]
        minw = header['CRVAL1']
        dw = header['CD1_1']

        # define wavelengths
        flux = np.arange(0, len(self.y_data))
        self.x_wavelengths = flux * dw + minw

        self.base1 = []
        self.base2 = []
        self.prewave = []
        self.postwave = []

    def find_feature_region(self, line):
        """
        """

        initial_region = {}
        feature_width = {}
        # this defines how many pixels out from the line to look
        lower, upper = setrange(line, 7)
        # store the region near feature
        for i, j in zip(self.x_wavelengths, self.y_data):
            if lower <= i <= upper:
                initial_region[i] = j

        if len(initial_region) == 0:
            raise Exception('poor flux calibration!')


        line_min = []
        line_max = []

        # find feature pixels near feature
        for i, j in zip(self.x_wavelengths, self.y_data):
            if j == min(initial_region.values()):
                line_min = [i, j]
            elif j == max(initial_region.values()):
                line_max = [i, j]

        # choose max or min feature based on depth
        if np.abs(np.mean(initial_region.values()) - line_min[1]) > np.abs(
                np.mean(initial_region.values()) - line_max[1]):
            feature_peak = line_min
        else:
            feature_peak = line_max

        # this defines how many pixels out from the feature to look
        lower, upper = setrange(feature_peak[0], 7)
        for i, j in zip(self.x_wavelengths, self.y_data):
            if lower <= i <= upper:
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

        self.prewave = prewave
        self.postwave = postwave  #TODO: not quite working yet! need to fix declaring selfs

        discription, start_point, end_point = base_point(prewave, postwave, feature_peak)

        self.base1 = start_point
        self.base2 = end_point

        return start_point, end_point, discription, prewave, postwave

    def find_continuum(self, base1, base2, width):
        """
        Define region just outside of halpha to be continuum
        base1 and base2 are 1D array of wavelength,flux values
        :return:
        """
        # save the data in before and after pockets
        cont_range_preline = {}
        cont_range_postline = {}

        for i, j in zip(self.x_wavelengths, self.y_data):
            if (int(base1[0]) - width) <= i <= int(base1[0]):
                cont_range_preline[i] = j
            elif int(base2[0]) <= i <= (int(base2[0]) + width):
                cont_range_postline[i] = j

        # pre, post = remove_cosmic_ray(cont_range_preline, cont_range_postline)

        # define continuum
        c1 = np.mean(cont_range_preline.values())
        c2 = np.mean(cont_range_postline.values())
        continuum = np.mean([c1, c2])

        return continuum, c1, c2

    def measure_ew(self, feature_width, peak, continuum, c1, c2, prewave, postwave, base1, base2, plot_region = True):
        '''measure area under feature'''

        #attempt 1:
        # take all data points in ha region and subtract continuum
        haflux = {}
        for wavelengths, fluxx in feature_width.iteritems():
            feature = fluxx - continuum
            haflux[wavelengths] = feature

        # sum all the feature data points
        area = trapz(haflux.values(), haflux.keys())
        fluxerr = np.std([c1, c2])  # std in continuum noise

        x = np.ones(len(haflux.values()))
        linerr = np.sqrt(trapz(haflux.values(), (fluxerr * x) ** 2))
        err = np.sqrt((linerr / (continuum)) ** 2 + (area / (continuum) ** 2 * np.std([c1, c2])) ** 2)
        ew = area / continuum

        # attempt 2
        x = []
        y = []
        for i in prewave:
            #if i[0] < base1[0]:
                x.append(i[0])
                y.append(i[1])
        x.append(peak[0])
        y.append(peak[1])
        for i in postwave:
            #if i[0] > base2[0]:
                x.append(i[0])
                y.append(i[1])
        n = np.ones(len(x))
        narea = trapz(y-continuum*n, x)
        newew = narea / continuum


        # attempt 3
        # def PolyArea(x,y):
        #     return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
        # print PolyArea(continuum, zip(x,y))


        #print y-continuum*n, haflux.values()

        if plot_region:

            #print continuum
            #print y - (continuum*n)
            plt.plot(x, y, linestyle = '-', marker = ',', color = 'b')
            plt.plot(x, continuum*n, linestyle = '-', color = 'r')
            plt.show()

        return newew, err

    def plot_ew(self, line, peak, base1, base2, continuum, c1, c2):
        """ plot to see what's being calculated   """

        # plot halpha
        plt.axvline(x=line, color='pink')

        # plot ha feature based on min/max value
        plt.axvline(x=peak[0], color='b')

        # plot boundries of halpha measurements
        plt.axvline(x=base1[0], color='g')
        plt.axvline(x=base2[0], color='g')
        plt.axvline(x=base1[0] - 5, linestyle='--', color='pink')
        plt.axvline(x=base2[0] + 5, linestyle='--', color='pink')

        # plot continuum
        plt.axhline(c1, color='orange')  # go to get_ew2.py
        plt.axhline(c2, color='orange')
        plt.axhline(continuum, color='red')

        plt.plot(self.x_wavelengths, self.y_data, color='black', linestyle='-', marker=',')
        plt.xlim(6530, 6600)
        plt.show()


    def avg_flux_per_region(self, start, end):
        """
        """
        # find feature pixels near feature
        region = {}
        for i, j in zip(self.x_wavelengths, self.y_data):
                if start <= i <= end:
                    region[i] = j

        if len(region.values()) > 0:
            average = sum(region.values()) / len(region.values())
            return average
########################################################################################################################
########################################################################################################################


def base_point(prewave, postwave, featurepeak):
    """ Looks at the two halfs of a line feature and finds the turn point at the bases    """

    # find the slopes by comparing the first point in the wave to the peak in the feature
    preslope = slope(prewave[0], featurepeak)
    postslope = slope(featurepeak, postwave[-1])

    dip = True
    str_feature = 'not defined'

    if preslope < 0 < postslope:
        str_feature = 'absorption'
        dip = True
    elif preslope > 0 > postslope:
        str_feature = 'emission!'
        dip = False
    else:
        pass
        print 'borked'

    # now check where the turning point is : where its no longer + or -
    greenline1 = prewave[0]
    greenline2 = postwave[-1]

    if dip:
        for i, e in list(enumerate(postwave)):
            if i < len(postwave) - 2 and postwave[i][1] > postwave[i + 1][1]:
                greenline2 = postwave[i]
                break
        for i, e in reversed(list(enumerate(prewave))):
            if i > 0 and prewave[i][1] > prewave[i - 1][1]:
                greenline1 = prewave[i]
                break
    else:
        for i, e in list(enumerate(postwave)):
            if i < len(postwave) - 2 and postwave[i][1] < postwave[i + 1][1]:
                greenline2 = postwave[i]
                break
        for i, e in reversed(list(enumerate(prewave))):
            if i > 0 and prewave[i][1] < prewave[i - 1][1]:
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



def slope(list1, list2):
    """  change in y / change in x , comparing start to mid and mid to end points"""
    w1, f1 = list1
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
