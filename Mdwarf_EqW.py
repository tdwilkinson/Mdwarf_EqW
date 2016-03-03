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
from scipy import asarray as ar,exp
import itertools
from numpy.polynomial.polynomial import polyfit, polyval

# global variables = which lines to measure
halpha = 6562.801
tio5 = 0

# set up output file and header
newfile = open('finalEW.txt', 'w')
newfile.write('kicnumber' + '\t' + 'ha_ew' + '\t' + 'ha_ew_err' + '\n')


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

def find_feature_region(r, data, line):

    region = {}
    lowerbound, upperbound = setrange(line, 7)
    # store the region near feature
    for i, j in zip(r, data):
        if i >= lowerbound and i <= upperbound:
            region[i] = j
    return region , lowerbound, upperbound

def find_peak_of_feature(r, data, line, plot):
    """
     define halpha
     halpha = 6563 as defined above. value is to find the min/max point of features nearest the line
    :param data: the flux values for the data
    :param plot: if you want to plot these values
    :return:variables at the moment
    """
    region, lowerbound, upperbound = find_feature_region(r, data, line)
    # find feature pixels near feature
    for i, j in zip(r, data):
        if j == min(region.values()):
            line_min = [i,j]
        elif j == max(region.values()):
            line_max = [i, j]
    # choose max or min feature based on depth
    if np.abs(np.mean(region.values()) - line_min[1]) > np.abs(np.mean(region.values()) - line_max[1]):
        value = line_min
    else:
        value = line_max
    return value

def map_feature(r, data,line, plot):
    """
    maps out the line feature to find where the curve turns over
    returns: start and end pixel arrays of feature

    """
    # TODO: work on this. refactor.

    value = find_peak_of_feature(r, data, line, plot)
    region, lowerbound, upperbound = find_feature_region(r, data, value[0])

    wave = []
    flux = []
    weight = []
    prewave = []
    postwave = []
    option1 = []
    option2 = []
    option3 = []
    option4 = []

    n = len(region.keys())
    region_keys = sorted(region)
    for i in region_keys:
        wave.append(i)
        flux.append(region[i])

    # add weights to feature peak
    for i in xrange(len(wave)):
        if int(wave[i]) == int(value[0]):
            w = 10
        else:
            w = 1
        weight.append(w)

    # calculate polynomial
    fit = polyfit(wave, flux, 2, w = weight)
    val = polyval(wave, fit)

    # mark areas before the wave and after the wave
    region_keys = sorted(region)
    for i in region_keys:
        if i <= value[0]:
            prewave.append([i, region[i]])
        elif i >= value[0]:
            postwave.append([i, region[i]])

    # sorry repetition but follow each side of wave until turn over
    for i in xrange(len(prewave)-1, -1, -1):
        if i >= 2:
            first = prewave[i][1]
            second = prewave[i-1][1]
            third = prewave[i-2][1]
            if first > second and second > third:
                option1.append(prewave[i])
            elif first < second and second < third:
                option2.append(prewave[i])

    for i in xrange(0,len(postwave), 1):
        if i < len(postwave)-2:
            first = postwave[i][1]
            second = postwave[i+1][1]
            third = postwave[i+2][1]
            if first > second and second > third:
                option3.append(postwave[i])
            elif first < second and second < third:
                option4.append(postwave[i])

    print option1
    print option2
    print option3
    print option4


    # determine where turnover is despite positive or negative slope of feature
    if len(option1) < len(option2):
        print len(option1), len(option2)
        for i in option1:
            print '1', i
            start = option1[0][0]
    else:
        print len(option1), len(option2)
        for i in option2:
            print '2', i
            start = option2[0][0]

    if len(option3) < len(option4):
        end = option3[0][0]

    else:
        end = option4[0][0]
    return  start, end, wave,val


def find_continuum(r, data, line, plot):
    """
    Define region just outside of halpha to be continuum
    Currently set to 10 angstroms wide
    :return:
    """

    cont_range_preline = {}
    cont_range_postline = {}

    # get feature region and the edges of it based on map feature

    start_feature, end_feature, wave, val = map_feature(r, data,line, plot)
    region, lowerbound, upperbound = find_feature_region(r, data, start_feature[0])

    #print lowerbound, start_feature[0], end_feature[0],  upperbound

    # define continuum region
    for i, j in zip(r, data):
        if i >= lowerbound and i <= start_feature[0]:
            cont_range_preline[i] = j
        elif i >= end_feature[0] and i <= upperbound:
            cont_range_postline[i] = j

    pre, post = remove_cosmic_ray(cont_range_preline, cont_range_postline)

    # define continuum
    c1 = np.mean(pre.values())
    c2 = np.mean(post.values())
    continuum = np.mean([c1, c2])

    return continuum, c1, c2, post, start_feature, end_feature, wave, val

def remove_cosmic_ray(dict_lines_region1, dict_lines_region2):
    """ bounds to test what is good. outputs two dictionaries
    """
    upperbound = 1.25
    lowerbound = 0.75

    preline_nocr = {}
    postline_nocr = {}

    # remove cosmic rays from the continuum regions by setting limits
    for n, (w, flux) in enumerate(dict_lines_region1.iteritems()):
        if dict_lines_region1.values()[n - 1] < flux * upperbound and dict_lines_region1.values()[n - 1] > flux * lowerbound:
            preline_nocr[w] = flux
    for n, (w, flux) in enumerate(dict_lines_region2.iteritems()):
        if dict_lines_region2.values()[n - 1] < flux * upperbound and dict_lines_region2.values()[n - 1] > flux * lowerbound:
            postline_nocr[w] = flux
    return preline_nocr, postline_nocr

def measure_ew(r, data, line, continuum, c1, c2, ind2_nocr):
        # measure area under feature

        haflux = {}
        region, lowerbound, upperbound = find_feature_region(r, data, line)

        for wavelengths, fluxx in region.iteritems():  # take all data points in ha region and subtract continuum
            feature = fluxx - continuum
            haflux[wavelengths] = feature

        area = trapz(haflux.values(), haflux.keys())  # sum all the feature data points
        fluxerr = np.std(ind2_nocr.values())  # std in continuum noise

        x = np.ones(len(haflux.values()))
        linerr = np.sqrt(trapz(haflux.values(), (fluxerr*x)**2))
        ha_err = np.sqrt((linerr / (continuum)) ** 2 + (area / (continuum) ** 2 * np.std([c1,c2]))**2)
        ew = [area / continuum, ha_err]

        return ew, area

def save_ew(fitsimage, ew):
    """
    output: text file with tab separated values of each kic number followed by EqW measurements
    """

    ew_dict = {}

    fi = fitsimage.find('kic')
    tmp = fitsimage[fi + 3:]
    dot = tmp.find('.')
    kic = int(tmp[:dot])  # gets just the kic number

    if kic not in ew_dict.keys():
        ew_dict[kic] = ew

    # TODO: this vstack works, just need to apply weighted-average? to print best to finalEW.txt
    # else:
    #     stack = np.vstack((ew))
    #     ew_dict[kic] = stack


    for k, v in sorted(ew_dict.items()):
        newfile.write(str(k) + '\t' + str(v[0]) + '\t' + str(v[1]) + '\n')


def plot_ew(halpha, continuum,c1, c2, area, r, data, value, start, end, wave, val):

    # plot halpha
    plt.axvline(x= halpha, color='orange')
    # plot ha feature based on min/max value
    plt.axvline(x= value[0], color = 'b')

    # plot boundries of halpha measurements
    plt.axvline(x= start[0], color = 'g')
    plt.axvline(x= end[0], color = 'g')

    # plot continuum

    plt.axhline(c1, color='orange')  # go to get_ew2.py
    plt.axhline(c2, color='orange')
    plt.axhline(continuum, color='red')
    lower, upper = setrange(halpha, 5)

    # plt.axvline(x=i, linestyle='--', color='pink')
    # plt.axvline(x=i, linestyle='--', color='pink')

    #plt.fill_between(halpha, continuum)

    # for fitting a gaussian
    #a = setrange(halpha)
    #plt.plot(a ,gaus(a,*popt),'b*:',label='fit')
    # plt.plot(wave,flux,'o', x_new, y_new)

    ## This one works but is not well fitted! plt.plot(wave, val, 'm')
    # plt.xlim([wave[0]-1, x[-1] + 1 ])


    plt.plot(area, color = 'red' )
    plt.plot(r, data, color='black', linestyle='-', marker=',')
    plt.xlim(6530, 6600)

    plt.show()

def roll(array, roll_amount):

    parameters = [array]
    ar = np.roll(array,1)
    return ar

def gaus(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

def func(x, a, b, c):
     return a * np.exp(-b * x) + c

def setrange(line, amount):
    # for measuring continuum
    lower, upper = [line - amount, line + amount]
    return lower, upper

def ew_per_directory(parent_directory):
    """
    output: a text file in the output_directory for each .fits file starting with dkic
            found in the parent directory
    """

    path_list = get_wavelength_calibrated_fits_files(parent_directory)
    for path in path_list:
        ew_per_image(path)


def ew_per_image(fitsimage, checkmode = False, plot = True):
    """
     :return: output file with ew's
    """
    # get data:
    image, header = pyfits.getdata(fitsimage, header=True)
    data = image[0][0, :]
    minw = header['CRVAL1']
    dw = header['CD1_1']

    #define wavelengths
    flux = np.arange(0, len(data))
    r = flux * dw + minw

    # clean investigation region of cosmic rays and define continuum region
    continuum, c1, c2, ind2_nocr , start, end, wave,val= find_continuum(r, data,halpha, plot)

    # define region to investigate
    value = find_peak_of_feature(r, data, halpha, plot)

    # take the area under the curve
    ew, area = measure_ew(r, data, halpha, continuum, c1, c2, ind2_nocr)

    save_ew(fitsimage, ew)

    if plot:
        plot_ew(halpha, continuum, c1, c2, area, r, data, value, start, end, wave,val)


def final_plot():
    """

    :return: plot output values of ew per star
    """

    kicnumber, ha, ha_err = np.loadtxt('finalEW.txt', dtype = (int, float, float), skiprows = 1, unpack = True )

    print kicnumber, ha ,ha_err
########################################################################################################################
########################################################################################################################


# TODO: make sure multiple files for one star or stored or appended in list (see todo above)
# TODO: separate ones with magnetic activity
# TODO: add in looking at other bands
# TODO: plot output values (maybe in new file)

# 1) create output file finalEW.txt

# change this line to search a directory that contains directories or files.
ew_per_directory('/home/tessa/astronomy/mdwarf/highres_data')
newfile.close()

