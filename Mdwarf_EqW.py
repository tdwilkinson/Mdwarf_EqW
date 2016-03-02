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
    lowerbound, upperbound = setrange(line, 5)

    # store the region near feature
    for i, j in zip(r, data):
        if i >= lowerbound and i <= upperbound:
            region[i] = j

    return region, lowerbound, upperbound

def find_peak_of_feature(r, data,line , continuum, plot):
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
    if np.abs(continuum - line_min[1]) > np.abs(continuum - line_max[1]):
        value = line_min
    else:
        value = line_max

    # look for turnover near continuum (# TODO: gaussian fit?)
    for n, flux in enumerate(region.values()):
        if int(flux) == int(value[1]):
            index = n

    return value

    # n = len(region_data)
    # print n
    # z = [region, region_data]
    # mean = np.sum(region * region_data)/n
    # m = np.mean(region_data)
    # print mean, m
    # sigma = np.sqrt(np.sum((region * region_data-mean)**2))
    # print sigma
    # popt,pcov = curve_fit(gaus,region_data,region,p0=[1,mean,sigma])
    # print popt

def gaus(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

def setrange(line, amount):
    # for measuring continuum
    lower, upper = [line - amount, line + amount]
    return lower, upper

def find_continuum(r, data, line):
    """
    Define region just outside of halpha to be continuum
    Currently set to 10 angstroms wide
    :return:
    """
    # TODO: set ha check to get max or min value in range
    # TODO: then set continuum to be much closer to actual ha dip/peak

    cont_range_preline = {}
    cont_range_postline = {}

    region, lowerbound, upperbound = find_feature_region(r, data, line)

    # define continuum outside of feature region
    for i, j in zip(r, data):

        if i >= lowerbound - 10 and i <= lowerbound:
            cont_range_preline[i] = j
        elif i >= upperbound and i <= upperbound + 10:
            cont_range_postline[i] = j

    pre, post = remove_cosmic_ray(cont_range_preline, cont_range_postline,1.25, 0.75)

    # define continuum
    c1 = np.mean(pre.values())
    c2 = np.mean(post.values())
    continuum = np.mean([c1, c2])
    contupper = continuum * 1.25  # limits for graphing
    contlower = continuum * .75  # limits for graphing

    return continuum, contupper, contlower, c1, c2, post

def remove_cosmic_ray(dict_lines_region1, dict_lines_region2, upperbound, lowerbound):
    """ bounds to test what is good. outputs two dictionaries
    """
    #
    # upperbound = 1.25
    # lowerbound = 0.75

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


def plot_ew(halpha, contupper, contlower, continuum, area, r, data, value):

    # plot halpha
    plt.axvline(x= halpha, color='orange')
    # plot ha feature based on min/max value
    plt.axvline(x= value[0], color = 'b')

    # plot boundries of halpha measurements
    for i in setrange(halpha, 5):
        plt.axvline(x=i, color='r', linestyle='--')

    # plot continuum
    plt.axhline(contupper, color='orange')  # go to get_ew2.py
    plt.axhline(contlower, color='orange')
    plt.axhline(continuum, color='red')
    lower, upper = setrange(halpha, 5)
    for i, j in zip(r, data):
        if i >= lower-10 and i <= lower:
            plt.axvline(x=i, linestyle='--', color='pink')

        elif i >= upper and i <= upper+10:
            plt.axvline(x=i, linestyle='--', color='pink')

    #plt.fill_between(halpha, continuum)

    # for fitting a gaussian
    #a = setrange(halpha)
    #plt.plot(a ,gaus(a,*popt),'b*:',label='fit')

    #plt.plot(area, color = 'red' )
    plt.plot(r, data, color='black', linestyle='-', marker=',')
    plt.xlim(6530, 6600)

    plt.show()



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
    continuum, contupper, contlower, c1, c2, ind2_nocr = find_continuum(r, data,halpha)

    # define region to investigate
    value = find_peak_of_feature(r, data, halpha, continuum, plot)


    # take the area under the curve
    ew, area = measure_ew(r, data, halpha, continuum, c1, c2, ind2_nocr)

    save_ew(fitsimage, ew)

    if plot:
        plot_ew(halpha, contupper, contlower, continuum, area, r, data, value)


def ew_per_directory(parent_directory):
    """
    output: a text file in the output_directory for each .fits file starting with dkic
            found in the parent directory
    """

    path_list = get_wavelength_calibrated_fits_files(parent_directory)
    for path in path_list:
        ew_per_image(path)




def roll(array, roll_amount):

    parameters = [array]
    ar = np.roll(array,1)

    return ar


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

