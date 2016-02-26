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

def find_halpha(r, data, plot):
    """
     define halpha
    :param r:
    :param data: the flux values for the data
    :param plot: if you want to plot these values
    :return:variables at the moment
    """

    ha_region = {}
    lower, upper = setrange(halpha)

    for i, j in zip(r, data):
        if i >= lower and i <= upper:
            ha_region[i] = j

    return ha_region

def find_continuum(r, data, flux):
    """
    Define region just outside of halpha to be continuum
    Currently set to 10 angstroms wide
    :return:
    """
    # TODO: set ha check to get max or min value in range
    # TODO: then set continuum to be much closer to actual ha dip/peak

    cont_range_preline = {}
    cont_range_postline = {}
    preline_nocr = {}
    postline_nocr = {}


    lower, upper = setrange(halpha)
    for i, j in zip(r, data):
        if i >= lower-10 and i <= lower:
            cont_range_preline[i] = j
        elif i >= upper and i <= upper+10:
            cont_range_postline[i] = j

    # remove cosmic rays from the continuum regions by setting limits
    for n, (w, flux) in enumerate(cont_range_preline.iteritems()):
        if cont_range_preline.values()[n - 1] < flux * 1.25 and cont_range_postline.values()[n - 1] > flux * 0.75:
            preline_nocr[w] = flux
    for n, (w, flux) in enumerate(cont_range_postline.iteritems()):
        if cont_range_postline.values()[n - 1] < flux * 1.25 and cont_range_postline.values()[n - 1] > flux * 0.75:
            postline_nocr[w] = flux

    # define continuum
    c1 = np.mean(preline_nocr.values())
    c2 = np.mean(postline_nocr.values())
    continuum = np.mean([c1, c2])
    contupper = continuum * 1.25  # limits for graphing
    contlower = continuum * .75  # limits for graphing

    return continuum, contupper, contlower, c1, c2, postline_nocr

def measure_ew(ha_region, continuum, c1, c2, ind2_nocr):
    # measure area under feature
    haflux = {}
    for wavelengths, fluxx in ha_region.iteritems():  # take all data points in ha region and subtract continuum
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


def plot_ew(halpha, contupper, contlower, continuum, area, r, data,):

    # plot halpha
    plt.axvline(x=halpha, color='g')
    # plot boundries of halpha measurements
    for i in setrange(halpha):
        plt.axvline(x=i, color='g', linestyle='--')
    # plot continuum
    plt.axhline(contupper, color='orange')  # go to get_ew2.py
    plt.axhline(contlower, color='orange')
    plt.axhline(continuum, color='red')
    lower, upper = setrange(halpha)
    for i, j in zip(r, data):
        if i >= lower-10 and i <= lower:
            plt.axvline(x=i, linestyle='--', color='pink')

        elif i >= upper and i <= upper+10:
            plt.axvline(x=i, linestyle='--', color='pink')

    #plt.fill_between(halpha, continuum)
    #plt.plot(area, color = 'red' )
    plt.plot(r, data, color='black', linestyle='-', marker=',')
    plt.xlim(6530, 6600)

    plt.show()



def ew_per_image(fitsimage, checkmode = False, plot = False):
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

    # define region to investigate
    ha_region= find_halpha(r, data, plot)
    # clean investigation region of cosmic rays and define continuum region
    continuum, contupper, contlower, c1, c2, ind2_nocr = find_continuum(r, data, flux)
    # take the area under the curve
    ew, area = measure_ew(ha_region, continuum, c1, c2, ind2_nocr)

    save_ew(fitsimage, ew)

    if plot:
        plot_ew(halpha, contupper, contlower, continuum, area, r, data)


def ew_per_directory(parent_directory):
    """
    output: a text file in the output_directory for each .fits file starting with dkic
            found in the parent directory
    """

    path_list = get_wavelength_calibrated_fits_files(parent_directory)
    for path in path_list:
        ew_per_image(path)


def setrange(line):
    # for measuring continuum
    lower, upper = [line - 5, line + 5]
    return lower, upper

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

