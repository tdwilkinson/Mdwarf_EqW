__author__ = 'Tessa Wilkinson'

# copy of get_ew3.py
# run in python console to get
# outfile: finalEW.txt

import os
import numpy as np
import pyfits
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import polyfit, polyval
from scipy.interpolate import UnivariateSpline as spline
import pandas as pd
from numpy import trapz  # trapezoidal sum of the area


def get_wavelength_calibrated_fits_files(input_directory):

    # fits_files = []
    # for file in os.listdir(input_directory):
    #     if file.startswith('dkic') and file.endswith('.fits'):
    #         fits_files.append(file)
    # return fits_files

    check = []
    fits_path = []
    for root, dirs, files in os.walk(highpath, topdown=False):
        for name in files:
            if (os.path.join(root, name)) not in check and name.startswith('dkic') and name.endswith('.fits'):
                # print(os.path.join(root,name)) # to get whole path to fits file
                check.append(os.path.join(root, name))

                fits_path.append(os.path.join(root, name))

    return fits_path


def find_halpha(r, data, plot):
    # define halpha
    halpha = 6562.801
    a = [halpha - 5, halpha + 5]
    ind1 = {}
    ind2 = {}
    ha_region = {}
    for i, j in zip(r, data):
        if i >= a[0]-10 and i <= a[0]:
            ind1[i] = j
            if plot:
                plt.axvline(x=i, linestyle='--', color='pink')

        elif i >= a[1] and i <= a[1]+10:
            ind2[i] = j
            if plot:
                plt.axvline(x=i, linestyle='--', color='pink')

        elif i >= a[0] and i <= a[1]:
            ha_region[i] = j

    return ha_region, ind1, ind2, halpha, a

def find_continuum(ind1, ind2, flux):
    # remove cosmic rays from the continuum regions by setting
    ind1_nocr = {}
    ind2_nocr = {}
    for n, (w, flux) in enumerate(ind1.iteritems()):
        if ind1.values()[n - 1] < flux * 1.25 and ind1.values()[n - 1] > flux * 0.75:
            ind1_nocr[w] = flux
    for n, (w, flux) in enumerate(ind2.iteritems()):
        if ind2.values()[n - 1] < flux * 1.25 and ind2.values()[n - 1] > flux * 0.75:
            ind2_nocr[w] = flux

        # define continuum
    c1 = np.mean(ind1_nocr.values())
    c2 = np.mean(ind2_nocr.values())
    continuum = np.mean([c1, c2])
    contupper = continuum * 1.25  # limits for graphing
    contlower = continuum * .75  # limits for graphing

    return continuum, contupper, contlower, c1, c2, ind2_nocr

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
    ew_dict = {}
    # single out kic to compare with prot kic
    fi = fitsimage.find('kic')
    tmp = fitsimage[fi + 3:]
    dot = tmp.find('.')
    kic = int(tmp[:dot])  # gets just the kic number

    if kic not in ew_dict.keys():
        ew_dict[kic] = []
    ew_dict[kic].append(ew)

    newfile = open('finalEW.txt', 'w')
    newfile.write('kicnumber' + '\t' + 'ha_ew' + '\t' + 'ha_ew_err' + '\n')
    for k, v in sorted(ew_dict.iteritems()):
        newfile.write(str(k) + '\t' + str(v[0][0]) + '\t' + str(v[0][1]) + '\n')
    newfile.close()

def plot_ew(halpha, a, contupper, contlower, continuum, area, r, data):

    figure = plt.figure()
    # plot halpha
    plt.axvline(x=halpha, color='g')
    # plot boundries of halpha measurements
    for i in a:
        plt.axvline(x=i, color='g', linestyle='--')
    # plot continuum
    plt.axhline(contupper, color='orange')  # go to get_ew2.py
    plt.axhline(contlower, color='orange')
    plt.axhline(continuum, color='red')
    #plt.fill_between(halpha, continuum)
    #plt.plot(area, color = 'red' )
    plt.plot(r, data, color='black', linestyle='-', marker=',')
    plt.xlim(6530, 6600)

    figure.show()

# TODO: currently pops up all ~80, may want to fix that before turning plot on! :)

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
    ha_region, ind1, ind2, halpha, a = find_halpha(r, data, plot)

    # clean investigation region of cosmic rays and define continuum region
    continuum, contupper, contlower, c1, c2, ind2_nocr = find_continuum(ind1,ind2, flux)

    ew, area = measure_ew(ha_region, continuum, c1, c2, ind2_nocr)

    if plot:
        plot_ew(halpha, a, contupper, contlower, continuum, area, r, data)

    save_ew(fitsimage, ew)

    # print ew_dict
    # if checkmode:
    #     print 'kicnumber:', kic
    #     print 'continuum mean:', continuum
    #     print 'area:', area
    #     print 'fluxerr:', fluxerr
    #     print 'ew:', ew


def ew_per_directory(parent_directory, output_directory):
    """
    output: a text file in the output_directory for each .fits file starting with dkic
            found in the parent directory
    """

    path_list = get_wavelength_calibrated_fits_files(parent_directory)
    for path in path_list:
        ew_per_image(path)


# TODO: plot only one at a time
# TODO: correct values printing to finalEW.txt
# TODO: make sure multiple files for one star or stored or something
# TODO: separate ones with magnetic activity
# TODO: add in looking at other bands



directory = (os.listdir('/home/tessa/astronomy/mdwarf/'))
highpath = '/home/tessa/astronomy/mdwarf/highres_data'
lowpath = '/home/tessa/astronomy/mdwarf/lowres_data'

ew_per_directory(highpath, os.path.join(highpath, 'disp_spectra_data'))
