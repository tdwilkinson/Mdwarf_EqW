__author__ = 'tdwilkinson'


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


class mdew_vs_spt(object):

    def __init__(self, path_to_ew, path_to_spt):
        """
        input: path_to_ew: a text file, the first column is the kic numbers, the second column is a dictionary with
            each line measured as the keys and the value is the measurement and error of the line
        """

        # load in data
        ew_df = pd.DataFrame.from_csv(path_to_ew, sep='\t')
        ha = ew_df['ha']
        #ha_err = ew_df['ha_ew_err']
        tio1 = ew_df['tio1']
        #tio1_err = ew_df['tio1_err']
        tio2 = ew_df['tio2']
        #tio2_err = ew_df['tio2_err']
        tio3 = ew_df['tio3']
        #tio3_err = ew_df['tio3_err']
        tio4 = ew_df['tio4']
        #tio4_err = ew_df['tio4_err']
        tio5 = ew_df['tio5']
        #tio5_err = ew_df['tio5_err']
        cah2 = 	ew_df['cah2']
        #cah2_err = ew_df['cah2_err']
        caoh = ew_df['caoh']
        #caoh_err = ew_df['caoh_err']
        o2 = ew_df['o2']
        #o2_err = ew_df['o2_err']
        n0 = ew_df['n0']
        #no_err = ew_df['no_err']

        # self.ewfile_contents = [ha, ha_err, tio1, tio1_err, tio2, tio2_err, tio3, tio3_err, tio4, tio4_err, tio5, tio5_err, cah2,\
        #         cah2_err, caoh, caoh_err, o2, o2_err, n0, no_err]
        self.linelist = [ha, tio1, tio2, tio3, tio4, tio5, cah2, caoh, o2, n0]


        spt_df = pd.DataFrame.from_csv(path_to_spt,sep= '\t')

        spt = spt_df['spectype']
        # TODO: add in error column too!

        self.spectypes = ['K5', 'K7', 'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', '?']
        #store data
        sptdata  = {}
        ewdata = {}

        for kic, sptype in spt.iteritems():
            if str(sptype) != 'bad':
                sptdata[int(kic)] = [sptype]
        self.sptype = sptdata
        self.xvalue = {}


    def spt_index(self):
        """ outputs dictionary of spectra type index's
        :rtype : dictionary
        """

        xvalue = {}
        for kic, value in self.sptype.items():
            for st in self.spectypes:
                if st == value[0]:
                    xpos = self.spectypes.index(st)
                    xvalue[int(kic)] = xpos
        self.xvalue = xvalue

        return xvalue

    def line_values(self, line, line_err):

        yvalue = {}
        for kicnumber, measurement in line:
            for kic, error in line_err:
                if kicnumber == kic:
                    yvalue[kic] = [line, line_err]


    def plot(self):

            linelist = ['ha', 'tio1', 'tio2', 'tio3', 'tio4', 'tio5', 'cah2', 'caoh', 'o2', 'n0']

            def get_line(something, index):
                return self.linelist[index].get(something) # .get gets a string

            spt = self.spt_index()
            for i in xrange(len(self.linelist)):


                for kic, index in spt.items():
                    for k,v in self.linelist[i].iteritems():
                        if k == kic:
                            try:
                                msmt = v[1:v.find(',')]
                            except AttributeError:
                                msmt = v
                            if msmt > 0 :
                                plt.scatter(index, msmt, marker = '*', alpha = 0.7)
                            #plt.errorbar(index, msmt, yerr= err, marker = '*', alpha = 0.5)
                        # if get_line(kic, i) != np.nan and get_line(kic,i) != None:
                        # try:
                        #     print get_line(kic,i), type(get_line(kic,i))
                        #     msmt = get_line(kic,i).split()
                        #     print msmt[0]
                        #     err = get_line(kic,i)
                        #     plt.errorbar(index, msmt, yerr= err, marker = '*', alpha = 0.5)
                        # except:
                        #     print 'nope'
                        #     plt.scatter(index, get_line(kic, i), marker= '*', alpha = 0.5, color = 'b')

                plt.xticks(np.arange(len(self.spectypes)), self.spectypes)
                plt.ylabel(linelist[i])
                plt.xlabel('Spectral Types')
                plt.title('Magnetic Activity per Spectral Type')
                plt.show()

#############################################################################################
#############################################################################################

mdew_vs_spt('finalEW.txt', '/home/tessa/PycharmProjects/mdwarf/finalSpT.txt').plot()