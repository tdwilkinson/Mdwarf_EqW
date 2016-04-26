__author__ = 'tdwilkinson'

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def combine_values(dictionary1, dictionary2, *dictionary3):
    """ appends the value from dictionary2 to an array"""

    out = {}
    for key, value in dictionary1.iteritems():
        for key1, value1 in dictionary2.iteritems():
            key = str(key)
            key1 = str(key1)
            if type(value) == str:
                value = float(value.split()[0][1:-1])
            if type(value1) == str:
                value1 = float(value1.split()[0][1:-1])
            if key == key1 and key not in out.keys() and value != np.nan != value1:
                stack = np.hstack((value, value1))

                # optional keyword argument
                if dictionary3:

                    if type(dictionary3) == tuple:
                        dictionary3 = dict(dictionary3[0])
                    for key3, value3 in dictionary3.iteritems():
                        key3 = str(key3)
                        if type(value3) == str:
                            value3 = float(value3.split()[0][1:-1])
                        if key3 == key and value3 != np.nan:
                            appended_stack = np.hstack((stack, value3))
                            out[str(key)] = appended_stack
                else:
                    out[str(key)] = stack

    assert len(out) != 0  # combine_values found no matched keys!
    return out


def append_values(dictionary1, dictionary2, *dictionary3):
    """
    """
    out = {}
    # dictionary 1
    for key, value in dictionary1.iteritems():
        if type(value) == str:
            value = float(value.split()[0][1:-1])
        if key not in out.keys() and 'non' != value != np.nan:
            out[key] = value
    # dictionary 2
    for key2, value2 in dictionary2.iteritems():
        key2 = str(key2)
        if type(value2) == str:
            value2 = float(value1.split()[0][1:-1])
        if key2 not in out.keys() and 'non' != value2 != np.nan:
            out[key2] = value2

    # optional keyword argument
    if dictionary3:
        if type(dictionary3) == tuple:
            dictionary3 = dict(dictionary3[0])
        for key3, value3 in dictionary3.iteritems():
            key3 = str(key3)
            if type(value3) == str:
                value3 = float(value3.split()[0][1:-1])
            if key3 not in out.keys() and 'non' != value3 != np.nan:
                out[key3] = value3

    assert len(out) != 0  # combine_values found no matched keys!
    return out


def show_plot(title='', x_label='', y_label='', saveaspng='', save=True):
    """
    """
    assert type(title) == str
    assert type(x_label) == str
    assert type(y_label) == str

    plt.ylabel(y_label)
    plt.xlabel(x_label)
    plt.title(title)
    plt.show()

    if save:
        plt.savefig(saveaspng + '.png')


class plot():
    def __init__(self, indices_out_file, lamost=False, lowres=False, calstar=False):
        """
        input: path_to_ew: a text file, the first column is the kic numbers, the second column is a dictionary with
            each line measured as the keys and the value is the measurement and error of the line
        """

        # predefined lists:
        self.spectypes = ['K5', 'K7', 'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6']
        self.color = ['pink', 'red', 'orange', 'yellow', 'green', 'blue', 'maroon', 'black']

        # load in data
        ew_df = pd.DataFrame.from_csv(indices_out_file, sep='\t')
        ha = ew_df['ha']
        tio1 = ew_df['tio1']
        tio2 = ew_df['tio2']
        tio3 = ew_df['tio3']
        tio4 = ew_df['tio4']
        tio5 = ew_df['tio5']
        cah2 = ew_df['cah2']
        cah3 = ew_df['cah3']
        caoh = ew_df['caoh']
        o2 = ew_df['o2']
        n0 = ew_df['n0']
        self.linelist = [ha, tio1, tio2, tio3, tio4, tio5, cah2, cah3, caoh, o2, n0]
        self.linelist_name = ['ha', 'tio1', 'tio2', 'tio3', 'tio4', 'tio5', 'cah2', 'cah3', 'caoh', 'o2', 'n0']
        self.ha = ha
        self.tio5 = tio5
        self.cah2 = cah2
        self.cah3 = cah3

        cah = {}
        for starid, (cah2, cah3) in (combine_values(self.cah2, self.cah3)).iteritems():
            cah[starid] = float(cah2) + float(cah3)
        self.cah = cah

        if lamost:
            xvalue = {}
            self.sptype = ew_df['SpT']
            for id, type in self.sptype.iteritems():
                for spectype in self.spectypes:
                    if str(type) == str(spectype):
                        xvalue[id] = self.spectypes.index(spectype)
            self.spt_index = xvalue
        elif lowres:
            spt_df = pd.DataFrame.from_csv('/home/tessa/PycharmProjects/mdwarf/finalSpT.txt', sep='\t')
            spt = spt_df['spectype']
            sptdata = {}
            for kic, sptype in spt.iteritems():
                if str(sptype) != 'bad':
                    sptdata[int(kic)] = [sptype]
            self.sptype = sptdata

            xvalue = {}
            for kic, value in self.sptype.iteritems():
                for st in self.spectypes:
                    if st == value[0]:
                        xpos = self.spectypes.index(st)
                        xvalue[int(kic)] = xpos
            self.spt_index = xvalue
        elif calstar:
            spt_df = pd.DataFrame.from_csv('Master_Data_42216.csv')
            spt = spt_df['Spectral Types - Hammer output']
            xvalue = {}
            for name, sptype in spt.iteritems():
                for spt in self.spectypes:
                    if spt == sptype != 'non':
                        xpos = self.spectypes.index(spt)
                        xvalue[name] = xpos
            self.spt_index = xvalue

        if lamost:
            zeta_values = pd.read_csv("lamost_zeta.csv", delimiter=',', header=0, index_col=0)
        elif lowres:
            zeta_values = pd.read_csv("lowres_zeta.csv", delimiter=',', header=0, index_col=0)
        elif calstar:
            zeta_values = pd.read_csv("calstar_zeta.csv", delimiter=',', header=0, index_col=0)
        self.zeta = zeta_values['zeta']

        kic_values = pd.DataFrame.from_csv(
            '/home/tessa/astronomy/mdwarf/hammer_output/Copy of APO_MDwarf_properties 1015 - updated - May 25, 2015.csv',
            sep=',')
        self.prot_mc14 = kic_values['Prot_Mc14']

        gmag = kic_values['g Mag']
        # print kic_values.keys()
        kmag = kic_values['K Mag']
        mag = {}
        for gkey, g in gmag.iteritems():
            for kkey, k in kmag.iteritems():
                if gkey == kkey:
                    g_k = g - k
                    mag[gkey] = g_k
        self.g_k = mag

        if lamost:
            logll = pd.DataFrame.from_csv('lamost_mdwarf_logll.csv')
        elif lowres:
            logll = pd.DataFrame.from_csv('mdwarf_logll.csv')
        elif calstar:
            logll = pd.DataFrame.from_csv('calstar_logll.csv')
        self.logll = logll['llog']

        if lamost:
            feh = pd.DataFrame.from_csv('lamost_fe_h.csv')
        elif lowres:
            feh = pd.DataFrame.from_csv('lowres_fe_h.csv')
        elif calstar:
            feh = pd.DataFrame.from_csv('calstar_feh.csv')
        self.feh = feh['Fe/H']

    def define_zeta(self, outfilename, save=True):
        """
        """
        out = {}

        # solar data
        tio_cah = combine_values(self.tio5, self.cah)
        for key, (tio5, cah) in tio_cah.iteritems():
            solar_tio5 = (-0.050) - 0.164 * (cah ** 3.) + 0.670 * (cah ** 2.) - 0.118 * cah
            zeta = (1.0 - tio5) / (1.0 - solar_tio5)
            out[key] = zeta

        assert type(outfilename) == str
        assert len(out) != 0
        if save:
            df = pd.DataFrame(out.values(), index=out.keys(), columns=['zeta'])
            df.to_csv(outfilename + '.csv', sep=',')

    def define_logll(self, outfilename):
        """
        """

        chi = [10 ** 0, 10 ** (-3.93438), 10 ** (-4.01516), 10 ** (-4.13192), 10 ** (-4.19592), 10 ** (-4.56243),
               10 ** (-4.75570), 10 ** (-5.28066), 10 ** (-5.21965), 10 ** (-5.41719)]
        logll_dict = {}

        values = combine_values(self.spt_index, self.ha)
        for kic, (index, ha) in values.iteritems():
            logll = np.log10(float(ha) * float(chi[int(index)]))
            logll_dict[kic] = logll

        assert type(outfilename) == str
        assert len(logll_dict) != 0
        df = pd.DataFrame(logll_dict.values(), index=logll_dict.keys(), columns=['llog'])
        df.to_csv(outfilename + '.csv', sep=',')

    def define_feh(self, outfilename):
        """ derived from Mann 2012"""

        out_feh = {}
        # calculate fe/h per star
        for key, z in self.zeta.iteritems():
            feh = 1.55 * z - 1.62
            out_feh[key] = feh

        assert type(outfilename) == str
        assert len(out_feh) != 0
        df = pd.DataFrame(out_feh.values(), index=out_feh.keys(), columns=['Fe/H'])
        df.to_csv(outfilename + '.csv', sep=',')

    def _vs_(self, xvalue={}, yvalue={}, marker='*', xlog=False, ylog=False, xmagnitude=False, ymagnitude=False):
        """
        """

        xy = combine_values(xvalue, yvalue, self.spt_index)
        for starid, (x, y, index) in xy.iteritems():
            plt.scatter(x, y, marker=marker, color=self.color[int(index)], alpha=0.5)
        if ymagnitude:
            plt.gca().invert_yaxis()
        if xmagnitude:
            plt.gca().invert_xaxis()
        if ylog:
            plt.yscale('log')
        if xlog:
            plt.xscale('log')

    def linelist_vs_spt(self, yvalue1 = {}, *yvalue2 = {}, *yvalue3 = {}):
        ''' automatically shows plot
        '''

    for i in xrange(len(self.linelist)):
        y = combine_values(self.linelist[i], yvalue.spt_index)
        for k, (line, index) in y.iteritems():
            try:
                msmt = line[1:line.find(',')]
            except AttributeError:
                msmt = line
            if msmt > 0:
                plt.scatter(index, msmt, marker='o', alpha=0.6)
        if yvalue2:
            y = combine_values(self.linelist[i], yvalue.spt_index)
            for k, (line, index) in y.iteritems():
                try:
                    msmt = line[1:line.find(',')]
                except AttributeError:
                    msmt = line
                if msmt > 0:
                    plt.scatter(index, msmt, marker='o', alpha=0.6)
        if magnitude:
            plt.gca().invert_yaxis()
        plt.ylabel(self.linelist_name[i])
        plt.xticks(np.arange(len(self.spectypes)), self.spectypes)
        #plt.margins(0.2)
        # if lamost:
        #     show_plot('LAMOST ' + self.linelist_name[i] + ' vs SpT', 'Spectral Types', self.linelist_name[i], \
        #               saveaspng='LAMOST_' + str(self.linelist_name[i]) + '.png')
        # elif lowres:
        #     show_plot('Lowres ' + self.linelist_name[i] + ' vs SpT', 'Spectral Types', self.linelist_name[i], \
        #               saveaspng='lowres_' + str(self.linelist_name[i]) + '.png')
        # elif highres:
        #     show_plot('Highres ' + self.linelist_name[i] + ' vs SpT', 'Spectral Types', self.linelist_name[i], \
        #               saveaspng='highres' + str(self.linelist_name[i]) + '.png')
        # elif calstar:
        #     show_plot('Calibration Stars ' + self.linelist_name[i]+ ' vs SpT', 'Spectral Types', self.linelist_name[i],\
        #               saveaspng='calstar' + str(self.linelist_name[i]) + '.png')

    def _vs_spt(self, yvalue={}, marker='*', color='blue', log=False, magnitude=False):

        xy = combine_values(self.spt_index, yvalue)
        for key, (index, y) in xy.iteritems():
            plt.scatter(index, y, marker=marker, alpha=0.6, color=color)
        if log:
            plt.yscale('log')
        if magnitude:
            plt.gca().invert_yaxis()
        plt.xticks(np.arange(len(self.spectypes)), self.spectypes)  # and i[0] are the sorted keys
        plt.margins(0.2)

    def histogram(self, dictionary, log=False, spt=False):

        array = []

        for spt_index in xrange(len(self.spectypes)):
            self.spectypes[spt_index] = []
            count = 0

            for k, index in dictionary.iteritems():
                # if log:
                #     array.append(np.log(index))
                if spt_index == index:
                    count += 1
                self.spectypes[spt_index] = count
            array.append(self.spectypes[spt_index])

        print array
        assert len(array) != 0

        # plt.hist(array,  align = 'mid')
        plt.bar(range(len(self.spectypes)), array, align='center')
        if log:
            plt.xscale('log')
        if spt:
            plt.xticks(np.arange(len(self.spectypes)), self.spectypes)
            # plt.margins(0.2)


#############################################################################################
#############################################################################################

lowres_filepath = 'Mdwarf_lowres_outfile.txt'
highres_filepath = 'finalEW.txt'
standard_filepath = 'calstar_outfile.txt'
lamost = 'lamost_lowres_outfile.txt'

# CHANGE as needed
prefix_title = 'All_data_lowres_'

l = plot(lamost, lamost=True)
m = plot(lowres_filepath, lowres=True)
s = plot(standard_filepath, calstar=True)

# define values (need to comment out in _init_ to define self.values)
# m.define_zeta('lowres_zeta')
# m.define_logll('lowres_logll')
# m.define_feh('lowres_feh')
# s.define_zeta('calstar_zeta')
# s.define_logll('calstar_logll')
# s.define_feh('calstar_feh')




# linelist vs spectral type
# l.linevs_spt(linelist=True, lamost = True)
# m._vs_spt(linelist=True, lowres = True)
s._vs_spt(linelist=True,calstar = True )


# histogram of data counts
mls = append_values(m.spt_index, l.spt_index, s.spt_index)
m.histogram(mls, spt=True)
# l.histogram(l.spt_index, spt = True)
show_plot(title=prefix_title + ' M Dwarfs per SpT', y_label='count', x_label='Spectral Types',
          saveaspng=str(prefix_title) + '_spt')


# # log LHa / Lbol
m._vs_spt(linelist=False, yvalue=m.logll, marker='o', color='blue')
l._vs_spt(linelist=False, yvalue=l.logll, marker='s', color='green')
s._vs_spt(linelist=False, yvalue=s.logll, marker='^', color='red')
show_plot(title=prefix_title + ' vs spt', y_label='log L Ha / L bol', saveaspng=prefix_title + '_logll_spt')




# # rotational period
# m._vs_spt(linelist=False, yvalue=m.prot_mc14, log=True)
# show_plot(title= prefix_title + ' Rotational Period vs Spectral Type',y_label='Prot')
#
# m._vs_(m.prot_mc14, m.logll, xlog=True)
# show_plot(title= prefix_title + ' logll vs Prot', x_label='Prot', y_label='log L Ha / L bol')
#
# m._vs_(m.prot_mc14, m.feh)
# show_plot(title = prefix_title + 'FeH vs Prot', xlog = True, x_label= ' Prot', y_label= ' Fe / H ')
#
# m.histogram(m.prot_mc14, log=True, title='lowres Prot Histogram', x_label='log Rotational Period')

# m._vs_(m.prot_mc14, m.zeta)
# show_plot(title= prefix_title + ' Zeta vs Prot', xlog=True, x_label='Prot', y_label='Zeta')



# tio5 vs zeta
m._vs_(m.tio5, m.zeta, marker='s')
l._vs_(l.tio5, l.zeta, marker='o')
s._vs_(s.tio5, s.zeta, marker='^')
show_plot(title=prefix_title + ' Zeta vs Tio5', x_label='Tio5', y_label='Zeta', saveaspng=prefix_title + '_zeta_tio5')


# # feh vs zeta
m._vs_(m.feh, m.zeta, xlog=True, marker='s')
l._vs_(l.feh, l.zeta, xlog=True, marker='o')
s._vs_(s.feh, s.zeta, xlog=True, marker='^')
show_plot(title=prefix_title + 'FeH vs Zeta', x_label='log Fe / H', y_label='Zeta',
          saveaspng=prefix_title + '_feh_zeta')

# zeta vs spt
m._vs_spt(linelist=False, yvalue=m.zeta, marker='s', color='blue')
l._vs_spt(linelist=False, yvalue=l.zeta, marker='o', color='green')
s._vs_spt(linelist=False, yvalue=s.zeta, marker='^', color='red')
show_plot(title=prefix_title + 'zeta vs SpT', y_label='Zeta', saveaspng=prefix_title + '_zeta_spt')


# cah vs tio5
m._vs_(m.tio5, m.cah, marker='s')
l._vs_(l.tio5, l.cah, marker='o')
s._vs_(s.tio5, s.cah, marker='^')
show_plot(title=prefix_title + 'CaH vs TiO5', y_label='CaH2 + CaH3', x_label='TiO5',
          saveaspng=prefix_title + '_cah_tio5')




# color
# m._vs_(m.tio5, m.g_k, ymagnitude=True)
# show_plot(title= prefix_title + ' Tio5 vs g-K', x_label='Tio5', y_label='g - K', saveaspng= prefix_title+'_tio5_gk')
#
# m._vs_(m.logll, m.g_k, ymagnitude=True)
# show_plot(title= prefix_title + 'logll vs g-K', x_label='log L Ha / L bol', y_label='g - K', saveaspng= prefix_title+'_logll_gk')
#
# m._vs_(m.cah2, m.g_k, ymagnitude=True)
# show_plot(title= prefix_title + 'cah2 vs g-K', x_label='CaH2', y_label='g-K',saveaspng= prefix_title+'_cah2_gk' )
#
# m._vs_spt(linelist=False, yvalue=m.g_k, magnitude=True)
# show_plot(title= prefix_title + 'g-K vs SpT', y_label='g-K', saveaspng= prefix_title+'_gk_spt')
