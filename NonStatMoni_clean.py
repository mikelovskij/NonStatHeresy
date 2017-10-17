#!/usr/bin/env python
# @author Gabriele Vajente
# /// Report generator for NonStatMoni BRMS monitor

import time
import numpy as np

from argparse import ArgumentParser
import markup
from virgotools import gps2str
import tables as tb
import matplotlib.gridspec as grsp
from functions import brms_reader, ipsh, extractbands, Parameters
import os
import fnmatch
from data_processing import DataProcessing
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

start = time.time()
# **************************** Initialization ******************************
ntop = 5
max_length = 10000

# directories where plots and html files must be saved
pdir = 'plots/'

style = 'general.css'
tabstyle = 'table.css'
source = ''
# brmsfile = '/users/valentin/PycharmProjects/Long Spectra/psd statistics 1186876818 1187222418 cutfreq=4000 resolution=2.5/brms.hdf5'
brmsfile = '/users/valentin/PycharmProjects/Long Spectra/psd statistics 1187913618 1187928018 cutfreq=4000 resolution=2.5/brms.hdf5'
hdir = '/users/valentin/PycharmProjects/Nonstatmoni_reader/Results/{}/'.format(
    brmsfile.split('/')[-2])

try:
    os.mkdir(hdir)
except OSError:
    print 'Results Directory already exists'
try:
    os.mkdir(hdir + pdir)
except OSError:
    print 'Plots Directory already exists'


# set up configuration file


# read configuration options and process
parser = ArgumentParser()
parser.add_argument("-i", "--init", dest="initialization",
                    default='NonStatMoni.ini',
                    help="set initialization FILE", metavar="cfgFILE")
parser.add_argument("-b", "--brms", dest='brms_file', default=brmsfile,
                    help="set brms data file", metavar="brmsFILE")
parser.add_argument("-d", "--dir", dest="hdir",
                    help="output directory", metavar="OutDir")

args = vars(parser.parse_args())
par = Parameters(args['initialization'])
# TODO: non e' piu' semplice semplicemente settare il default per opt.hdir?
if args['hdir']:
    hdir = args['hdir']

# ##### Loop over each channel ##############################
# Rebuild the band list as a list
for group, g_dict in par.group_dict.iteritems():
    extractbands(g_dict)
# Read the BRMS data stored in the file.
gpsb, gpse, fs, segments, times = brms_reader(args['brms_file'], par.group_dict)
par.extract_aux_channels(gpsb)
print "Analyzing lock between %d and %d" % (gpsb, gpse)
# TODO: trasformo questo pezzetto in una funzione o perdo in leggibilita'?
if max(times) <= 300:
    tunits = 's'
    t = times
if (max(times) > 300) and (max(times) <= 60 * 100):
    t = times / 60
    tunits = 'min'
if max(times) > (60 * 100):
    t = times / 3600
    tunits = 'h'


# TODO: - che cazzo e' dt, non fa niente!

# TODO: - capire come plottare timeseries e psd dei gruppi nella nuova struttura




processor = DataProcessing(segments, fs, par.nav, par.group_dict, times)
processor.cumulative_psd_computation()
for (aux_name, aux_groups), i in zip(par.aux_dict.iteritems(),
                                     xrange(len(par.aux_dict))):
    processor.auxillary_psd_csd_correlation(aux_name, aux_groups,
                                            par.aux_source)

processor.pearson_cohefficient_computation()
processor.coherence_computation()

ipsh()

for channel, units, dt, nav, cohe, nchan in zip(par.channel, par.units, par.dt,
                                                par.nav,
                                                par.coherences,
                                                xrange(len(par.channel))):
    band_data = brms_data[channel + '_' + str(nchan)]
    bands = channelsandbands[channel + '_' + str(nchan)]
    print("Channel " + str(nchan + 1) + " of " + str(len(par.channel)))
    for b in bands:
        chbandsall.append('%s_%d@%s' % (channel, nchan, b))

    # ####### Plot time series
    fig = plt.figure(figsize=(10, 6))
    plt.axes([0.1, 0.1, 0.65, 0.8])
    for y in band_data.itervalues():
        plt.semilogy(t[::50], y[::50], linewidth=0.5)

    plt.grid(True)
    plt.axis(xmax=max(t))
    plt.xlabel('Time [' + tunits + ']')
    plt.ylabel('Normalized BRMS [' + units + '/rHz]')
    plt.legend(bands, loc=(1.02,
                           0.2))  # Todo: uso invece le bande in band_data, o almeno faccio un confronto per vedere se van bene?
    plt.title('BRMS %s - %s' % (gps2str(gpsb), gps2str(gpse)), fontsize=10)
    bmin = bands[0].split('_')[0]
    bmax = bands[-1].split('_')[1].split('Hz')[0]
    plt.savefig(hdir + pdir + channel + '_' + str(
        nchan) + '_' + dt + 's_' + bmin + '_' + bmax + 'Hz_time.png')
    plt.close()
    # ####### Compute and plot spectra
 # all bands should have the same length, according to 9 out of 10 statisticians

    sp = []
    cohs = []
    ccfs = []
    histograms = []
    seg_mask = np.ones(len(segments), dtype=bool)
    for j in xrange(len(segments) - 1):
        if (segments[j][1] == segments[j + 1][0]) and (
                    (segments[j + 1][1] - segments[j][0]) < max_length):
            segments[j + 1][0] = segments[j][0]
            seg_mask[j] = False
    segments = segments[seg_mask]
    for _ in cohe:
        cohs.append([])
        ccfs.append([])
        histograms.append([])
    data_storage = None
    for y in band_data.itervalues():
        for psd, csd, aux_mean, square_aux_mean, prod_mean, h, k in zip(s2,
                                                                        csds,
                                                                        aux_means,
                                                                        square_aux_means,
                                                                        prod_means,
                                                                        hist,
                                                                        xrange(
                                                                            len(
                                                                                s2))):
            cohs[k].append((np.abs(csd)) ** 2 / (psd * s1))
            ccfs[k].append((prod_mean - mean_1 * aux_mean) / np.sqrt(
                (square_mean_1 - mean_1 ** 2) * (
                    square_aux_mean - aux_mean ** 2)))
            histograms[k].append(h)
    plt.figure(figsize=(10, 6))
    plt.axes([0.1, 0.1, 0.65, 0.8])
    for s in sp:
        try:
            plt.loglog(f, np.sqrt(s), linewidth=0.5)
        except ValueError:
            print "Issue in plotting psd of channel {0} {1}, may there be no positive data?".format(
                channel, s)
    plt.axis(xmax=0.5)
    plt.axis(xmin=1 / (2 * processor.n_points))
    plt.grid(True)
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Spectrum of normalized BRMS [' + units + '/rHz/rHz]')
    plt.legend(bands, loc=(1.02, 0.2))
    plt.title('BRMS %s - %s' % (gps2str(gpsb), gps2str(gpse)), fontsize=10)
    plt.savefig(hdir + pdir + channel + '_' + str(
        nchan) + '_' + dt + 's_' + bmin + '_' + bmax + 'Hz_psd.png')
    plt.close()

    for C, corrcoef, hist, cname, naux in zip(cohs, ccfs, histograms, cohe,
                                              xrange(len(cohs))):
        # C = []
        nplot = 0
        hicorrlist = []
        hicorrbands = []

        # do this for the first main channel, for each aux channel
        if nchan == 0:
            meancoh.append([])
            ccf.append([])

        for b, j in zip(band_data.iteritems(), xrange(len(band_data.items()))):
            meancoh[naux].append(np.mean(
                C[j]))  # TODO: Isn't the maximum coherence more interesting?
            # TODO: use better and clearer names for the variables.
            ccf[naux].append(corrcoef[j])
            # TODO: surely some distributions can give a significance estimate of the corrcoef
            if abs(corrcoef[j]) >= 0.15:
                nplot += 1
                hicorrlist.append(hist[j])
                hicorrbands.append(b)
        plt.figure(figsize=(15, 6))
        # axes([0.1, 0.1, 0.65, 0.8])
        if nplot > 0:
            m = int(np.floor(np.sqrt(nplot)))
            n = int(np.ceil(np.sqrt(nplot)))
            gs = grsp.GridSpec(n, 2 * m)
            gs.update(wspace=0.2, hspace=0.5)
            for j in xrange(n):
                for (k, h, b) in zip(xrange(m), hicorrlist[m * j:m * (j + 1)],
                                     hicorrbands[m * j:m * (j + 1)]):
                    ax = plt.subplot(gs[j, m + k])
                    a = ax.pcolormesh(h[1], h[2], h[0])
                    ax.tick_params(axis='both', which='minor', labelsize=6)
                    ax.tick_params(axis='both', which='major', labelsize=6)
                    ax.ticklabel_format(style="sci", scilimits=(0, 0),
                                        axis="both")
                    plt.title(b[0])
            ax1 = plt.subplot(gs[0:n, 0:m])

        else:
            ax1 = plt.subplot2grid((1, 1), (0, 0))
        for c in C:
            try:
                ax1.semilogx(f, c, linewidth=0.5)
                ax1.axis(ymin=0, ymax=1)
                ax1.axis(xmin=1 / (2 * processor.n_points))
                ax1.grid(True)
                ax1.set_xlabel('Frequency [Hz]')
                ax1.set_ylabel('Coherence')
                ax1.legend(bands, loc=(1.02, 0.2))
            except ValueError:
                print "Issue in plotting psd of channel {0} {1}, may there be no positive data?".format(
                    cname, channel)
        plt.title('Coherence with ' + cname + ' %d - %d' % (gpsb, gpse))
        plt.savefig(hdir + pdir + channel + '_' + str(
            nchan) + '_' + dt + 's_' + bmin + '_' + bmax + 'Hz_cohe_' + cname + '.png')
        plt.close()

# ##### Create summary web page ##############################################
print "Generating Web Page..."
page = markup.page()
page.init(title="NonStatMoni", css=(style, tabstyle),
          footer="(2017)  <a href=mailto:snao20@hotmail.com>Michele Valentini</a>")
page.a(name="top")

# brms summary table generation
ch = [par.channel[0]]
num = 0
chnum = []
# obtain an unique channel list, and a vector containing how many times each channel is repeated
# works only if channels are listed "in order"
for channel in par.channel:
    if channel != ch[-1]:
        ch.append(channel)
        chnum.append(num)
        num = 1
    else:
        num += 1
chnum.append(num)

page.table(id="t01")
prenum = 0
for c, num in zip(ch, chnum):
    page.th()
    page.td(style="width:50%")
    page.h4(c)
    page.th.close()
    for channel, dt, n in zip(par.channel[prenum: prenum + num],
                              par.dt[prenum: prenum + num], xrange(num)):
        bands = channelsandbands[channel + '_' + str(prenum + n)]
        page.td(style="width:100%")
        bmin = bands[0].split('_')[0]
        bmax = bands[-1].split('_')[1].split('Hz')[0]
        chbands = bmin + '_' + bmax
        page.add(
            '<a href=#%s>%s</a><br>(dt=%s)' % (channel + chbands, chbands, dt))
        page.td.close()
    page.tr.close()
    prenum += num
page.table.close()

page.br()

# best coherences and ccfs summary table generation
# TODO: this should be inside the channel loop! maybe, kind of (atm uses only the aux_channels of the last loop...)
page.table(id="t01")
for (cname, naux) in zip(cohe, xrange(len(cohs))):
    page.tr()
    page.td()
    page.h4(cname)
    page.td.close()
    ccfsort = np.argsort(np.abs(ccf[naux]))[::-1]
    mncohsort = np.argsort(meancoh[naux])[::-1]
    for chan in ccfsort[0:4]:
        pagename = 'cohe_%s_Hz.html' % chbandsall[chan].split('@')[0]
        page.td(style="width:10%; font-size:11px;")

        page.add(
            "<a target=_blank href={}{}>{}</a><br>(ccf={:3f})".format(
                hdir.replace(' ', '%20'), pagename,
                chbandsall[chan], ccf[naux][chan]))
        page.td.close()

    page.td()
    page.add('<a>           ********         </a>')
    page.td.close()

    for chan in mncohsort[0:4]:
        #   (channel + '%s_%d@%s' + b % (channel, nchan, b))
        pagename = 'cohe_%s_Hz.html' % chbandsall[chan].split('@')[0]
        page.td(style="width:10%; font-size:11px;")
        # page.add(
        #     '<a target=_blank href=http://www.virgo.infn.it/MonitoringWeb/Noise/NonStatMoni/%s>%s</a><br>(mncoh=%3f)' % (
        #         pagename, chbandsall[chan], meancoh[naux][chan]))
        page.add(
            '<a target=_blank href=%s/%s>%s</a><br>(mncoh=%3f)' % (
                hdir.replace(' ', '%20'),
                pagename, chbandsall[chan], meancoh[naux][chan]))
        page.td.close()

    page.tr.close()
page.table.close()
# rest of the page generation
nbands = 0
for (channel, dt, cohe, nchan) in zip(par.channel, par.dt, par.coherences,
                                      xrange(len(par.channel))):
    bands = channelsandbands[channel + '_' + str(nchan)]
    bmin = bands[0].split('_')[0]
    bmax = bands[-1].split('_')[1].split('Hz')[0]

    page.h1("NonStatMoni BRMS for " + channel + " GPS %d - %d"
            % (gpsb, gpse))

    ccftab = np.zeros(ntop)
    cohtab = np.zeros(ntop)
    idccftab = np.zeros(ntop, dtype='i')
    idcohtab = np.zeros(ntop, dtype='i')
    for naux in xrange(len(cohs)):
        #  best coh /ccf per channel table
        for (b, corr, coh, i) in zip(bands,
                                     ccf[naux][nbands:(nbands + len(bands))],
                                     meancoh[naux][
                                     nbands:(nbands + len(bands))],
                                     xrange(len(bands))):
            if corr > min(ccftab):
                tccftab = ccftab
                tidccftab = idccftab
                tcorr = np.concatenate((tccftab, [corr]))
                tidccf = np.concatenate((tidccftab, [naux]))
                ii = tcorr.argsort()
                ii = ii[1:]
                ii = ii[::-1]
                ccftab = tcorr[ii]
                idccftab = tidccf[ii]
            if coh > min(cohtab):
                tcohtab = cohtab
                tidcohtab = idcohtab
                tcoh = np.concatenate((tcohtab, [coh]))
                tidcoh = np.concatenate((tidcohtab, [naux]))
                ii = tcoh.argsort()
                ii = ii[1:]
                ii = ii[::-1]
                cohtab = tcoh[ii]
                idcohtab = tidcoh[ii]

    nbands += len(bands)

    page.h2("Normalized BRMS time series")
    chbands = channel + bmin + '_' + bmax
    page.a(name=chbands)
    page.a("Top", href="#top")

    tab = []
    row_0 = [" CCF ", "Coherence"]
    tab.append(row_0)
    for i in xrange(ntop):
        row = [
            "{0}<br>Highest CCFs = {1}".format(cohe[idccftab[i]], ccftab[i]),
            "{0}<br>Highest Coher. = {1}".format(cohe[idcohtab[i]], cohtab[i])]
        tab.append(row)
    tabstr = tb.tabmaker(tab, True, False)

    page.div(
        "<div style=\"float:left\"><img src=\"" + pdir + channel + '_' + str(
            nchan) + '_' + dt + "s_" + bmin + "_" + bmax + "Hz_time.png\" /></div><div>" + tabstr + "</div>")

    page.h2("Spectrum of BRMS time series")
    page.img(src=(pdir + channel + '_' + str(
        nchan) + '_' + dt + 's_' + bmin + '_' + bmax + 'Hz_psd.png'),
             alt="Plots")
    page.h2("<a href=#top> ^ Top ^ </a> ")
    # create coherence subpage
    page2 = markup.page()
    page2.init(title="NonStatMoni", css='../style/style.css')
    page2.h1("Coherences with slow channels")
    for cname in cohe:
        page2.img(src=(pdir + channel + '_' + str(
            nchan) + '_' + dt + 's_' + bmin + '_' + bmax + 'Hz_cohe_'
                       + cname + '.png'), alt="Plots")
    pagename = 'cohe_%s_%d_Hz.html' % (channel, nchan)
    # 'cohe_' + channel + '_' + dt + 's_' + bmin + '_' + bmax + 'Hz.html'
    page2.savehtml(hdir + pagename)
    page.h2(
        "<a target=_blank href=" + hdir.replace(' ',
                                                '%20') + pagename + ">Coherences with slow channels</a>")
    page.br()

page.br()
page.h2("Contacts")

page.savehtml(hdir + 'index.html')

print "Elapsed time %d seconds" % int(time.time() - start)
