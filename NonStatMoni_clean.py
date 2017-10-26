#!/usr/bin/env python
# @author Gabriele Vajente
# /// Report generator for NonStatMoni BRMS monitor

import time
import numpy as np

from argparse import ArgumentParser
import markup
from markup import oneliner as ol
from virgotools import gps2str
import tables as tb
import matplotlib.gridspec as grsp
from functions import brms_reader, ipsh, extractbands, Parameters
import os
import fnmatch
from data_processing import DataProcessing
import matplotlib
from multiprocessing import Pool

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from collections import OrderedDict

start = time.time()
# **************************** Initialization ******************************
ntop = 5


# directories where plots and html files must be saved
pdir = 'plots/'
style = 'general.css'
tabstyle = 'table.css'
modalstyle = 'modal.css'
script = 'scripts.js'
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
parser.add_argument("-p", "--proc", dest="n_proc", metavar="NumberOfProcesses",
                    help='number of parallel processes to run', default=4)

args = vars(parser.parse_args())
par = Parameters(args['initialization'])
# TODO: non e' piu' semplice semplicemente settare il default per opt.hdir?
if args['hdir']:
    hdir = args['hdir']

# ##### Loop over each channel ##############################
# Rebuild the band list as a list
for _, g_dict in par.group_dict.iteritems():
    extractbands(g_dict)
# Read the BRMS data stored in the file.
gpsb, gpse, fs, segments, times = brms_reader(args['brms_file'],
                                              par.group_dict)
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


# Coalesce the segments if they are too short.
approx_fft_duration = int(np.floor(float(gpse - gpsb) / par.nav))
seg_mask = np.ones(len(segments), dtype=bool)
for j in xrange(len(segments) - 1):
    if (segments[j][1] == segments[j + 1][0]) and (
                (segments[j + 1][1] - segments[j][0]) < (
                10 * approx_fft_duration)):
        segments[j + 1][0] = segments[j][0]
        seg_mask[j] = False
segments = segments[seg_mask]

proc = DataProcessing(segments, fs, par.nav, par.group_dict, times)
proc.cumulative_psd_computation()

pool = Pool(processes=args['n_proc'])
pool_args = []
for (aux_name, aux_groups), i in zip(par.aux_dict.iteritems(),
                                     xrange(len(par.aux_dict))):
    #proc.auxillary_psd_csd_correlation(aux_name, aux_groups, par.aux_source)
    pool_args.append([aux_name, aux_groups, par.aux_source])

pool.map(proc.auxillary_psd_csd_correlation, pool_args)
proc.pearson_cohefficient_computation()
proc.coherence_computation()


# for channel, units, dt, nav, cohe, nchan in zip(par.channel, par.units, par.dt,
#                                                par.nav,
#                                                par.coherences,
#                                                xrange(len(par.channel))):
for group, g_dict in par.group_dict.iteritems():
    band_data = g_dict['brms_data']
    # bands = channelsandbands[channel + '_' + str(nchan)]
    # for b in bands:
    #   chbandsall.append('%s_%d@%s' % (channel, nchan, b))

    # ####### Plot time series
    fig = plt.figure(figsize=(10, 6))
    plt.axes([0.1, 0.1, 0.65, 0.8])
    for band_name, data in band_data.iteritems():
        plt.semilogy(t[::50], data[::50], linewidth=0.5, label=band_name)
    plt.grid(True)
    plt.axis(xmax=max(t))
    plt.xlabel('Time [' + tunits + ']')
    plt.ylabel('Normalized BRMS [' + g_dict['units'] + '/rHz]')
    plt.legend(loc=(1.02, 0.2))
    plt.title('{} BRMS {} - {}'.format(g_dict['channel'], gps2str(gpsb),
                                       gps2str(gpse)), fontsize=10)
    plt.savefig(hdir + pdir + group + '_time.png')
    plt.close()

    # ####### Plot  band spectra
    plt.figure(figsize=(10, 6))
    plt.axes([0.1, 0.1, 0.65, 0.8])
    for band, sp in proc.brms_psd[group].iteritems():
        try:
            plt.loglog(proc.freqs, np.sqrt(sp), linewidth=0.5, label=band)
        except ValueError:
            print "Issue in plotting psd of {0}_{1}," \
                  " may there be no positive data?".format(g_dict['channel'],
                                                           band)
    plt.axis(xmax=0.5)
    plt.axis(xmin=1 / (2 * proc.n_points))
    plt.grid(True)
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Spectrum of normalized BRMS [' + g_dict['units'] + '/rHz/rHz]')
    plt.legend(loc=(1.02, 0.2))
    plt.title('{} BRMS {} - {}'.format(g_dict['channel'], gps2str(gpsb),
                                       gps2str(gpse)), fontsize=10)
    plt.savefig(hdir + pdir + group + '_psd.png')
    plt.close()

    # for C, corrcoef, hist, cname, naux in zip(cohs, ccfs, histograms, cohe,
    #                                          xrange(len(cohs))):
    for aux_name, aux_groups in par.aux_dict.iteritems():
        if group in aux_groups:
            aux_results = proc.aux_results[aux_name]
            # hicorrlist = []
            hi_corr = []
            # hicohlist = []
            hi_coh = []
            # plotting every coherence is probably too much.
            #  fo una cosa tipo coherence?
            for (band_name, data), j in zip(band_data.iteritems(),
                                               xrange(len(band_data.items()))):
                ccf = proc.ccfs[aux_name][group + '_' + band_name]
                mean_coh = proc.mean_cohs[aux_name][
                    group + '_' + band_name]
                # TODO: surely some distributions can give a significance estimate of the corrcoef
                if abs(ccf) >= 0.15:
                    # hicorrlist.append(hist[j])
                    hi_corr.append((band_name, j))
                if np.mean(mean_coh) >= 0.03:
                    # hichlist
                    hi_coh.append((band_name, j))
            plt_path = hdir + pdir + group + '_cohe_' + aux_name.split(':')[1] + '.png'
            if (len(hi_corr) + len(hi_coh)) > 0:
                plt.figure(figsize=(15, 6))
                # axes([0.1, 0.1, 0.65, 0.8])
                if len(hi_corr) > 0:
                    n_rows = int(np.floor(np.sqrt(len(hi_corr))))
                    n_cols = int(np.ceil(np.sqrt(len(hi_corr))))
                    gs = grsp.GridSpec(n_cols, 2 * n_rows)
                    gs.update(wspace=0.2, hspace=0.5)
                    for j in xrange(n_cols):
                        for k, (band_name, band_num) in zip(xrange(n_rows),
                                                            hi_corr[
                                                            n_rows * j:n_rows * (
                                                                j + 1)]):
                            h = proc.aux_results[aux_name]['histogram'][band_num]
                            ax = plt.subplot(gs[j, n_rows + k])
                            a = ax.pcolormesh(h[1], h[2], h[0])
                            ax.tick_params(axis='both', which='minor', labelsize=6)
                            ax.tick_params(axis='both', which='major', labelsize=6)
                            ax.ticklabel_format(style="sci", scilimits=(0, 0),
                                                axis="both")
                            plt.title(band_name)
                    ax1 = plt.subplot(gs[0:n_cols, 0:n_rows])

                else:
                    ax1 = plt.subplot2grid((1, 1), (0, 0))
                for (band_name, band_num) in hi_coh:
                    try:
                        coh = proc.cohs[aux_name][group + '_' + band_name]
                        ax1.semilogx(proc.freqs, coh, linewidth=0.5,
                                     label=band_name)
                        ax1.axis(ymin=0, ymax=1)
                        ax1.axis(xmin=1 / (2 * proc.n_points))
                        ax1.grid(True)
                        ax1.set_xlabel('Frequency [Hz]')
                        ax1.set_ylabel('Coherence')
                        ax1.legend(loc=(1.02, 0.2))
                    except ValueError:
                        print "Issue in plotting cpherence of channel {0} _ {1}" \
                              ", could there be no positive data?" \
                            .format(g_dict['channel'], band_name)
                ax1.set_title('Coherence with {}, {:d} - {:d}'.
                              format(aux_name, gpsb, gpse))
                plt.savefig(plt_path)
                plt.close()
            else:
                try:
                    os.remove(plt_path)
                except OSError:
                    pass


# ##### Create summary web page ##############################################
print "Generating Web Page..."
page = markup.page()
page.init(title="NonStatMoni", css=(style, tabstyle, modalstyle),
          footer="(2017)" +
                 ol.a("Michele Valentini", href='mailto:snao20@hotmail.com'))

# best coherences and ccfs summary table generation

page.button("Auxillary channels summary table", class_="accordion")
page.div(class_="panel")
page.input(type="text", id="myInput", onkeyup="myFunction()",
           placeholder="Search for aux_channel names..")
page.table(id="t01")
for aux_name, aux_groups in par.aux_dict.iteritems():
    page.tr()
    page.td()
    page.h4(aux_name)
    page.td.close()
    # Sort the dictionary according to the abs of the values,
    #  from highest to lowest
    ccf_sorted = OrderedDict(sorted(proc.ccfs[aux_name].iteritems(),
                                    key=lambda x: abs(x[1]), reverse=True))
    mean_coh_sorted = OrderedDict(sorted(proc.mean_cohs[aux_name].iteritems(),
                                         key=lambda x: x[1], reverse=True))

    for chan_band, ccf in ccf_sorted.items()[0:4]:
        pagename = 'cohe_{}.html'.format(chan_band.split('_')[0])
        page.td(style="width:10%; font-size :11px;")

        page.add("<a target=_blank href={}>{}</a><br>(ccf={:3f})"
                 .format(pagename, chan_band, ccf))
        page.td.close()

    page.td()
    page.add('<a>      ********     </a>')
    page.td.close()

    for chan_band, coh in mean_coh_sorted.items()[0:4]:
        pagename = 'cohe_{}.html'.format(chan_band.split('_')[0])
        page.td(style="width:10%; font-size:11px;")
        page.add('<a target=_blank href={}>{}</a><br>(mncoh={:3f})'
                 .format(pagename, chan_band, coh))
        page.td.close()

    page.tr.close()
page.table.close()
page.div.close()

# build the page menu
onclick_gen = ("openGroup(event, '{}')".format(group)
               for group in par.group_dict.keys())
page.div(ol.button(par.group_dict.keys(), class_='tablinks',
                   onclick=onclick_gen), class_='tab')

# build each group subpage
for group, g_dict in par.group_dict.iteritems():
    # Build the highest ccf and coherence for this group table
    ccftab = np.zeros(ntop)
    cohtab = np.zeros(ntop)
    ccftab_names = np.zeros(ntop, dtype='string')
    cohtab_names = np.zeros(ntop, dtype='string')
    for aux_name, aux_groups in par.aux_dict.iteritems():
        if group in aux_groups:
            for band in g_dict['band_list']:
                ccf = proc.ccfs[aux_name][group + '_' + band]
                if abs(ccf) > min(np.abs(ccftab)):
                    ccftab = np.concatenate((ccftab, [ccf]))
                    ccftab_names = np.concatenate((ccftab_names, [aux_name]))
                    best_indexes = np.abs(ccftab).argsort()[1:][::-1]
                    ccftab = ccftab[best_indexes]
                    ccftab_names = ccftab_names[best_indexes]
                coh = proc.mean_cohs[aux_name][group + '_' + band]
                if coh > min(cohtab):
                    cohtab = np.concatenate((cohtab, [coh]))
                    cohtab_names = np.concatenate((cohtab_names, [aux_name]))
                    best_indexes = cohtab.argsort()[1:][::-1]
                    cohtab = cohtab[best_indexes]
                    cohtab_names = cohtab_names[best_indexes]
    tab = [[" CCF ", "Coherence"]]
    for i in xrange(ntop):
        row = [
            "{0}<br>Highest CCFs = {1}".format(ccftab_names[i], ccftab[i]),
            "{0}<br>Mean Coher. = {1}".format(cohtab_names[i], cohtab[i])]
        tab.append(row)
    tabstr = tb.tabmaker(tab, True, False)

    # build the rest of the frame
    frame = ol.h1("NonStatMoni BRMS for {} GPS {:d} - {:d}".format(
                  g_dict['channel'], gpsb, gpse))
    frame += ol.h2("Normalized BRMS time series")
    # todo: normalized by what?
    img = ol.img(class_="myImg", src=pdir + group + "_time.png",
                 alt=group + "Time plot", width="400")
    frame += ol.div(img, style="float:left")
    frame += ol.div(tabstr)
    frame += ol.h2("Spectrum of BRMS time series")
    frame += ol.img(class_="myImg", src=pdir + group + "_psd.png",
                    alt=group + "PSD plot", width="400")
    cohe_page_name = 'cohe_{}.html'.format(group)
    frame += ol.h2(ol.a("Coherences with slow channels", target='_blank',
                        href=cohe_page_name))
    page.div(frame, id=group, class_='tabcontent')

    # create coherence subpage
    page2 = markup.page()
    page2.init(title="NonStatMoni", css='../style/style.css')
    page2.h1("Coherences with slow channels")
    for aux_name, aux_groups in par.aux_dict.iteritems():
        if group in aux_groups:
            page2.img(src=(pdir + group + '_' + 'cohe_' + aux_name.split(':')[1] + '.png'),
                      alt=(" No plots for" + aux_name))
    page2.savehtml(hdir + cohe_page_name)

# create the modal for the plot images
modal = ol.span("&times;", class_="close")
modal += ol.img(class_="modal-content", id="img01")
modal += ol.div('', id="caption")
page.div(modal, id='myModal', class_="modal")

page.br()
page.h2("Contacts")
page.scripts({script: 'javascript'})
page.savehtml(hdir + 'index.html')

print "Elapsed time %d seconds" % int(time.time() - start)
