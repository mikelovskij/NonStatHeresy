#!/usr/bin/env python
# @author Gabriele Vajente
# /// Report generator for NonStatMoni BRMS monitor


import numpy as np

from argparse import ArgumentParser
import markup
from markup import oneliner as ol
from virgotools import gps2str
import tables as tb
from functions import brms_reader, ipsh, extractbands, Parameters,\
    string_repeater
import os
from data_processing import DataProcessing
import matplotlib
from pathos import pools
from time import time, ctime
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from collections import OrderedDict
from plot_generation import brms_time_plots, brms_asd_plot, auxiliary_plots

start = time()
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
                    help='number of parallel processes to run', default=4,
                    type=int)

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

pool = pools.ProcessPool(args['n_proc'])
pool_args = []


aux_results = {}
t0 = time()
n_aux = len(par.aux_dict)
for j, res in enumerate(pool.uimap(proc.auxillary_psd_csd_correlation,
                        par.aux_dict.iterkeys(), par.aux_dict.itervalues(),
                        string_repeater(par.aux_source, n_aux))):
    aux_results[res['aux_name']] = res
    t1 = time()
    tend = (t1 - t0) / (j + 1) * (n_aux - j)
    print 'Estimated completion {0}, in {1:.2f} minutes'.format(
        ctime(t1 + tend), (tend / 60))
    print "completed channel {} of {}".format(j, n_aux)

proc.pearson_cohefficient_computation(aux_results)
proc.coherence_computation(aux_results)
print "Post-processing complete, saving the results in {}".format(???)
# todo: un banale dump di tutta la classe?
# todo: oppure salvo come testo le variabili principali?
# todo: e magari in un altro testo i parametri usati?
# quali sono le cose importanti da salvare:
  # time
  # freqs
  # cohs
  # ccfs
  # mean cohs
  # histograms
  # group dict? hmm non basta il cfg e il cfg e hdf5 reader?
  # segmenti coaletti?
  # main psd
  # forse e' meglio iniziare dall'"incapsulaggio " dei plot e della generaz pagine
  # saving paths

# Generate the plots, one group at a time.
print "starting the plot generation."
for group, g_dict in par.group_dict.iteritems():
    # Plot time series of the brms
    brms_time_plots(times, g_dict, hdir + pdir + group + '_time.png')
    # Plot ASDs of the brms
    brms_asd_plot(proc.brms_psd[group], proc.freqs, g_dict,
                  hdir + pdir + group + '_psd.png', gpsb, gpse)

    # Plot the auxilliary channels coherences and histograms if their
    # mean coherence or correlation cohefficient surpasses a threshold
    auxiliary_plots(group, par.aux_dict, g_dict, proc.freqs, proc.ccfs,
                    proc.cohs, proc.mean_cohs, aux_results, gpsb, gpse,
                    hdir + pdir, ccf_thresh=0.15, mn_coh_thresh=0.05)


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
            page2.img(src=(pdir + group + '_' + 'cohe_' +
                           aux_name.split(':')[1] + '.png'),
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

print "Elapsed time %d seconds" % int(time() - start)
