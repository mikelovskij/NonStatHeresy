#!/usr/bin/env python
# @author Michele Valentini
# /// Report generator for NonStatMoni BRMS monitor

import numpy as np
from argparse import ArgumentParser
from functions import brms_reader, ipsh, extractbands, Parameters,\
    string_repeater
import os
import cPickle
from data_processing import DataProcessing
from pathos import pools
from time import time, ctime
import plot_generation
import report_generator

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# **************************** Initialization ******************************* #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
start = time()

# read configuration options and process
parser = ArgumentParser()
parser.add_argument("-i", "--init", dest="initialization",
                    default='NonStatMoni.ini',
                    help="set initialization FILE", metavar="cfgFILE")
parser.add_argument("-b", "--brms", dest='brms_file',
                    help="set brms data file", metavar="brms_file")
parser.add_argument("-d", "--dir", dest="hdir",
                    help="output directory", metavar="OutDir")
parser.add_argument("-p", "--proc", dest="n_proc", metavar="NumberOfProcesses",
                    help='number of parallel processes to run', default=4,
                    type=int)
parser.add_argument("-n", "--ntop", dest="ntop",
                    default='10', type=int,
                    help="number of top coherences and correlation "
                         "cohefficients to be shown in the result table",
                    metavar="nTop")

args = vars(parser.parse_args())
if args['hdir']:
    hdir = args['hdir']
else:
    # set up a default result directory
    basedir = os.path.dirname(os.path.abspath(__file__))
    hdir = basedir + '/Results/{}/'.format(args['brms_file'].split('/')[-2])
    print "No result directory provided, results will be saved in" + hdir
# relative path (to hdir) of the plots folder
pdir = 'plots/'

# Create the results and plot directory if it does not already exist
try:
    os.mkdir(hdir)
except OSError:
    print 'Results Directory already exists'
try:
    os.mkdir(hdir + pdir)
except OSError:
    print 'Plots Directory already exists'

# set up and read the configuration file
par = Parameters(args['initialization'])

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# ++++++++++++++ Prepare the BRMS data for the post-processing ++++++++++++++ #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# Rebuild the band list as a list
for _, g_dict in par.group_dict.iteritems():
    extractbands(g_dict)
# Read the BRMS data stored in the hdf5 file and put it in the group_dict
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


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# +++++++++++++++++++++++ Do the actual post-processing +++++++++++++++++++++ #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# Allocate the post processing class
proc = DataProcessing(segments, fs, par.nav, par.group_dict, times)

# Process the brms bands psd and mean and square mean.
proc.cumulative_psd_computation()

# Create parallel pool and use it to compute the auxillary channels data.
pool = pools.ProcessPool(args['n_proc'])
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

# print "Post-processing complete, saving the results in {}".format(???)
# todo: un banale dump di tutta la classe?
# todo: oppure salvo come testo le variabili principali?
# todo: e magari in un altro testo i parametri usati?

# Build a smaller version of aux_results in order to use less space when saving
save_results = {}
for aux_name, result in aux_results.iteritems():
    save_results[aux_name] = {'brms_bands': result['brms_bands'],
                              'aux_psd': result['aux_psd'],
                              'histogram': result['histogram'],
                              'aux_name': result['aux_name']}

savingdict = {'aux_dict': par.aux_dict,
              'group_dict': par.group_dict,
              'times': times,
              'freqs': proc.freqs,
              'ccfs': proc.ccfs,
              'mean_cohs': proc.mean_cohs,
              'brms_psd': proc.brms_psd,
              'cohs': proc.cohs,
              'aux_results': save_results}
par.save_extended_config(**{'gpsb': gpsb,
                          'gpse': gpse,
                          'hdir': hdir,
                          'pdir': pdir,
                          'brms_file': args['brms_file']})

with open(hdir + 'post_proc_data.dat', mode='w') as f:
    cPickle.dump(savingdict, f)


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++++++++++++++++++++++ Generate the report +++++++++++++++++++++++++++++++ #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
print "Elapsed time %d seconds" % int(time() - start)
# Generate the plots, one group at a time.
print "starting the plot generation."
plot_generation.main(hdir + 'config.ini')
print 'starting the report generation'
report_generator.main(hdir + 'config.ini', args['ntop'])
print "Done!, results are located in" + hdir

print "Elapsed time %d seconds" % int(time() - start)
