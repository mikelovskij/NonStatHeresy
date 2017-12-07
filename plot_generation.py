import cPickle
from config_manager import Parameters
from argparse import ArgumentParser
from plot_functions import brms_asd_plot, brms_time_plots, auxiliary_plots


class Data:
    def __init__(self, filename):
        with open(filename) as f:
            savingdict = cPickle.load(f)
        self.ccfs = savingdict['ccfs']
        self.mean_cohs = savingdict['mean_cohs']
        self.group_dict = savingdict['group_dict']
        self.times = savingdict['times']
        self.freqs = savingdict['freqs']
        self.aux_dict = savingdict['aux_dict']
        self.brms_psd = savingdict['brms_psd']
        self.cohs = savingdict['cohs']
        self.aux_results = savingdict['aux_results']


# Plot generating function
def main(initialization):
    par = Parameters(initialization)
    gpsb = par.res_param['gpsb']
    gpse = par.res_param['gpse']
    pdir = par.res_param['pdir']
    hdir = par.res_param['hdir']
    data = Data(hdir + 'post_proc_data.dat')
    for group, g_dict in data.group_dict.iteritems():
        # Plot time series of the brms
        brms_time_plots(data.times, g_dict, hdir + pdir + group + '_time.png')
        # Plot ASDs of the brms
        brms_asd_plot(data.brms_psd[group], data.freqs, g_dict,
                      hdir + pdir + group + '_psd.png', gpsb, gpse)

        # Plot the auxilliary channels coherences and histograms if their
        # mean coherence or correlation cohefficient surpasses a threshold
        auxiliary_plots(group, data.aux_dict, g_dict, data.freqs, data.ccfs,
                        data.cohs, data.mean_cohs, data.aux_results, gpsb,
                        gpse, hdir + pdir, ccf_thresh=0.5, mn_coh_thresh=0.05)

if __name__ == "__main__":
    # read configuration options and dataess
    parser = ArgumentParser()
    parser.add_argument("-i", "--init", dest="initialization",
                        default='NonStatMoni.ini',
                        help="set initialization FILE", metavar="cfgFILE")
    args = vars(parser.parse_args())
    # launch main function
    main(args['initialization'])
