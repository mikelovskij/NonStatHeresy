from virgotools import gps2str
import os
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import gridspec as grsp


# BRMS time series plots
def brms_time_plots(times, g_dict, savepath):
    # check the length of the plot and scale the time unit of measurement
    if max(times) <= 300:
        tunits = 's'
        t = times
    if (max(times) > 300) and (max(times) <= 60 * 100):
        t = times / 60
        tunits = 'min'
    if max(times) > (60 * 100):
        t = times / 3600
        tunits = 'h'

    # Plot the stuff
    _ = plt.figure(figsize=(10, 6))
    plt.axes([0.1, 0.1, 0.65, 0.8])
    for b_name, data in g_dict['brms_data'].iteritems():
        # todo: shouldn't I filter the data before this "downsampling"?
        plt.semilogy(t[::50], data[::50], linewidth=0.5, label=b_name)
    plt.grid(True)
    plt.axis(xmax=max(t))
    plt.xlabel('Time [' + tunits + ']')
    plt.ylabel('Normalized BRMS [' + g_dict['units'] + ']')
    plt.legend(loc=(1.02, 0.2))
    plt.title('{} BRMS {} - {}'.format(g_dict['channel'], gps2str(times[0]),
                                       gps2str(times[-1])), fontsize=10)
    plt.savefig(savepath)
    plt.close()


# Plot  band spectra
def brms_asd_plot(brms_psds, freqs, g_dict, savepath, gpsb, gpse):
    # TODO: can't i get n_points from the length of the vectors?
    plt.figure(figsize=(10, 6))
    plt.axes([0.1, 0.1, 0.65, 0.8])
    for band, sp in brms_psds.iteritems():
        try:
            plt.loglog(freqs, np.sqrt(sp), linewidth=0.5, label=band)
        except ValueError:
            print "Issue in plotting psd of {0}_{1}," \
                  " may there be no positive data?".format(g_dict['channel'],
                                                           band)
    plt.axis(xmax=freqs[-1])
    plt.axis(xmin=freqs[0])
    plt.grid(True)
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Spectrum of normalized BRMS [' + g_dict['units'] + '/rHz]')
    plt.legend(loc=(1.02, 0.2))
    plt.title('{} BRMS \n {} - {}'.format(g_dict['channel'], gps2str(gpsb),
                                          gps2str(gpse)), fontsize=10)
    plt.savefig(savepath)
    plt.close()


def auxiliary_plots(group, aux_dict, g_dict, freqs, ccfs, cohs, mean_cohs,
                    aux_results, gpsb, gpse, plt_path,
                    ccf_thresh=0.15, mn_coh_thresh=0.05):
    for aux_name, aux_groups in aux_dict.iteritems():
        if group in aux_groups:
            # hicorrlist = []
            hi_corr = []
            # hicohlist = []
            hi_coh = []
            # plotting every coherence is probably too much.
            #  fo una cosa tipo coherence?
            for j, b_name in enumerate(g_dict['brms_data'].iterkeys()):
                ccf = ccfs[aux_name][group + '_' + b_name]
                mean_coh = mean_cohs[aux_name][group + '_' + b_name]
                # TODO: surely some distributions can give a significance estimate of the corrcoef
                if abs(ccf) >= ccf_thresh:
                    hi_corr.append((b_name, j))
                if np.mean(mean_coh) >= mn_coh_thresh:
                    hi_coh.append((b_name, j))
            plt_path += group + '_cohe_' + aux_name.split(':')[1] + '.png'
            if (len(hi_corr) + len(hi_coh)) > 0:
                plt.figure(figsize=(15, 6))
                if len(hi_corr) > 0:
                    n_rows = int(np.floor(np.sqrt(len(hi_corr))))
                    n_cols = int(np.ceil(np.sqrt(len(hi_corr))))
                    gs = grsp.GridSpec(n_cols, 2 * n_rows)
                    gs.update(wspace=0.2, hspace=0.5)
                    for j in xrange(n_cols):
                        for k, (b_name, b_num) in zip(xrange(n_rows),
                                                      hi_corr[n_rows * j:
                                                              n_rows *
                                                              (j + 1)]):
                            h = aux_results[aux_name]['histogram'][b_num]
                            ax = plt.subplot(gs[j, n_rows + k])
                            ax.pcolormesh(h[1], h[2], h[0])
                            ax.tick_params(axis='both', which='minor',
                                           labelsize=6)
                            ax.tick_params(axis='both', which='major',
                                           labelsize=6)
                            ax.ticklabel_format(style="sci", scilimits=(0, 0),
                                                axis="both")
                            plt.title(b_name)
                    ax1 = plt.subplot(gs[0:n_cols, 0:n_rows])

                else:
                    ax1 = plt.subplot2grid((1, 1), (0, 0))
                for (b_name, b_num) in hi_coh:
                    try:
                        coh = cohs[aux_name][group + '_' + b_name]
                        ax1.semilogx(freqs, coh, linewidth=0.5,
                                     label=b_name)
                        ax1.axis(ymin=0, ymax=1)
                        ax1.axis(xmin=freqs[0])
                        ax1.grid(True)
                        ax1.set_xlabel('Frequency [Hz]')
                        ax1.set_ylabel('Coherence')
                        ax1.legend(loc=(1.02, 0.2))
                    except ValueError:
                        print "Issue in plotting coherence of  {0} _ {1}" \
                              ", could there be no positive data?" \
                            .format(g_dict['channel'], b_name)
                ax1.set_title('Coherence with {}, \n {:d} - {:d}'.
                              format(aux_name, gps2str(gpsb), gps2str(gpse)))
                plt.savefig(plt_path)
                plt.close()
            else:
                try:
                    os.remove(plt_path)
                except OSError:
                    pass
