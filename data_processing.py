import numpy as np
import scipy.signal as sig
from time import time, ctime
import virgotools as vrg
from functions import decimate, tryfivetimes, ipsh # TODO: aggiornare un po' sto decimate dai


# Class for the post_processing of the brms data with auxillary channels
class DataProcessing:
    def __init__(self, segments, down_freq, n_averages, group_dic, times):
        self.segments = segments
        self.nsegments = len(self.segments)
        self.down_freq = down_freq
        self.n_averages = n_averages
        self.group_dic = group_dic
        self.brms_psd = {}
        self.brms_mean = {}
        self.brms_sqmean = {}
        self.times = times
        self.freqs = []
        self.aux_results = {}
        self.cohs = {}
        self.ccfs = {}
        self.mean_cohs = {}
        # Rice estimator for nbins
        estimated_points = 0
        for seg in self.segments:
            estimated_points += (seg[1] - seg[0]) * self.down_freq
        self.nbins = int(np.ceil(2 * (estimated_points ** (1.0 / 3))))
        self.n_points = int(pow(2, np.floor(np.log((2 * estimated_points) / (self.n_averages + 1)) / np.log(2))))

    # returns the indexes corresponding to the gpse and gpsb times specified
    # and the number of fft windows that can be fitted in the segment
    def segment_indexes(self, gpsb, gpse):
        fft_per_segment = int(np.floor((gpse - gpsb) *
                                       self.down_freq / self.n_points))
        gpsb_index = np.argmin(abs(self.times - gpsb))
        gpse_index = np.argmin(abs(self.times - gpse))
        return fft_per_segment, gpsb_index, gpse_index

    # checks that the conditions required for the decimation are met
    # and if true calculates the decimation factor required for the wanted
    # downsampling_frequency and then calls the decimation function
    def decimator_wrapper(self, ch):
        if ch.fsample < self.down_freq:
            print "Channels must have a larger or equal sampling frequency " \
                  "than the desired downsampled freq"
            raise ValueError
        else:
            if ch.fsample % self.down_freq != 0:
                print "Sampling frequency must be equal or an integer " \
                      "multiple of the downsampled frequency"
                raise ValueError
            else:
                if ch.fsample > self.down_freq:
                    decimfactor = ch.fsample / self.down_freq
                    # print "Decimation factor {0}".format(decimfactor)
                    c = decimate(ch.data, int(decimfactor))
                else:
                    c = ch.data
        return c

    # computes the psds for the various brmss and stores them in a dictionary
    def cumulative_psd_computation(self):
        for group, dic in self.group_dic.iteritems():
            self.brms_psd[group] = {}
            self.brms_mean[group] = {}
            self.brms_sqmean[group] = {}
            for band, brms in dic['brms_data'].iteritems():
                self.brms_psd[group][band] = np.zeros(
                    (int(self.n_points / 2) + 1), dtype='float64')
                self.brms_mean[group][band] = 0
                self.brms_sqmean[group][band] = 0
                n_fft = 0
                for (gpsb, gpse), j in zip(self.segments,
                                           xrange(self.nsegments)):
                    fft_per_segm, start, end = self.segment_indexes(gpsb, gpse)
                    for i in xrange(fft_per_segm):
                        self.freqs, spec = sig.welch(
                            brms[start + i * self.n_points:
                                 start + (i + 1) * self.n_points],
                            fs=self.down_freq,
                            window='hanning', nperseg=self.n_points)
                        self.brms_psd[group][band] += spec
                        self.brms_mean[group][band] += np.sum(
                            brms[start + i * self.n_points:
                                 start + (i + 1) * self.n_points])
                        self.brms_sqmean[group][band] += np.sum(
                            brms[start + i * self.n_points:
                                 start + (i + 1) * self.n_points] ** 2)
                        n_fft += 1
                self.brms_psd[group][band] /= n_fft
                self.brms_mean[group][band] /= (n_fft * self.n_points)
                self.brms_sqmean[group][band] /= (n_fft * self.n_points)

    def auxillary_psd_csd_correlation(self, aux_channel, aux_groups,
                                      data_source):
        brms_bands = []
        for group in aux_groups:
            for band in self.group_dic[group]['brms_data'].iterkeys():
                brms_bands.append((group, band))
        aux_psd = np.zeros(int(self.n_points / 2 + 1), dtype='float64')
        csds = np.zeros((len(brms_bands), int(self.n_points / 2) + 1),
                        dtype='complex128')
        aux_sum = 0
        prod_sum = np.zeros(len(brms_bands), dtype='float64')
        aux_square_sum = 0
        t0 = time()
        nfft = 0
        hist = []
        for ((gpsb, gpse), j) in zip(self.segments, xrange(self.nsegments)):
            print 'Processing channel {0}, step {1} of {2} ...'.format(
                aux_channel, j, self.nsegments)
            fft_per_segm, start, end = self.segment_indexes(gpsb, gpse)
            aux_data = self.get_channel_data(data_source,
                                             aux_channel, gpsb, gpse)
            for i in xrange(fft_per_segm):
                _, psdtemp = sig.welch(aux_data[i * self.n_points
                                                : (i + 1) * self.n_points],
                                       fs=self.down_freq,
                                       window='hanning',
                                       nperseg=self.n_points)
                aux_psd += psdtemp
                aux_sum += np.sum(aux_data[i * self.n_points:
                                           (i + 1) * self.n_points])
                aux_square_sum += np.sum(aux_data[i * self.n_points
                                                  :(i + 1) * self.n_points]**2)
                for (group, band), k in zip(brms_bands,
                                            xrange(len(brms_bands))):
                    band_data = self.group_dic[group]['brms_data'][band][
                                start + i * self.n_points:
                                start + (i + 1) * self.n_points]

                    _, csdtemp = sig.csd(band_data,
                                         aux_data[i * self.n_points: (i + 1) *
                                                  self.n_points],
                                         fs=self.down_freq,
                                         window='hanning',
                                         nperseg=self.n_points)
                    csds[k] += csdtemp
                    prod_sum[k] += np.sum(
                        aux_data[i * self.n_points: (i + 1) * self.n_points] *
                        band_data)
                    if nfft == 0:
                        h, x_edges, y_edges = np.histogram2d(
                            band_data,
                            aux_data[i * self.n_points:(i + 1) * self.n_points]
                            , bins=self.nbins)
                        hist.append([h, x_edges, y_edges])
                    else:
                        h, x_edges, y_edges = np.histogram2d(
                            band_data,
                            aux_data[i * self.n_points:(i + 1) * self.n_points]
                            , bins=[hist[k][1], hist[k][2]])
                        hist[k][0] += h
                nfft += 1
            t1 = time()
            tend = (t1 - t0) / (j + 1) * (self.nsegments - j)
            print 'Estimated completion {0}, in {1:.2f} minutes'.format(
                ctime(t1 + tend), (tend / 60))
        aux_psd /= nfft
        csds /= nfft
        aux_sum /= (nfft * self.n_points)
        aux_square_sum /= (nfft * self.n_points)
        prod_sum /= (nfft * self.n_points)
        self.aux_results[aux_channel] = {'brms_bands': brms_bands,
                                         'aux_psd': aux_psd,
                                         'aux_mean': aux_sum,
                                         'aux_square_mean': aux_square_sum,
                                         'prod_mean': prod_sum,
                                         'csds': csds,
                                         'histogram': hist}

    def coherence_computation(self):
        for aux_name, aux_dict in self.aux_results.iteritems():
            self.cohs[aux_name] = {}
            self.mean_cohs[aux_name] = {}
            brms_dict = aux_dict['brms_bands']
            aux_psd = aux_dict['aux_psd']
            for (group, band), k in zip(brms_dict, xrange(len(brms_dict))):
                csd = aux_dict['csds'][k]
                band_psd = self.brms_psd[group][band]
                self.cohs[aux_name][group + '_' + band] = (np.absolute(csd) ** 2) \
                    / (band_psd * aux_psd)
                self.mean_cohs[aux_name][group + '_' + band] = np.mean(self.cohs[aux_name][group + '_' + band])

    def pearson_cohefficient_computation(self):
        for aux_name, aux_dict in self.aux_results.iteritems():
            self.ccfs[aux_name] = {}
            brms_dict = aux_dict['brms_bands']
            a_mn = aux_dict['aux_mean']
            a_sq_mn = aux_dict['aux_square_mean']
            for (group, band), k in zip(brms_dict, xrange(len(brms_dict))):
                p_mn = aux_dict['prod_mean'][k]
                b_mn = self.brms_mean[group][band]
                b_sq_mn = self.brms_sqmean[group][band]
                self.ccfs[aux_name][group + '_' + band] = (p_mn - b_mn * a_mn) \
                    / np.sqrt((b_sq_mn - b_mn ** 2) * (a_sq_mn - a_mn ** 2))

    @tryfivetimes
    def get_channel_data(self, data_source, channel, gpsb, gpse):
        # TODO: extract other stuff from ch i.e. the units?
        with vrg.getChannel(data_source, channel, gpsb, gpse - gpsb) as ch:
            data = self.decimator_wrapper(ch)
        return data

