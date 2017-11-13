import numpy as np
import scipy.signal as sig
import virgotools as vrg
from functions import decimate, tryfivetimes, ipsh # TODO: aggiornare un po' sto decimate dai
import scipy.sparse as sp

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
        nfft = 0
        hist = []
        for ((gpsb, gpse), j) in zip(self.segments, xrange(self.nsegments)):
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
                        h, x_edges, y_edges = self.sparse_histogram(
                            band_data,
                            aux_data[i * self.n_points:(i + 1) * self.n_points]
                            , n_bins=self.nbins)
                        hist.append([h, x_edges, y_edges])
                    else:
                        h, x_edges, y_edges = self.update_histogram(hist[k][0],
                            hist[k][1], hist[k][2],
                            band_data,
                            aux_data[i * self.n_points:(i + 1) * self.n_points]
                            )
                nfft += 1
        aux_psd /= nfft
        abs_csds = np.absolute(csds / nfft) ** 2
        aux_sum /= (nfft * self.n_points)
        aux_square_sum /= (nfft * self.n_points)
        prod_sum /= (nfft * self.n_points)
        return{'brms_bands': brms_bands,
               'aux_psd': aux_psd,
               'aux_mean': aux_sum,
               'aux_square_mean': aux_square_sum,
               'prod_mean': prod_sum,
               'abs_csds': abs_csds,
               'histogram': hist,
               'aux_name': aux_channel}

    def coherence_computation(self, aux_results):
        for aux_name, aux_dict in aux_results.iteritems():
            self.cohs[aux_name] = {}
            self.mean_cohs[aux_name] = {}
            brms_dict = aux_dict['brms_bands']
            aux_psd = aux_dict['aux_psd']
            for (group, band), k in zip(brms_dict, xrange(len(brms_dict))):
                abs_csd = aux_dict['abs_csds'][k]
                band_psd = self.brms_psd[group][band]
                self.cohs[aux_name][group + '_' + band] = abs_csd / (band_psd *
                                                                     aux_psd)
                self.mean_cohs[aux_name][group + '_' + band] = np.mean(self.cohs[aux_name][group + '_' + band])

    def pearson_cohefficient_computation(self, aux_results):
        for aux_name, aux_dict in aux_results.iteritems():
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

    @staticmethod
    def sparse_histogram(x, y, n_bins):
        h, x_edges, y_edges = np.histogram2d(x, y, bins=n_bins)
        return sp.bsr_matrix(np.int16(h)), x_edges, y_edges

    @staticmethod
    def update_histogram(hist, x_edges, y_edges, x, y, border_fraction=0.1):
        # check  the data borders and enlarge them if necessary
        enlarged_edges = []
        enlargement_size = []
        for edges, data in zip([x_edges, y_edges], [x, y]):
            bin_size = edges[1] - edges[0]
            added_size_min = 0
            added_size_max = 0
            interval = (max(data) - min(data)) * border_fraction
            if max(data) > edges[-1]:
                added_edges = np.arange(edges[-1] + bin_size,
                                        max(data) + interval,
                                        bin_size)
                added_size_max = len(added_edges)
                new_edges = np.concatenate((edges, added_edges))
            else:
                new_edges = edges
            if min(data) < new_edges[0]:
                added_edges = np.arange(new_edges[0] - bin_size,
                                        min(data) - interval,
                                        -bin_size)[::-1]
                added_size_min = len(added_edges)
                new_edges = np.concatenate((added_edges, new_edges))
            enlarged_edges.append(new_edges)
            enlargement_size.append([added_size_min, added_size_max])
        h, x_edges, y_edges = np.histogram2d(x, y, bins=enlarged_edges)

        # enlarge x size of the old histogram by zero padding
        if enlargement_size[0][0]:
            top = sp.bsr_matrix(np.zeros([enlargement_size[0][0],
                                          hist.shape[1]]))
            hist = sp.vstack([top, hist])
        if enlargement_size[0][1]:
            bot = sp.bsr_matrix(np.zeros([enlargement_size[0][1],
                                          hist.shape[1]]))
            hist = sp.vstack([hist, bot])

        # enlarge y size of the old histogram by zero padding
        if enlargement_size[1][0]:
            left = sp.bsr_matrix(np.zeros([hist.shape[0],
                                           enlargement_size[1][0]]))
            hist = sp.hstack([left, hist])
        if enlargement_size[1][1]:
            right = sp.bsr_matrix(np.zeros([hist.shape[0],
                                            enlargement_size[1][1]]))
            hist = sp.hstack([hist, right])

        return sp.bsr_matrix(np.int16(h)) + hist, x_edges, y_edges


    @tryfivetimes
    def get_channel_data(self, data_source, channel, gpsb, gpse):
        # TODO: extract other stuff from ch i.e. the units?
        with vrg.getChannel(data_source, channel, gpsb, gpse - gpsb) as ch:
            data = self.decimator_wrapper(ch)
        return data

