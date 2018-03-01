import numpy as np
import scipy.signal as sig
import virgotools as vrg
from functions import retry, ipsh, decimator_wrapper # TODO: aggiornare un po' sto decimate dai
import scipy.sparse as sp


# Class for the post_processing of the brms data with auxillary channels
class DataProcessing:
    def __init__(self, segments, down_freq, n_averages, group_dic, times,
                 outliers_frac, overlap):
        self.segments = segments
        self.nsegments = len(self.segments)
        self.down_freq = down_freq
        self.n_averages = n_averages
        self.group_dic = group_dic
        self.brms_psd = {}
        self.brms_mean = {}
        self.overlap = overlap
        self.outliers_quantile = 1 - (1 - outliers_frac) / 2
        self.brms_sqmean = {}
        self.times = times
        self.freqs = []
        self.cohs = {}
        self.ccfs = {}
        self.mean_cohs = {}
        # Rice estimator for nbins
        estimated_len = 0
        for seg in self.segments:
            estimated_len += (seg[1] - seg[0])
        estimated_points = int(np.floor(float(estimated_len * self.down_freq) * outliers_frac))
        self.nbins = int(np.ceil((estimated_points ** (1.0 / 3))))
        # estimate the number of points for each fft
        est_points_w_overlap = float(estimated_points / (1 - self.overlap))
        self.n_points = int(pow(2, np.floor(np.log2(est_points_w_overlap /
                                                    self.n_averages))))
        # shift between each FFT window
        self.n_over_points = int(np.rint(self.n_points * (1 - self.overlap)))

    # returns the indexes corresponding to the gpse and gpsb times specified
    # and the number of fft windows that can be fitted in the segment
    def segment_indexes(self, gpsb, gpse):
        # todo: the fft_per_segments computed here is no longer used, remove it
        fft_per_segment = int(np.floor((gpse - gpsb) *
                                       self.down_freq / self.n_points))
        gpsb_index = int(np.argmin(abs(self.times - gpsb)))
        gpse_index = int(np.argmin(abs(self.times - gpse)))
        return fft_per_segment, gpsb_index, gpse_index


    # computes the psds for the various brmss and stores them in a dictionary
    def cumulative_psd_computation(self):
        # todo: implement outlier removal in this computation too
        # but then I should remove the brms outliers from the cross computation too
        # maybe it is better to clean the brmss when they are generated
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
                total_points = 0
                for (gpsb, gpse), j in zip(self.segments,
                                           xrange(self.nsegments)):
                    _, start, end = self.segment_indexes(gpsb, gpse)
                    #print "fftperseg"
                    #print fft_per_segm
                    #print "npoints"
                    #print self.n_points
                    #print "end-start"
                    #print end - start
                    seg_points = end - start
                    total_points += seg_points
                    fft_per_segm = (seg_points - self.n_points) / self.n_over_points + 1
                    for i in xrange(fft_per_segm):
                        s_start = start + i * self.n_over_points
                        s_end = s_start + self.n_points
                        self.freqs, spec = sig.welch(
                            brms[s_start: s_end],
                            fs=self.down_freq,
                            window='hanning', nperseg=self.n_points)
                        self.brms_psd[group][band] += spec
                        n_fft += 1
                    self.brms_mean[group][band] += np.sum(
                            brms[start:end])
                    self.brms_sqmean[group][band] += np.sum(
                            brms[start:end] ** 2)

                self.brms_psd[group][band] /= n_fft
                self.brms_mean[group][band] /= total_points
                self.brms_sqmean[group][band] /= total_points

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
        total_points = 0
        first = True
        hist = []
        for ((gpsb, gpse), j) in zip(self.segments, xrange(self.nsegments)):
            aux_data = self.get_channel_data(data_source,
                                             aux_channel, gpsb, gpse)
            min_thr, max_thr = np.percentile(aux_data,
                                             [100 * (1 - self.outliers_quantile),
                                              100 * self.outliers_quantile])
            bad_idxs = np.nonzero(np.logical_or((aux_data < min_thr),
                                              (aux_data > max_thr)))[0]
            total_points += len(aux_data) - len(bad_idxs)
            id_segments = []
            previous_bad = -1
            for bad in bad_idxs.sort():
                id_segments.append([previous_bad + 1, bad])
                previous_bad = bad
            id_segments.append([previous_bad + 1, len(aux_data)])
            _, start, end = self.segment_indexes(gpsb, gpse)
            for seg in id_segments:
                segment_points = int(seg[1] - seg[0])
                fft_per_segm = (segment_points - self.n_points) / self.n_over_points + 1
                for i in xrange(fft_per_segm):
                    s_start = seg[0] + i * self.n_over_points
                    s_end = s_start + self.n_points
                    _, psdtemp = sig.welch(aux_data[start:end],
                                           fs=self.down_freq,
                                           window='hanning',
                                           nperseg=self.n_points)
                    aux_psd += psdtemp

                    for (group, band), k in zip(brms_bands,
                                                xrange(len(brms_bands))):
                        band_data = self.group_dic[group]['brms_data'][band][
                                    start + s_start: start + s_end]
                        _, csdtemp = sig.csd(band_data,
                                             aux_data[
                                             i * self.n_points: (i + 1) *
                                                                self.n_points],
                                             fs=self.down_freq,
                                             window='hanning',
                                             nperseg=self.n_points)
                        csds[k] += csdtemp
                    nfft += 1
            good_mask = np.ones(len(aux_data), dtype=bool)
            good_mask[bad_idxs] = False

            aux_sum += np.sum(aux_data[good_mask])
            aux_square_sum += np.sum(aux_data[good_mask]**2)
            for (group, band), k in zip(brms_bands,
                                        xrange(len(brms_bands))):
                band_data = self.group_dic[group]['brms_data'][band][
                            start:end]
                ipsh()
                # todo: controllo che sia della lunghezza ggiusta
                prod_sum[k] += np.sum(aux_data[good_mask] * band_data[good_mask])
                # TODO: since they are linear, I could save only the edges of the edges
                if first:
                    h, x_edges, y_edges = self.sparse_histogram(
                        np.log(band_data[good_mask]), aux_data[good_mask],
                        n_bins=self.nbins)
                    hist.append([h, x_edges, y_edges])
                else:
                    h, x_edges, y_edges = self.update_histogram(hist[k][0],
                        hist[k][1], hist[k][2],
                        np.log(band_data[good_mask]),
                        aux_data[good_mask])

                    hist[k] = [h, x_edges, y_edges]

        aux_psd /= nfft
        abs_csds = np.absolute(csds / nfft) ** 2
        aux_sum /= total_points
        aux_square_sum /= total_points
        prod_sum /= total_points
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
        x_width = x_edges[1] - x_edges[0]
        y_width = y_edges[1] - y_edges[0]
        # todo: check that linspace as used in update histogram gives the same results
        return (sp.coo_matrix(h.astype('int16')),
                [x_edges[0], x_edges[-1], x_width, len(x_edges)],
                [y_edges[0], y_edges[-1], y_width, len(y_edges)])

    @staticmethod
    def update_histogram(hist, x_lim, y_lim, x, y, border_fraction=0.05):
        # check  the data borders and enlarge them if necessary
        enlarged_edges = []
        enlargement_size = []
        for limits, data in zip([x_lim, y_lim], [x, y]):
            old_edges = np.linspace(limits[0], limits[1], limits[3])
            added_size_min = 0
            added_size_max = 0
            mx = np.max(data)
            mn = np.min(data)
            interval = (mx - mn) * border_fraction
            if mx > limits[1]:
                added_edges = np.arange(limits[1] + limits[2],
                                        mx + interval,
                                        limits[2])
                added_size_max = len(added_edges)
                new_edges = np.concatenate((old_edges, added_edges))
            else:
                new_edges = old_edges
            if mn < new_edges[0]:
                added_edges = np.arange(new_edges[0] - limits[2],
                                        mn - interval,
                                        -limits[2])[::-1]
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

        return (sp.coo_matrix(h.astype('int16')) + hist,
                [x_edges[0], x_edges[-1], x_lim[2], len(x_edges)],
                [y_edges[0], y_edges[-1], y_lim[2], len(y_edges)])

    @retry
    def get_channel_data(self, data_source, channel, gpsb, gpse):
        # TODO: extract other stuff from ch i.e. the units?
        with vrg.getChannel(data_source, channel, gpsb, gpse - gpsb) as ch:
            data = decimator_wrapper(self.down_freq, ch)
        return data

