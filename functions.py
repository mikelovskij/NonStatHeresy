import h5py
import numpy
from scipy.signal import cheby1, lfilter, firwin
import inspect
from IPython.terminal.embed import InteractiveShellEmbed
import numpy as np
import virgotools as vrg
from time import time, ctime
import scipy.signal as sig


# ipshell debugger during run
def ipsh():
    ipshell = InteractiveShellEmbed(banner1="risolvitore di problemi attivato")

    frame = inspect.currentframe().f_back
    msg = 'Stopped at {0.f_code.co_filename} at line {0.f_lineno}'.format(frame)

    # Go back one level!
    # This is needed because the call to ipshell is inside the function ipsh()
    ipshell(msg, stack_depth=2)


def brms_reader(file, channelsandbands):
    with h5py.File(file, 'r') as f:
        brms_data = {}
        [gpsb, gpse] = f.attrs['gps']
        fsample = f.attrs['fsample']
        segments = f.attrs['segments']
        for (channel, bands) in channelsandbands.iteritems():
            brms_data[channel] = {}
            real_channel = ''
            split_channel = channel.split('_')
            for s in split_channel[0:-2]:
                real_channel += (s + '_')
            real_channel += split_channel[-2]
            times = f[real_channel]['times'][:]
            for b in bands:
                brms_data[channel][b] = f[real_channel][b][:]
    return gpsb, gpse, fsample, segments, times, brms_data


# come gestisco la segmentazione dei dati? Per come funziona al momento il calcolatore di psd, perdo l'effettiva corrispondenza temporale
# delle fft calcolate, ma forse potrebbe risultare non eccessivamente complesso il recuperarla. In questo modo le brsm avrebbero anche
# il loro asse temporale. Come gestisco quindi le eventuali discontinuita' di cotale asse? Se non son troppo brutte
# cioe' se i segmenti non sono troppo lunghi, probabilmente posso ignorarle, quando si calcolan correlazioni e coerenze


# problema 2: uso i trend, o segnali piu' veloci, tipo gli rms magari decimati a 10 Hz? Il problema e' che in ogni caso
# la risoluz in freq della psd e' 1/dt quindi se voglio brms a frequenza alta, voglio un dt basso e quindi perdo risoluz
# in frequenza.


# factor a number
def factors(n):
    return numpy.sort(reduce(list.__add__,
                             ([i, n // i] for i in range(1, int(n ** 0.5) + 1) if n % i == 0)))[1:-1]


# Custom decimation function, copied a long time ago from somewhere on the web
def decimate(x, q, n=None, ftype='iir', axis=-1):
    """downsample the signal x by an integer factor q, using an order n filter
    By default, an order 8 Chebyshev type I filter is used or a 30 point FIR
    filter with hamming window if ftype is 'fir'.
    (port to python of the GNU Octave function decimate.)
    Inputs:
        x -- the signal to be downsampled (N-dimensional array)
        q -- the downsampling factor
        n -- order of the filter (1 less than the length of the filter for a
             'fir' filter)
        ftype -- type of the filter; can be 'iir' or 'fir'
        axis -- the axis along which the filter should be applied
    Outputs:
        y -- the downsampled signal
    """

    if type(q) != type(1):
        print "q should be an integer"
        raise

    # check if the user is asking for too large decimation factor
    if q > 10:
        # compute factors
        qf = factors(q)
        if len(qf) != 0:
            # find the largest factor smaller than ten and decimate using it
            qnew = int(next(x for x in qf[::-1] if x <= 10))
            # decimate first using the cofactor (recursive call)
            x = decimate(x, q / qnew, n=n, ftype=ftype, axis=axis)
            # use what's left for the next step
            q = qnew

    if n is None:
        if ftype == 'fir':
            n = 30
        else:
            n = 4
    if ftype == 'fir':
        b = firwin(n + 1, 1. / q, window='hamming')
        y = lfilter(b, 1., x, axis=axis)
    else:
        (b, a) = cheby1(n, 0.05, 0.8 / q)

        y = lfilter(b, a, x, axis=axis)
    return y.swapaxes(0, axis)[::q].swapaxes(0, axis)


# TODO: add overlap support
# TODO: add correlation computation support
def simple_cumulative_psd_csd_correlation(main_data, main_times, auxchannels, n_points, segments, downfreq, data_source,
                                          nbins):
    s1 = np.zeros((int(n_points / 2) + 1), dtype='float64')
    freqs = np.zeros((int(n_points / 2) + 1), dtype='float64')
    s2 = np.zeros((len(auxchannels), int(n_points / 2) + 1), dtype='float64')
    csd = np.zeros((len(auxchannels), int(n_points / 2) + 1), dtype='complex128')
    sum_1 = 0
    square_sum_1 = 0
    aux_sum = np.zeros(len(auxchannels), dtype='float64')
    prod_sum = np.zeros(len(auxchannels), dtype='float64')
    square_aux_sum = np.zeros(len(auxchannels), dtype='float64')
    t0 = time()
    nfft = 0
    nsegments = len(segments)
    hist = []
    data_storage = {}  # some dictionary?
    for ((gpsb, gpse), j) in zip(segments, xrange(nsegments)):
        print 'Data acquisition and fft computation in progress, step {0} of {1} ...'.format(j, nsegments)
        fftperacquisition = int(np.floor((gpse - gpsb) * downfreq / n_points))
        caux = []
        start = np.argmin(abs(main_times - gpsb))
        end = np.argmin(abs(main_times - gpse))
        c = main_data[start:end]
        for ch in auxchannels:
            attempts = 0
            while attempts < 5:
                try:
                    with vrg.getChannel(data_source, ch, gpsb, gpse - gpsb) as ch2:
                        caux.append(decimator_wrapper(downfreq, ch2))
                    break
                except IOError, e:
                    attempts += 1
                    print e
        for i in xrange(fftperacquisition):
            freqs, s1temp = sig.welch(c[i * n_points: (i + 1) * n_points], fs=downfreq,
                                      window='hanning', nperseg=n_points)
            s1 += s1temp
            sum_1 += np.sum(c[i * n_points: (i + 1) * n_points])
            square_sum_1 += np.sum(c[i * n_points: (i + 1) * n_points] ** 2)
            for (c2, k) in zip(caux, xrange(len(caux))):
                _, csdtemp = sig.csd(c[i * n_points: (i + 1) * n_points],
                                     c2[i * n_points: (i + 1) * n_points], fs=downfreq,
                                     window='hanning', nperseg=n_points)
                _, s2temp = sig.welch(c2[i * n_points: (i + 1) * n_points], fs=downfreq,
                                      window='hanning', nperseg=n_points)
                s2[k] += s2temp
                csd[k] += csdtemp
                aux_sum[k] += np.sum(c2[i * n_points: (i + 1) * n_points])
                square_aux_sum[k] += np.sum(c2[i * n_points: (i + 1) * n_points] ** 2)
                prod_sum[k] += np.sum(c2[i * n_points: (i + 1) * n_points] * c[i * n_points: (i + 1) * n_points])
                # TODO: maybe calculating the histogram for every band and for every aux channel is a bit too expensive?
                if nfft == 0:
                    h, x_edges, y_edges = np.histogram2d(c[i * n_points: (i + 1) * n_points], c2[i * n_points: (i + 1) * n_points],
                                                         bins=nbins)
                    hist.append([h, x_edges, y_edges])
                else:
                    h, x_edges, y_edges = np.histogram2d(c[i * n_points: (i + 1) * n_points], c2[i * n_points: (i + 1) * n_points],
                                                         bins=[hist[k][1], hist[k][2]])
                    hist[k][0] += h

            nfft += 1
        t1 = time()
        tend = (t1 - t0) / (j + 1) * (nsegments - j)
        print 'Estimated completion {0}, in {1:.2f} minutes'.format(ctime(t1 + tend), (tend / 60))
    return s1 / nfft, s2 / nfft, csd / nfft, sum_1 / (nfft * n_points), square_sum_1 / (nfft * n_points), aux_sum / (nfft * n_points), square_aux_sum / (nfft * n_points), prod_sum / (nfft * n_points), freqs, hist


def less_simple_cumulative_psd_csd_correlation(main_data, main_times, auxchannels, n_points, segments, downfreq, data_source,
                                          nbins, data_storage=None):
    s1 = np.zeros((int(n_points / 2) + 1), dtype='float64')
    freqs = np.zeros((int(n_points / 2) + 1), dtype='float64')
    s2 = np.zeros((len(auxchannels), int(n_points / 2) + 1), dtype='float64')
    csd = np.zeros((len(auxchannels), int(n_points / 2) + 1), dtype='complex128')
    sum_1 = 0
    square_sum_1 = 0
    aux_sum = np.zeros(len(auxchannels), dtype='float64')
    prod_sum = np.zeros(len(auxchannels), dtype='float64')
    square_aux_sum = np.zeros(len(auxchannels), dtype='float64')
    t0 = time()
    nfft = 0
    nsegments = len(segments)
    hist = []
    acquire = False
    if not data_storage:
        acquire = True
        data_storage = {}
        for ch in auxchannels:
            data_storage[ch] = []
    for ((gpsb, gpse), j) in zip(segments, xrange(nsegments)):
        print 'Data acquisition and fft computation in progress, step {0} of {1} ...'.format(j, nsegments)
        fftperacquisition = int(np.floor((gpse - gpsb) * downfreq / n_points))
        caux = []
        start = np.argmin(abs(main_times - gpsb))
        end = np.argmin(abs(main_times - gpse))
        c = main_data[start:end]
        for ch in auxchannels:
            attempts = 0
            if acquire:
                while attempts < 5:
                    try:
                        with vrg.getChannel(data_source, ch, gpsb, gpse - gpsb) as ch2:
                            c_data = decimator_wrapper(downfreq, ch2)
                            data_storage[ch].append(c_data) # todo: forse qui non fa una copia e quindi quando vien modificato cambia?
                            caux.append(c_data)  # todo : no, in teoria append fa una copia.
                        break
                    except IOError, e:
                        attempts += 1
                        print e
                        if attempts == 5:
                            raise
            else:
                caux.append(data_storage[ch][j])
        for i in xrange(fftperacquisition):
            freqs, s1temp = sig.welch(c[i * n_points: (i + 1) * n_points], fs=downfreq,
                                      window='hanning', nperseg=n_points)
            s1 += s1temp
            sum_1 += np.sum(c[i * n_points: (i + 1) * n_points])
            square_sum_1 += np.sum(c[i * n_points: (i + 1) * n_points] ** 2)
            for (c2, k) in zip(caux, xrange(len(caux))):
                _, csdtemp = sig.csd(c[i * n_points: (i + 1) * n_points],
                                     c2[i * n_points: (i + 1) * n_points], fs=downfreq,
                                     window='hanning', nperseg=n_points)
                _, s2temp = sig.welch(c2[i * n_points: (i + 1) * n_points], fs=downfreq,
                                      window='hanning', nperseg=n_points)
                s2[k] += s2temp
                csd[k] += csdtemp
                aux_sum[k] += np.sum(c2[i * n_points: (i + 1) * n_points])
                square_aux_sum[k] += np.sum(c2[i * n_points: (i + 1) * n_points] ** 2)
                prod_sum[k] += np.sum(c2[i * n_points: (i + 1) * n_points] * c[i * n_points: (i + 1) * n_points])
                # TODO: maybe calculating the histogram for every band and for every aux channel is a bit too expensive?
                if nfft == 0:
                    h, x_edges, y_edges = np.histogram2d(c[i * n_points: (i + 1) * n_points], c2[i * n_points: (i + 1) * n_points],
                                                         bins=nbins)
                    hist.append([h, x_edges, y_edges])
                else:
                    h, x_edges, y_edges = np.histogram2d(c[i * n_points: (i + 1) * n_points], c2[i * n_points: (i + 1) * n_points],
                                                         bins=[hist[k][1], hist[k][2]])
                    hist[k][0] += h

            nfft += 1
        t1 = time()
        tend = (t1 - t0) / (j + 1) * (nsegments - j)
        print 'Estimated completion {0}, in {1:.2f} minutes'.format(ctime(t1 + tend), (tend / 60))
    return s1 / nfft, s2 / nfft, csd / nfft, sum_1 / (nfft * n_points), square_sum_1 / (nfft * n_points), aux_sum / (nfft * n_points), square_aux_sum / (nfft * n_points), prod_sum / (nfft * n_points), freqs, hist, data_storage


def even_less_simple_cumulative_psd_csd_correlation(main_data, main_times, auxchannels, n_points, segments, downfreq, data_source,
                                          nbins, data_storage=None):
    s1 = np.zeros((int(n_points / 2) + 1), dtype='float64')
    freqs = np.zeros((int(n_points / 2) + 1), dtype='float64')
    s2 = np.zeros((len(auxchannels), int(n_points / 2) + 1), dtype='float64')
    csd = np.zeros((len(auxchannels), int(n_points / 2) + 1), dtype='complex128')
    sum_1 = 0
    square_sum_1 = 0
    aux_sum = np.zeros(len(auxchannels), dtype='float64')
    prod_sum = np.zeros(len(auxchannels), dtype='float64')
    square_aux_sum = np.zeros(len(auxchannels), dtype='float64')
    t0 = time()
    nfft = 0
    nsegments = len(segments)
    hist = []
    acquire = False
    if not data_storage:
        acquire = True
        data_storage = {}
        for ch in auxchannels:
            data_storage[ch] = []
    for ((gpsb, gpse), j) in zip(segments, xrange(nsegments)):
        print 'Data acquisition and fft computation in progress, step {0} of {1} ...'.format(j, nsegments)
        fftperacquisition = int(np.floor((gpse - gpsb) * downfreq / n_points))
        start = np.argmin(abs(main_times - gpsb))
        end = np.argmin(abs(main_times - gpse))
        c = main_data[start:end]
        for ch, k in zip(auxchannels, xrange(len(auxchannels))):
            attempts = 0
            if acquire:
                while attempts < 5:
                    try:
                        with vrg.getChannel(data_source, ch, gpsb, gpse - gpsb) as ch2:
                            c_data = decimator_wrapper(downfreq, ch2)
                            data_storage[ch].append(c_data)
                            c2 = c_data
                        break
                    except IOError, e:
                        attempts += 1
                        print e
                        if attempts == 5:
                            raise
            else:
                c2 = data_storage[ch][j]
            for i in xrange(fftperacquisition):
                if k == 0:
                    freqs, s1temp = sig.welch(c[i * n_points: (i + 1) * n_points], fs=downfreq,
                                          window='hanning', nperseg=n_points)
                    s1 += s1temp
                    sum_1 += np.sum(c[i * n_points: (i + 1) * n_points])
                    square_sum_1 += np.sum(c[i * n_points: (i + 1) * n_points] ** 2)
                #for (c2, k) in zip(caux, xrange(len(caux))):
                _, csdtemp = sig.csd(c[i * n_points: (i + 1) * n_points],
                                     c2[i * n_points: (i + 1) * n_points], fs=downfreq,
                                     window='hanning', nperseg=n_points)
                _, s2temp = sig.welch(c2[i * n_points: (i + 1) * n_points], fs=downfreq,
                                      window='hanning', nperseg=n_points)
                s2[k] += s2temp
                csd[k] += csdtemp
                aux_sum[k] += np.sum(c2[i * n_points: (i + 1) * n_points])
                square_aux_sum[k] += np.sum(c2[i * n_points: (i + 1) * n_points] ** 2)
                prod_sum[k] += np.sum(c2[i * n_points: (i + 1) * n_points] * c[i * n_points: (i + 1) * n_points])
                # TODO: maybe calculating the histogram for every band and for every aux channel is a bit too expensive?
                if (i == 0) and (j == 0):
                    h, x_edges, y_edges = np.histogram2d(c[i * n_points: (i + 1) * n_points], c2[i * n_points: (i + 1) * n_points],
                                                         bins=nbins)
                    hist.append([h, x_edges, y_edges])
                else:
                    h, x_edges, y_edges = np.histogram2d(c[i * n_points: (i + 1) * n_points], c2[i * n_points: (i + 1) * n_points],
                                                         bins=[hist[k][1], hist[k][2]])
                    hist[k][0] += h
                if k == 0:
                    nfft += 1
        t1 = time()
        tend = (t1 - t0) / (j + 1) * (nsegments - j)
        print 'Estimated completion {0}, in {1:.2f} minutes'.format(ctime(t1 + tend), (tend / 60))
    return s1 / nfft, s2 / nfft, csd / nfft, sum_1 / (nfft * n_points), square_sum_1 / (nfft * n_points), aux_sum / (nfft * n_points), square_aux_sum / (nfft * n_points), prod_sum / (nfft * n_points), freqs, hist, data_storage


def decimator_wrapper(downfreq, ch):
    if ch.fsample < downfreq:
        print "Channels must have a larger or equal sampling frequency than the desired downsampled freq"
        raise ValueError
    else:
        if ch.fsample % downfreq != 0:
            print "Sampling frequency must be equal or an integer multiple of the downsampled frequency"
            raise ValueError
        else:
            if ch.fsample > downfreq:
                decimfactor = ch.fsample / downfreq
                #print "Decimation factor {0}".format(decimfactor)
                c = decimate(ch.data, int(decimfactor))
            else:
                c = ch.data
    return c


def get_channel_list(gpsb, source):
    channels = []
    with vrg.FrameFile(source) as ffl:
        with ffl.get_frame(gpsb) as frame:
            for adc in frame.iter_adc():
                if int(adc.contents.sampleRate) == 50:
                    channels.append(str(adc.contents.name))
    return np.array(channels)
