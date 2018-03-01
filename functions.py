import h5py
import numpy
import inspect
from IPython.terminal.embed import InteractiveShellEmbed
import scipy.signal as sig
from collections import OrderedDict
import virgotools as vrg
from time import sleep


# ipshell debugger during run
def ipsh():
    ipshell = InteractiveShellEmbed(banner1="risolvitore di problemi attivato")

    frame = inspect.currentframe().f_back
    msg = 'Stopped at {0.f_code.co_filename} at line {0.f_lineno}'.format(
        frame)

    # Go back one level!
    # This is needed because the call to ipshell is inside the function ipsh()
    ipshell(msg, stack_depth=2)


def string_repeater(string, n):
    u = 0
    while u < n:
        yield string
        u += 1


def brms_reader(file_name, group_dict):
    with h5py.File(file_name, 'r') as f:
        [gpsb, gpse] = f.attrs['gps']
        sampling_f = f.attrs['fsample']
        segments = f.attrs['segments']
        for (group, g_dic) in group_dict.iteritems():
            g_dic['brms_data'] = OrderedDict()
            channel = g_dic['channel']
            times = f[channel]['times'][:]
            for band in g_dic['band_list']:
                g_dic['brms_data'][band] = f[channel][band][:]
    return gpsb, gpse, sampling_f, segments, times


# factor a number
def factors(n):
    return numpy.sort(reduce(list.__add__,
                             ([i, n // i] for i in range(1, int(n ** 0.5) + 1)
                              if n % i == 0)))[1:-1]


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
        raise ValueError

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
        b = sig.firwin(n + 1, 1. / q, window='hamming')
        y = sig.lfilter(b, 1., x, axis=axis)
    else:
        (b, a) = sig.cheby1(n, 0.05, 0.8 / q)

        y = sig.lfilter(b, a, x, axis=axis)
    return y.swapaxes(0, axis)[::q].swapaxes(0, axis)

# checks that the conditions required for the decimation are met
# and if true calculates the decimation factor required for the wanted
# downsampling_frequency and then calls the decimation function
# TODO: update these functions
def decimator_wrapper(ds_freq, ch):
    f_sample = round(ch.fsample, 5)
    ds_freq = round(ds_freq, 5)
    if f_sample < ds_freq:
        print "Channels must have a larger or equal sampling frequency than " \
              "the desired downsampled freq"
        raise ValueError
    else:
        if f_sample % ds_freq != 0:
            err  = ("Sampling frequency {:f} must be equal or an integer " \
                  "multiple of the downsampled frequency {:f} {:f}".format(f_sample, ds_freq, f_sample % ds_freq))
            raise ValueError(err)
        else:
            if f_sample > ds_freq:
                decimfactor = f_sample / ds_freq
                # print "Decimation factor {0}".format(decimfactor)
                c = decimate(ch.data, int(decimfactor))
            else:
                c = ch.data
    return c

def retry(f):
    def wrapper(that, *args):
        attempts = 0
        out = None
        while attempts < 10:
            try:
                out = f(that, *args)
                break
            except (IOError, vrg.frame_lib.ChannelNotFound) as e:
                attempts += 1
                print e
                sleep(10)
                if attempts == 10:
                    raise
        return out

    return wrapper
