import h5py
import numpy
import inspect
from IPython.terminal.embed import InteractiveShellEmbed
import numpy as np
import virgotools as vrg
from time import sleep
import scipy.signal as sig
import ConfigParser
import fnmatch
from collections import OrderedDict


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


class Parameters:
    def __init__(self, initfile):
        self.cfg = ConfigParser.ConfigParser(
            {'aux_channel_source': '/virgoData/ffl/rds.ffl'},
            allow_no_value=True)
        self.cfg.read(initfile)

        self.aux_source = self.cfg.get('GENERAL', 'aux_channel_source')
        self.n_groups = self.cfg.getint('GENERAL', 'n_groups')
        self.all_aux = self.cfg.getboolean('GENERAL', 'all_aux')
        self.aux = self.cfg.get('GENERAL', 'aux_channels').split(
            '\n') if self.cfg.get('GENERAL', 'aux_channels') else None
        self.excluded = self.cfg.get('GENERAL', 'exclude').split(
            '\n') if self.cfg.get('GENERAL', 'exclude') else None
        self.nav = self.cfg.getint('GENERAL', 'averages')
        self.group_dict = OrderedDict()
        self.aux_dict = {}
        for group_n in xrange(self.n_groups):
            sctn = 'GROUP' + str(group_n + 1)
            self.group_dict[sctn] = {'channel': self.cfg.get(sctn, 'channel'),
                                     'units': self.cfg.get(sctn, 'units'),
                                     'bands': self.cfg.get(sctn, 'bands'),
                                     'aux': self.cfg.get(sctn,
                                                         'aux_channels').split(
                                         '\n') if self.cfg.get(sctn,
                                                               'aux_channels') else None,
                                     'excl': self.cfg.get(sctn,
                                                          'exclusions').split(
                                         '\n') if self.cfg.get(sctn,
                                                               'exclusions') else None
                                     }
        if self.cfg.has_section('RESULTS'):
            self.res_param = {}
            for option in self.cfg.options('RESULTS'):

                self.res_param[option] = self.cfg.get('RESULTS', option)
                # Try to convert them to integers
                # todo: is there a better way to check if they can be converted?
                try:
                    self.res_param[option] = int(self.res_param[option])
                except ValueError:
                    pass

    def save_extended_config(self, **kwargs):
        # todo: other interesting parameters, such as averages, npoints,
        if not self.cfg.has_section('RESULTS'):
            self.cfg.add_section('RESULTS')
        for key, arg in kwargs.iteritems():
            self.cfg.set('RESULTS', key, arg)
        with open(kwargs['hdir'] + 'config.ini', mode='w') as f:
            self.cfg.write(f)

    def extract_aux_channels(self, gpsb):
        if self.all_aux:
            chlist = get_channel_list(gpsb, self.aux_source)
            for aux in chlist:
                self.aux_dict[aux] = set(self.group_dict.keys())
        for excl_name in self.excluded or []:
            try:
                del (self.aux_dict[excl_name])
            except KeyError:
                for aux in self.aux_dict.iterkeys():
                    if (fnmatch.fnmatch(aux, 'V1:' + excl_name + '*') or
                       fnmatch.fnmatch(aux, excl_name + '*')):
                        self.excluded.append(aux)
        for aux in self.aux or []:
            self.aux_dict[aux] = set(self.group_dict.keys())
        for group_name, dic in self.group_dict.iteritems():
            for excl_name in dic['excl'] or []:
                try:
                    self.aux_dict[excl_name].remove(group_name)
                except KeyError:
                    for aux in self.aux_dict.iterkeys():
                        if (fnmatch.fnmatch(aux, 'V1:' + excl_name + '*') or
                           fnmatch.fnmatch(aux, excl_name + '*')):
                            try:
                                self.aux_dict[aux].remove(group_name)
                            except KeyError:
                                pass
            for aux in dic['aux'] or []:
                try:
                    self.aux_dict[aux].add(group_name)
                except KeyError:
                    self.aux_dict[aux] = {group_name}


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


def extractbands(g_dict):
    g_dict['band_list'] = []
    b = g_dict['bands'].split(':')
    for i in range(len(b) - 1):
        g_dict['band_list'].append(b[i] + '_' + b[i + 1] + 'Hz')


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
        b = sig.firwin(n + 1, 1. / q, window='hamming')
        y = sig.lfilter(b, 1., x, axis=axis)
    else:
        (b, a) = sig.cheby1(n, 0.05, 0.8 / q)

        y = sig.lfilter(b, a, x, axis=axis)
    return y.swapaxes(0, axis)[::q].swapaxes(0, axis)


def decimator_wrapper(downfreq, ch):
    if ch.fsample < downfreq:
        print "Channels must have a larger or equal sampling frequency than " \
              "the desired downsampled freq"
        raise ValueError
    else:
        if ch.fsample % downfreq != 0:
            print "Sampling frequency must be equal or an integer " \
                  "multiple of the downsampled frequency"
            raise ValueError
        else:
            if ch.fsample > downfreq:
                decimfactor = ch.fsample / downfreq
                # print "Decimation factor {0}".format(decimfactor)
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
