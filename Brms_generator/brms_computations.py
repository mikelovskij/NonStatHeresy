import numpy as np
import virgotools as vrg
import scipy.signal as sig
from functions import decimator_wrapper
__author__ = 'Michele Valentini snao20@hotmail.com'


# Class to store and process the parameters of a channel
class ChannelParameters:
    def __init__(self, channel, parameters):
        self.channel = channel
        self.ds_freq = 2 * parameters.max_freq
        self.n_points = int(parameters.n_points)
        self.f_res = float(self.ds_freq) / self.n_points
        self.overlap = parameters.overlap
        # Number of non-overlapped points in each FFT window
        self.n_over_points = int(np.rint(self.n_points * (1 - self.overlap)))
        self.times = np.empty(0, dtype='float64')
        self.steps = parameters.step_dict
        self.bands = parameters.channels_bands[channel]
        self.steps = parameters.step_dict
        self.steps["BrmsComputation"]['bands'] = self.bands
        self.steps["BrmsComputation"]['f_res'] = self.f_res

    def compute_n_fft(self, segment_points):
        segment_points = int(segment_points)
        # since the inputs are ints, the division also truncates the result
        # as intended.
        return (segment_points - self.n_points) / self.n_over_points + 1

    def compute_time_vector(self, fft_per_segment, gpsb):
        # index where the last fft window begins
        last_index = (fft_per_segment - 1) * self.n_over_points
        # gps corresponding to this index
        last_gps = gpsb + float(last_index) / self.ds_freq
        # todo: should times be part of the channel class or of the progbrms ?
        # todo: use the time of the middle of the windows instead of the begin?
        time_discontinuity = True
        try:
            dt = self.times[-1] - self.times[-2]
            if gpsb <= (self.times[-1] + dt):
                time_discontinuity = False
        except IndexError:
            time_discontinuity = True

        self.times = np.concatenate((self.times, np.linspace(gpsb, last_gps,
                                                             fft_per_segment,
                                                             endpoint=False)))
        return time_discontinuity


def process_channel(ch_p, src, segments):
    step_instances = []
    for step in ch_p.steps.itervalues():
        step_instances.append(step['class'](step))
    brms_buffer = []
    for (j, (gpsb, gpse)) in enumerate(segments):
        print 'segment {} of {}'.format(j, len(segments))
        with vrg.getChannel(src, ch_p.channel, gpsb, gpse - gpsb) as r_data:
            # Decimate the channel data
            ch_data = decimator_wrapper(ch_p.ds_freq, r_data)
        fft_per_segment = ch_p.compute_n_fft(len(ch_data))
        # Estimate the time vector of the brms samples for this seg
        t_disc = ch_p.compute_time_vector(fft_per_segment, gpsb)
        # Reset the buffers if there is a time discontinuity between
        # this segment and the previous one.
        if t_disc:
            for step in step_instances:
                step.reset()
        # make a function with the spectra as output?
        # il ciclo su fft_per segment lo terrei cmq qui. o no?
        for i in xrange(fft_per_segment):
            freqs, s = sig.welch(ch_data[i * ch_p.n_over_points:
                                         i * ch_p.n_over_points +
                                         ch_p.n_points],
                                 fs=ch_p.ds_freq,
                                 window='hanning',
                                 nperseg=ch_p.n_points)
            pipe = (freqs, s)
            for step in step_instances:

                pipe = step(pipe)
                # todo:use an hdf5 buffer between segments?probably unnecessary
            brms_buffer.append(pipe[1])
            real_bands = pipe[0]
    # transpose the buffer in order for it to be used
    brms_buffer = np.array(brms_buffer).T.tolist()
    return brms_buffer, real_bands


# factor a number
def factors(n):
    return np.sort(reduce(list.__add__,
                ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))[1:-1]

