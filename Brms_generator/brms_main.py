import h5py
from config_manager import Parameters
from segment_manipulation import check_state_vec, segment_files_reader,\
    segment_splitter
from brms_computations import ChannelParameters, process_channel
from argparse import ArgumentParser
import os


# todo: add an argument parser
def main(config_file, savedir, segment_params=None, segmentfiles=None,
         max_segment_length=1200):
    par = Parameters(config_file)
    par.merge_bands()
    # Load the segments to use from segment files or
    # generate them from a channel and a threshold
    if segment_params:
        gpsb = segment_params[0]
        gpse = segment_params[1]
        if len(segment_params) > 2:
            state_channel = segment_params[2]
            threshold = segment_params[3]
            try:
                condition = segment_params[4]
            except IndexError:
                condition = 'greater_equal'
            segments = check_state_vec(state_channel, gpsb, gpse, threshold,
                                       condition)
        else:
            if len(segment_params) == 2:
                segments = [(gpsb, gpse)]
            else:
                # todo: Explain better.
                raise ValueError('You should provide two or more segment '
                                 'parameters')
    else:
        if segmentfiles:
            segments = segment_files_reader(segmentfiles[0], segmentfiles[1:])
            gpsb = segments[0][0]
            gpse = segments[-1][-1]
        else:
            raise ValueError('Need to Provide segment parameters'
                             ' or segment files')
    segments = segments.tolist()
    segment_splitter(segments, max_segment_length)
    print len(segments)
    f_res = float(2 * par.max_freq) / par.n_points
    res_dir = savedir + "/ {:d} - {:d} - resolution {:.3f} - npoints {}, overlap {:.2f}".format(gpsb, gpse, f_res, par.overlap, par.n_points)
    try:
        print "Creating directory {}".format(res_dir)
        os.mkdir(res_dir)
    except OSError as er:
        print er
    fname = res_dir + '/brms.hdf5'
    for channel, bands in par.channels_bands.iteritems():
        ch_p = ChannelParameters(channel, par)
        results = process_channel(ch_p, par.data_source, segments)

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        # +++++++++++++++++++++++ Save the results ++++++++++++++++++++++++++ #
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        with h5py.File(fname, 'a') as f:
            f.attrs.create('gps', [gpsb, gpse])
            f.attrs.create('fsample', data=float(ch_p.ds_freq) / (ch_p.n_points * (1 - ch_p.overlap)))
            f.attrs.create('segments', data=segments)
            try:
                g = f.create_group(channel)
            except ValueError:  # delete channel data if already exists and create anew
                del f[channel]
                g = f.create_group(channel)
            g.create_dataset('times', data=ch_p.times)
            for band, j in zip(results[1], xrange(len(results[1]))):
                g.create_dataset(bands[j], data=results[0][j])

# todo: write the resulting bands in the cfg?

if __name__ == "__main__":
    # read configuration options and process
    parser = ArgumentParser()
    parser.add_argument("-i", "--init", dest="initialization",
                        default='NonStatMoni.ini',
                        help="set initialization FILE", metavar="cfgFILE")
    parser.add_argument("-d", "--dir", dest="hdir",
                        help="output directory", metavar="OutDir")
    parser.add_argument("-b", "--gpsb", dest='gpsb',
                        help="start gps (used if not using segment files)",
                        metavar="gpsb", type=int)
    parser.add_argument("-e", "--gpse", dest='gpse',
                        help="end gps (used if not using segment files)",
                        metavar="gpse", type=int)
    parser.add_argument("-sc", "--statechannel", dest="state_channel",
                        help="state channel that can be used to generate time "
                             "segments that will be analyzed",
                        metavar="State_Channel")
    parser.add_argument("-st", "--threshold", dest="state_threshold",
                        help="Threshold that should be applied to the state"
                             " channel in order to generate segments",
                        metavar="State_Threshold", type=float)
    parser.add_argument("-scn", "--statecondition", dest="state_condition",
                        help="How the state channel should be compared with "
                             "the threshold i.e. equal or lesser_equal",
                        metavar="State_Condition", default="greater_equal")
    parser.add_argument("-ss", "--sciencesegments",
                        dest="science_segments_file",
                        help="Json formatted science segments file."
                             "If provided, gpsb, gpse, statechannel threshold"
                             " and state condition arguments will be ignored",
                        metavar="science_segment_file")
    parser.add_argument("-dqfl", "--dqflags", dest="dq_flag_files", nargs='*',
                        help="Json formatted data quality flag files. "
                             "Periods where any of these flags are active "
                             "will be excluded from analysis",
                        metavar="DQ-Flag_files")
    parser.add_argument("-pp", "--postprocess", action="store_true",
                        dest="post_process",
                        help="Launch the nonstatmoni post processing of the "
                             "brmss after the computation")
    args = parser.parse_args()
    try:
        print "Creating directory {}".format(args.hdir)
        os.mkdir(args.hdir)
    except OSError as err:
        print err
    if args.science_segments_file:
        print "Loading segment files, gpsb, gpse and state channel" \
              " arguments will be ignored"
        if args.dq_flag_files:
            main(args.initialization, args.hdir,
                 segmentfiles=([args.science_segments_file] +
                               args.dq_flag_files))
        else:
            main(args.initialization, args.hdir,
                 segmentfiles=[args.science_segments_file])
    else:
        if args.gpsb:
            if args.state_channel:
                main(args.initialization, args.hdir,
                     segment_params=[args.gpsb,
                                     args.gpse,
                                     args.state_channel,
                                     args.state_threshold,
                                     args.state_condition])
            else:
                main(args.initialization, args.hdir,
                     segment_params=[args.gpsb,
                                     args.gpse])

    if args.post_process:
        # todo: lancio nonstatmoni_main, ma prima devo creare un nonstatmoni main
        pass


    # todo: write the used arguments in the cfg file?
