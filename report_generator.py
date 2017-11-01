from collections import OrderedDict
import markup
from markup import oneliner as ol
import tables as tb
import numpy as np
from functions import Parameters
import cPickle
from argparse import ArgumentParser
from shutil import copy2
import os


class Data:
    def __init__(self, filename):
        with open(filename) as f:
            savingdict = cPickle.load(f)
        self.ccfs = savingdict['ccfs']
        self.mean_cohs = savingdict['mean_cohs']


def main(initfile, ntop):
    # read configuration options and process
    par = Parameters(initfile)
    gpsb = par.res_param['gpsb']
    gpse = par.res_param['gpse']
    pdir = par.res_param['pdir']
    hdir = par.res_param['hdir']
    data = Data(hdir + 'post_proc_data.dat')

    # relative paths of the stilesheets and javascripts for the results page
    stylesheets = ('general.css', 'table.css', 'modal.css')
    script = 'scripts.js'

    # copy them in the results folder (that should already exist)
    for sheet in stylesheets:
        copy2(os.path.dirname(__file__) + sheet, hdir)
    copy2(os.path.dirname(__file__) + script, hdir)

    # Create summary web page
    print "Generating Web Page..."
    page = markup.page()
    page.init(title="NonStatMoni", css=stylesheets,
              footer="(2017)" +
                     ol.a("Michele Valentini", href='mailto:snao20@hotmail.com'))

    # best coherences and ccfs summary table generation
    page.button("Auxillary channels summary table", class_="accordion")
    page.div(class_="panel")
    page.input(type="text", id="myInput", onkeyup="myFunction()",
               placeholder="Search for aux_channel names..")
    page.table(id="t01")
    for aux_name, aux_groups in par.aux_dict.iteritems():
        page.tr()
        page.td()
        page.h4(aux_name)
        page.td.close()
        # Sort the dictionary according to the abs of the values,
        #  from highest to lowest
        ccf_sorted = OrderedDict(sorted(data.ccfs[aux_name].iteritems(),
                                        key=lambda x: abs(x[1]), reverse=True))
        mean_coh_sorted = OrderedDict(sorted(data.mean_cohs[aux_name].iteritems(),
                                             key=lambda x: x[1], reverse=True))

        for chan_band, ccf in ccf_sorted.items()[0:4]:
            pagename = 'cohe_{}.html'.format(chan_band.split('_')[0])
            page.td(style="width:10%; font-size :11px;")

            page.add("<a target=_blank href={}>{}</a><br>(ccf={:3f})"
                     .format(pagename, chan_band, ccf))
            page.td.close()

        page.td()
        page.add('<a>      ********     </a>')
        page.td.close()

        for chan_band, coh in mean_coh_sorted.items()[0:4]:
            pagename = 'cohe_{}.html'.format(chan_band.split('_')[0])
            page.td(style="width:10%; font-size:11px;")
            page.add('<a target=_blank href={}>{}</a><br>(mncoh={:3f})'
                     .format(pagename, chan_band, coh))
            page.td.close()

        page.tr.close()
    page.table.close()
    page.div.close()

    # build the page menu
    onclick_gen = ("openGroup(event, '{}')".format(group)
                   for group in par.group_dict.keys())
    page.div(ol.button(par.group_dict.keys(), class_='tablinks',
                       onclick=onclick_gen), class_='tab')

    # build each group subpage
    for group, g_dict in par.group_dict.iteritems():
        # Build the highest ccf and coherence for this group table
        ccftab = np.zeros(ntop)
        cohtab = np.zeros(ntop)
        ccftab_names = np.zeros(ntop, dtype='string')
        cohtab_names = np.zeros(ntop, dtype='string')
        for aux_name, aux_groups in par.aux_dict.iteritems():
            if group in aux_groups:
                for band in g_dict['band_list']:
                    ccf = data.ccfs[aux_name][group + '_' + band]
                    if abs(ccf) > min(np.abs(ccftab)):
                        ccftab = np.concatenate((ccftab, [ccf]))
                        ccftab_names = np.concatenate((ccftab_names, [aux_name]))
                        best_indexes = np.abs(ccftab).argsort()[1:][::-1]
                        ccftab = ccftab[best_indexes]
                        ccftab_names = ccftab_names[best_indexes]
                    coh = data.mean_cohs[aux_name][group + '_' + band]
                    if coh > min(cohtab):
                        cohtab = np.concatenate((cohtab, [coh]))
                        cohtab_names = np.concatenate((cohtab_names, [aux_name]))
                        best_indexes = cohtab.argsort()[1:][::-1]
                        cohtab = cohtab[best_indexes]
                        cohtab_names = cohtab_names[best_indexes]
        tab = [[" CCF ", "Coherence"]]
        for i in xrange(ntop):
            row = [
                "{0}<br>Highest CCFs = {1}".format(ccftab_names[i], ccftab[i]),
                "{0}<br>Mean Coher. = {1}".format(cohtab_names[i], cohtab[i])]
            tab.append(row)
        tabstr = tb.tabmaker(tab, True, False)

        # build the rest of the frame
        frame = ol.h1("NonStatMoni BRMS for {} GPS {:d} - {:d}".format(
                      g_dict['channel'], gpsb, gpse))
        frame += ol.h2("Normalized BRMS time series")
        # todo: normalized by what?
        img = ol.img(class_="myImg", src=pdir + group + "_time.png",
                     alt=group + "Time plot", width="400")
        frame += ol.div(img, style="float:left")
        frame += ol.div(tabstr)
        frame += ol.h2("Spectrum of BRMS time series")
        frame += ol.img(class_="myImg", src=pdir + group + "_psd.png",
                        alt=group + "PSD plot", width="400")
        cohe_page_name = 'cohe_{}.html'.format(group)
        frame += ol.h2(ol.a("Coherences with slow channels", target='_blank',
                            href=cohe_page_name))
        page.div(frame, id=group, class_='tabcontent')

        # create coherence subpage
        page2 = markup.page()
        page2.init(title="NonStatMoni", css='../style/style.css')
        page2.h1("Coherences with slow channels")
        for aux_name, aux_groups in par.aux_dict.iteritems():
            if group in aux_groups:
                page2.img(src=(pdir + group + '_' + 'cohe_' +
                               aux_name.split(':')[1] + '.png'),
                          alt=(" No plots for" + aux_name))
        page2.savehtml(hdir + cohe_page_name)

    # create the modal for the plot images
    modal = ol.span("&times;", class_="close")
    modal += ol.img(class_="modal-content", id="img01")
    modal += ol.div('', id="caption")
    page.div(modal, id='myModal', class_="modal")

    page.br()
    page.h2("Contacts")
    page.scripts({script: 'javascript'})
    page.savehtml(hdir + 'index.html')


# Parse the system arguments and execute the main if it is called from bash
if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--init", dest="initialization",
                        default='NonStatMoni.ini',
                        help="set initialization FILE", metavar="cfgFILE")
    parser.add_argument("-n", "--ntop", dest="ntop",
                        default='5', type=int,
                        help="number of top coherences and correlation "
                             "cohefficients to be shown in the result table",
                        metavar="nTop")
    args = vars(parser.parse_args())
    main(args['initialization'], args['ntop'])
