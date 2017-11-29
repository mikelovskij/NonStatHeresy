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
        self.aux_dict = savingdict['aux_dict']
        self.group_dict = savingdict['group_dict']


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
        copy2(os.path.dirname(os.path.abspath(__file__)) + '/' + sheet, hdir)
    copy2(os.path.dirname(os.path.abspath(__file__)) + '/' + script, hdir)

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
    for aux_name, aux_groups in data.aux_dict.iteritems():
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
                   for group in data.group_dict.keys())
    page.div(ol.button(data.group_dict.keys(), class_='tablinks',
                       onclick=onclick_gen), class_='tab')

    # build each group subpage
    for group, g_dict in data.group_dict.iteritems():
        # Build the highest ccf and coherence for this group table
        ccftab = {group: np.zeros(ntop)}
        cohtab = {group: np.zeros(ntop)}
        ccftab_names = {group: np.zeros(ntop, dtype='string')}
        cohtab_names = {group: np.zeros(ntop, dtype='string')}
        for aux_n, aux_name, aux_groups in enumerate(data.aux_dict.iteritems()):
            if group in aux_groups:
                for band in g_dict['band_list']:
                    if aux_n == 0:
                        ccftab[band] = np.zeros(ntop)
                        cohtab[band] = np.zeros(ntop)
                        ccftab_names[band] = np.zeros(ntop, dtype='string')
                        cohtab_names[band] = np.zeros(ntop, dtype='string')
                    ccf = data.ccfs[aux_name][group + '_' + band]
                    coh = data.mean_cohs[aux_name][group + '_' + band]
                    # todo: build single band tables too
                    if abs(ccf) > min(np.abs(ccftab[band])):
                        ccftab[band] = np.concatenate((ccftab[band], [ccf]))
                        ccftab_names[band] = np.concatenate((ccftab_names[band],
                                                             [aux_name]))
                        best_indexes = np.abs(ccftab[band]).argsort()[1:][::-1]
                        ccftab[band] = ccftab[band][best_indexes]
                        ccftab_names[band] = ccftab_names[band][best_indexes]
                    if coh > min(cohtab[band]):
                        cohtab[band] = np.concatenate((cohtab[band], [coh]))
                        cohtab_names[band] = np.concatenate((cohtab_names[band],
                                                             [aux_name]))
                        best_indexes = cohtab[band].argsort()[1:][::-1]
                        cohtab[band] = cohtab[band][best_indexes]
                        cohtab_names[band] = cohtab_names[band][best_indexes]

                    # Build the full group best tabs
                    if abs(ccf) > min(np.abs(ccftab[group])):
                        ccftab[group] = np.concatenate((ccftab[group], [ccf]))
                        ccftab_names[group] = np.concatenate((ccftab_names[group],
                                                             [aux_name +
                                                             ' with ' + band]))
                        best_indexes = np.abs(ccftab[group]).argsort()[1:][::-1]
                        ccftab[group] = ccftab[group][best_indexes]
                        ccftab_names[group] = ccftab_names[group][best_indexes]
                    if coh > min(cohtab[group]):
                        cohtab[group] = np.concatenate((cohtab[group], [coh]))
                        cohtab_names[group] = np.concatenate((cohtab_names[group],
                                                             [aux_name +
                                                             ' with ' + band]))
                        best_indexes = cohtab[group].argsort()[1:][::-1]
                        cohtab[group] = cohtab[group][best_indexes]
                        cohtab_names[group] = cohtab_names[group][best_indexes]
        tab = [[" CCF ", "Coherence"]]
        for i in xrange(ntop):
            row = [
                "{0}<br>CCFs = {1}".format(ccftab_names[group][i],
                                           ccftab[group][i]),
                "{0}<br>Mean Coher. = {1}".format(cohtab_names[group][i],
                                                  cohtab[group][i])]
            tab.append(row)
        tab_str = {group: tb.tabmaker(tab, True, False)}

        for band in g_dict['band_list']:
            tab = [[" CCF ", "Coherence"]]
            for i in xrange(ntop):
                row = [
                    "{0}<br>CCFs = {1}".format(ccftab_names[band][i],
                                               ccftab[band][i]),
                    "{0}<br>Mean Coher. = {1}".format(cohtab_names[band][i],
                                                      cohtab[band][i])]
                tab.append(row)
            tab_str[band] = tb.tabmaker(tab, True, False)

        # build the rest of the frame
        frame = ol.h1("NonStatMoni BRMS for {} GPS {:d} - {:d}".format(
                      g_dict['channel'], gpsb, gpse))
        frame += ol.h2("Normalized BRMS time series")
        # todo: normalized by what?
        img = ol.img(class_="myImg", src=pdir + group + "_time.png",
                     alt=group + "Time plot", width="400")
        frame += ol.div(img, style="float:left")
        onclick_gen = ("openBand(event, '{}')".format(band)
                       for band in g_dict['band_list'])

        frame += ol.div(ol.button("Group", class_='v_tablinks', id="defaultOpen",
                                  onclick="openBand(event, group)") +
                        ol.button(g_dict['band_list'], class_='v_tablinks',
                                  onclick=onclick_gen),
                        class_='vertical_tab')
        frame += ol.div(tab_str, class_="v_tabcontent",
                        id=([group]+g_dict['band_list']))
        #todo: inserire lle nuove tabelle
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
        for aux_name, aux_groups in data.aux_dict.iteritems():
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
