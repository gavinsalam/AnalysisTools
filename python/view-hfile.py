#!/usr/bin/env python3
"""Matplotlib template generated automatically with 
/Users/gsalam/scripts/mptemplate.py hfile-viewer.py
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from   matplotlib.backends.backend_pdf import PdfPages
from   matplotlib.ticker import ScalarFormatter
from   matplotlib.ticker import FuncFormatter
from   matplotlib import cycler
from   matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import numpy as np
import copy
import hfile as h
import re
import sys
import argparse
from shlex import quote

styles = [
    {'color':'#f00000'},
    {'color':'#0000e0'},
    {'color':'#00c000'},
    {'color':'#404040'},
    {'color':'#e0a040'},
    {'color':'#00a0a0'},
    ]
nstyles = len(styles)
dashstyles = ['solid', 'dashed', 'dashdot', 'dotted']
ndashstyles = len(dashstyles)

colors = cycler('color', [style['color'] for style in styles])
# see options for things to set with 
#     python -c "import matplotlib as mpl; print(mpl.rcParams)"
plt.rc('axes', prop_cycle=colors)
mpl.rcParams.update({"axes.grid" : True, "grid.linestyle": ":"})
plt.rc('figure', figsize=(5,3))

# turn off divison by zero warnings
np.seterr(divide='ignore', invalid='ignore')


def main():
    parser = argparse.ArgumentParser(description='Plot histograms from hfile(s)')
    # add argument that gives any number of filenames to plot
    parser.add_argument('filenames', metavar='filename', type=str, nargs='+',
                        help='hfile(s) to plot')
    parser.add_argument("--logy", metavar="regexp", type=str, default = "", 
                        help="use log scale for anything that matches the regexp")
    parser.add_argument("--norm", help="plot histograms normalised to the first file", action="store_true")
    parser.add_argument("--value-and-ratio","-r", help="plot histograms and a version normalised to the first file", action="store_true")

    parser.add_argument("--only", type=str, help='Plot only histograms whose name matches this regexp')
    parser.add_argument("--exclude", "-v", type=str, help='Do not plot histograms whose name matches this regexp')

    parser.add_argument("--yrange", metavar="min,max", type=str, default = "", 
                        help="y range for all plots")
    parser.add_argument("--ratio-range", metavar="min,max", type=str, default = "", 
                        help="y range in ratios for all plots")
    parser.add_argument("--xrange", metavar="min,max", type=str, default = "", 
                        help="x range for all plots")
    parser.add_argument("--xlabel", type=str, help="x-axis label")
    parser.add_argument("--ylabel", type=str, help="y-axis label")
    parser.add_argument("--pdfname", "-o","--out","-out", metavar='pdfname', type=str, default = "", 
                        help='Name of the output PDF file (if not specified, it is derived from first filename)')
    parser.add_argument("--steps", default=False, action="store_true", help="Plot histograms with steps instead of lines")
    parser.add_argument("--merge", default=None, type=str, help="regex of string to remove, e.g. '(good|bad)', "
                                                                "histograms with identical string post-removal are on the same plot")

    args = parser.parse_args()


    filenames = args.filenames
    hfiles = []
    for filename in filenames:
        parts = [f.strip() for f in filename.split(' + ')]
        hfile = h.HFile(parts[0])
        if len(parts) > 1:
            for part in parts[1:]:
                hfile += h.HFile(part)
        hfiles.append(hfile)
        print (f"{hfile.filename}: {hfile.histograms[0].value_array()[0]}")
    nfiles = len(filenames)

    #hfiles = [HFile(filename) for filename in filenames]

    print (f"===================  header from first file {hfiles[0].filename} ===================")
    print(hfiles[0].header)
    #for hfiles
    for ihfiles in range(len(hfiles)):
        print (f"===================  warnings file {hfiles[ihfiles].filename} ===================")
        print(hfiles[ihfiles].warnings)

    # get name of pdf
    if not args.pdfname: 
        pdffile = hfiles[0].filename + '.pdf'
    else:
        pdffile = args.pdfname 
    shfile = pdffile + ".sh"
    cmdline = " ".join([quote(arg) for arg in sys.argv])
    with open(shfile,'w') as sh: print(cmdline,file=sh)

    # handle the possibility of merging histograms onto a single page
    if not args.merge: histogram_lists = [[{'name': h.name,'common': h.name,'mergelabel': ''}] for h in hfiles[0].histograms]
    else:
        # if we are merging histograms onto a single page, we need to create
        # a new list of histograms with the merged names
        histogram_lists = []
        merged_idx = {}
        for ih, histogram in enumerate(hfiles[0].histograms):
            merge_common = re.sub(args.merge, "", histogram.name)
            merge_name = re.search(args.merge, histogram.name).group(0)
            hist_dict = {'name':histogram.name,'common':merge_common,'mergelabel': merge_name + ", "}
            if merge_common not in merged_idx:
                merged_idx[merge_common] = len(histogram_lists)
                histogram_lists.append([hist_dict])
            else:
                histogram_lists[merged_idx[merge_common]].append(hist_dict)

    # make plot
    with PdfPages(pdffile) as pdf: 
        norm = hfiles[0].by_name(histogram_lists[0][0]['name']) if args.norm or args.value_and_ratio else None

        for ih, histogram_list in enumerate(histogram_lists):
            h0name = histogram_list[0]['name']
            if args.only    and not re.search(args.only   , h0name): continue
            if args.exclude and     re.search(args.exclude, h0name): continue
            print (histogram_list)
            print("Plotting histograms", [h['name'] for h in histogram_list])
            if args.value_and_ratio:
                fig,(axh,ax) = plt.subplots(nrows=2, sharex = True, figsize = (5,6))
                # reduce space between subplots
                fig.subplots_adjust(hspace=0.05)
                if args.ratio_range: ax.set_ylim([float(y) for y in args.ratio_range.split(',')])
            else:
                fig,ax = plt.subplots()
                axh = ax # this alias helps avoid some ifs below

            # get minor ticks to show up with automatic spacing
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            # set the y range if requested
            if args.xrange: ax.set_xlim([float(y) for y in args.xrange.split(',')])
            if args.yrange: axh.set_ylim([float(y) for y in args.yrange.split(',')])
            hfiles[ihfiles].by_name(h0name).set_axes_data(ax)
            ax.set_title(histogram_list[0]['common'])
            ax.set_xlabel(re.sub(r'.*?:', '', histogram_list[0]['common']))
            if args.xlabel: ax.set_xlabel(args.xlabel)
            if args.ylabel: ax.set_ylabel(args.ylabel)

            if args.logy and re.search(args.logy, histogram.name): axh.set_yscale('log')

            
            for ihfiles in range(nfiles):
              for ih_entry, h_entry in enumerate(histogram_list):
                hh = hfiles[ihfiles].by_name(h_entry['name'])
                dashstyle=dashstyles[ih_entry%ndashstyles]
                hh.plot_args['label'] = h_entry['mergelabel'] + hh.plot_args['label']
                if args.norm: 
                    hh.plot_to_axes(ax, **styles[ihfiles%nstyles], ls=dashstyle, norm=norm) 
                elif args.value_and_ratio:
                    if ihfiles == 0 and ih_entry == 0: 
                        axh.set_title(ax.get_title())
                        axh.set_ylabel(ax.get_ylabel()+"")
                        ax.set_ylabel("ratio to first")
                        ax.set_title("")
                    hh.plot_to_axes(axh, **styles[ihfiles%nstyles], ls=dashstyle, steps=args.steps) 
                    hh.plot_to_axes(ax , **styles[ihfiles%nstyles], ls=dashstyle, norm=norm, steps=args.steps) 
                else:
                    hh.plot_to_axes(ax, **styles[ihfiles%nstyles], ls=dashstyle, steps=args.steps) 

            axh.legend(loc='best')
            pdf.savefig(fig,bbox_inches='tight')
            plt.close()
            #main(pdf)

    print("Wrote", pdffile)

    # now open the file in the default viewer
    # on mac use "open", on linux use "xdg-open"
    import subprocess
    # figure out which system we're on
    import platform
    system = platform.system()
    if system == 'Darwin':
        subprocess.run(["open",pdffile])
    else:
        subprocess.run(["xdg-open",pdffile])


if __name__ == "__main__": 
    main()