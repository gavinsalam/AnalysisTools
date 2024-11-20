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

    parser.add_argument("--only", type=str, help='Plot only histograms who name matches this regexp')

    parser.add_argument("--yrange", metavar="min,max", type=str, default = "", 
                        help="y range for all plots")
    parser.add_argument("--xrange", metavar="min,max", type=str, default = "", 
                        help="x range for all plots")
    parser.add_argument("--xlabel", type=str, help="x-axis label")
    parser.add_argument("--ylabel", type=str, help="y-axis label")
    parser.add_argument("--pdfname", "-o","--out","-out", metavar='pdfname', type=str, default = "", 
                        help='Name of the output PDF file (if not specified, it is derived from first filename)')

    args = parser.parse_args()


    filenames = args.filenames
    hfiles = []
    for filename in filenames:
        parts = [f.strip() for f in filename.split(' + ')]
        hfile = h.HFile(parts[0])
        if len(parts) > 1:
            for part in parts[1:]:
                hfile += HFile(part)
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
    
    # make plot
    with PdfPages(pdffile) as pdf: 
        for ih, histogram in enumerate(hfiles[0].histograms):
            if args.only and not re.search(args.only, histogram.name): continue
            print("Plotting histogram", histogram.name)
            fig,ax = plt.subplots()
            # get minor ticks to show up with automatic spacing
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            # set the y range if requested
            if args.xrange: ax.set_xlim([float(y) for y in args.xrange.split(',')])
            if args.yrange: ax.set_ylim([float(y) for y in args.yrange.split(',')])
            histogram.set_axes_data(ax)
            if args.xlabel: ax.set_xlabel(args.xlabel)
            if args.ylabel: ax.set_ylabel(args.xlabel)

            if args.logy and re.search(args.logy, histogram.name): ax.set_yscale('log')

            for ihfiles in range(nfiles):
                hh = hfiles[ihfiles].by_name(histogram.name)
                # print(hh)
                # hh = hfiles[ihfiles].histograms[ih]
                # if hh.name != histogram.name:
                #     raise ValueError(f"histogram {histogram.name} not found in file {filenames[ihfiles]}")
                if (args.norm): 
                    hh.plot_to_axes(ax, **styles[ihfiles%nstyles], norm=hfiles[0].histograms[ih]) 
                else:
                    hh.plot_to_axes(ax, **styles[ihfiles%nstyles]) 

            ax.legend(loc='best')
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