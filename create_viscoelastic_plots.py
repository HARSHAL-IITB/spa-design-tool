#! /usr/bin/env python
# The MIT License (MIT)
#
# Copyright (c) 2015, EPFL Reconfigurable Robotics Laboratory,
#                     Philip Moseley, philip.moseley@gmail.com
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

from matplotlib import pyplot,axes
import calc_viscoelastic_parameters as CVP
import model_utility as MU
import numpy as np
import argparse, numpy

#--------------------------------------------------------------------------------
# Main.
#--------------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description="Create a comparison plot from viscoelastic fits.",
            epilog="Example: ./create_viscoelastic_plots.py --ylim 0.8 1 plot.eps fit1.mat fit2.mat")
    parser.add_argument("img",help="Output image to create (eps or png).")
    parser.add_argument("datafile",help="Dataset to fit, consisting of (time, shear mod, bulk mod) columns.")
    parser.add_argument("matfiles",nargs="+",help="List of .mat files to plot.")
    parser.add_argument("--xlim",nargs=2,type=float,default=[0.0,60.0],help="Min,Max values on xscale (0-60 if undefined).")
    # parser.add_argument("--ylim",nargs=2,type=float,help="Min,Max values on yscale (AUTO if undefined).")
    parser.add_argument("--paper",action='store_true',help="Create plots designed for the paper, rather than for general use.")
    args = parser.parse_args()

    # Mathematical rendering can be enabled by enclosing the symbols in $'s
    pyplot.rc('savefig',dpi=300)
    pyplot.rc('font',size=(14 if args.paper else 8))
    pyplot.rc('mathtext',default='regular') # Don't use italics for mathmode.

    # Read in the given datasets.
    print '--------------------------------------------------------------------'
    print ' Importing dataset...'
    data = np.loadtxt(args.datafile,comments='#',delimiter=',',dtype=np.longdouble)
    print '  Imported',data.shape[0],'datapoints.'

    # Read and plot each matfile.
    ax = []
    fig = pyplot.figure()
    ax.append(fig.add_subplot(111))
    ax.append(ax[0].twinx())
    for fname in args.matfiles:
        # Read file.
        mat = MU.read_viscomatfile(fname)
        popt = list()
        for i in range(len(mat['tau'])):
            popt.append(mat['G0']*mat['g'][i])  # Re-dimensionalize.
            popt.append(mat['tau'][i])
        # Plot.
        CVP.create_data_plots(ax,data,[mat['G0'],0.0]+popt,args.xlim,args.paper)
    if not args.paper: ax.legend(loc='best',frameon=False,framealpha=0)
    pyplot.tight_layout()
    pyplot.savefig(args.img)
    pyplot.close()
