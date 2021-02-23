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

import argparse, sys, os, shutil, time
import os.path as path
import numpy as np
import utility
from matplotlib import pyplot,axes


#--------------------------------------------------------------------------------
# Create the plots from extracted experimental data files.
#   axis      = axis handle of figure.
#   fname     = file to plot.
#   test      = type of test to analyse.
#   line      = choose a linetype for all lines from this datafile.
#   labelpref = prefix for labels.
#--------------------------------------------------------------------------------
def plot_exp_data(axis, fname, test, line='', labelpref=''):
    # Load the file.
    try:
        if   fname.endswith('.csv'): alldata = np.loadtxt(fname,delimiter=',')
        elif fname.endswith('.txt'): alldata = np.loadtxt(fname)
    except:
        utility.print_error('Problem reading experimental data file: '+fname,True)

    cm = pyplot.get_cmap('gist_rainbow')
    if test=='linU':
        axis.plot([2.0,2.0], [0.0, 50.0], '-', color='0.35',linewidth=3)
        if labelpref!='': labelpref = labelpref+' - '
        npress = alldata.shape[1]/2
        colors = [cm(float(i)/npress) for i in range(npress)]
        for i in range(npress):
            # Note that we skip the first row, since it contains the pressure data.
            # label = labelpref + str(alldata[0,2*i]) + r' ${}^N\!/{}_{mm^2}$'
            label = labelpref + str(int(1000*alldata[0,2*i])) + ' kPa'
            axis.plot(alldata[1:,2*i], alldata[1:,2*i+1], line, color=colors[i], label=label)
    elif test=='linF' or test=='bndF':
        if labelpref!='': labelpref = labelpref+' - '
        pressures = np.unique(alldata[:,2])
        npress = len(pressures)-1       # Don't seperately plot the 0 pressure.
        colors = [cm(float(i)/npress) for i in range(npress)]
        T = []
        F = []
        P = alldata[0,2]
        for i in range(alldata.shape[0]):
            if alldata[i,2]==P or alldata[i,2]==0.0:
                P = max(alldata[i,2],P)
                T.append(alldata[i,0])
                F.append(alldata[i,1])
            else:
                label = labelpref + str(P) + r' ${}^N\!/{}_{mm^2}$'
                cidx = np.where(pressures==P)[0]-1
                axis.plot(T, F, line, color=colors[cidx], label=label)
                T = [alldata[i,0]]
                F = [alldata[i,1]]
                P = alldata[i,2]
        label = labelpref + str(P) + r' ${}^N\!/{}_{mm^2}$'
        cidx = np.where(pressures==P)[0]-1
        axis.plot(T, F, line, color=colors[cidx], label=label)
    elif test=='bndU':
        axis.plot(alldata[:,0], alldata[:,1], line, label=labelpref)
    else:
        utility.print_error('Invalid test type in plot_exp_data.',True)



#--------------------------------------------------------------------------------
# Main.
#--------------------------------------------------------------------------------
if __name__ == "__main__":
    tinit = time.time()
    # Handle user input.
    parser = argparse.ArgumentParser(description='Create plots of experimental actuator results versus time. Image can be png or eps (or others).',
                                     epilog='Example: create_experimental_plots.py --sort --ylim 0 10 linU linear-U.png ../../data/*.csv')
    parser.add_argument('test',choices=['linF','linU','bndF','bndU','unknown'],help='Type of test to analyse.')
    parser.add_argument('img',help='Output image to create.')
    parser.add_argument('data',nargs='+',help='Data file(s) to plot (CSV or TXT files).')
    parser.add_argument('--sort',action='store_true',help='Numerically sort inputs from smallest to largest.')
    parser.add_argument('--xlim',nargs=2,type=float,metavar=('MIN','MAX'),help='Optional min and max for x-axis.')
    parser.add_argument('--ylim',nargs=2,type=float,metavar=('MIN','MAX'),help='Optional min and max for y-axis.')
    parser.add_argument('--title',help='Optional title for plot.')
    parser.add_argument("--paper",action='store_true',help="Create plots designed for the paper, rather than for general use.")
    args = parser.parse_args()
    pyplot.rc('mathtext',default='regular') # Don't use italics for mathmode.

    # Try to sort the data.
    expdata = args.data
    if args.sort:
        try:
            from natsort import natsort
            expdata = natsort(expdata)
        except:
            print 'WARNING: no natsort module found, sorting not available.'

    # Prepare the plots.
    fig,ax = pyplot.subplots()

    # TODO - remove longest common prefix from label strings.
    # Plot the EXP data.
    lines = ['-', '--', '-.', ':']
    if args.test=='linU':
        xlabel = 'Time (s)'
        ylabel = 'Displacement (mm)'
        title = 'Displacements vs Time'
    elif args.test=='linF':
        xlabel = 'Time (s)'
        ylabel = 'Blocked Force (N)'
        title = 'Blocked Forces vs Time'
    elif args.test=='bndU':
        lines = ['-']
        xlabel = r'Pressure $\left({}^N\!/{}_{mm^2}\right)$'
        ylabel = 'Bending Angle ($^\circ$)'
        title = 'Bending Angle vs Pressure'
    elif args.test=='bndF':
        xlabel = 'Time (s)'
        ylabel = 'Blocked Force (N)'
        title = 'Blocked Forces vs Time'

    for i in range(len(expdata)):
        print 'Reading data for:',expdata[i]
        lidx = i%len(lines)
        if not args.paper: labelpref = path.splitext(expdata[i])[0]
        else:              labelpref = ''
        plot_exp_data(ax, expdata[i], args.test, lines[lidx], labelpref)

    # Write the plot to disk.
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid()
    if args.paper:
        lgd = ax.legend(loc='best',frameon=False,framealpha=0)
    else:
        lgd = ax.legend(loc=2,bbox_to_anchor=(1,1),frameon=False,framealpha=0)
    for legobj in lgd.legendHandles: legobj.set_linewidth(2.0)
    if args.xlim: ax.set_xlim(args.xlim)
    if args.ylim: ax.set_ylim(args.ylim)
    if args.title:       ax.set_title(args.title)
    elif not args.paper: ax.set_title(title)
    print 'Writing plot to:',args.img
    pyplot.savefig(args.img,bbox_extra_artists=(lgd,), bbox_inches='tight')
    pyplot.close()
