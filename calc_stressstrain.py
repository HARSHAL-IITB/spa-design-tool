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

import argparse, time, csv, sys, inspect
import os.path as path
import numpy as np
import utility
from matplotlib import pyplot,axes,collections


#--------------------------------------------------------------------------------
# Remove data from first NCC cycles due to Mullin's Effect.
# Remove data collected during strain relaxation.
#   data = experimental datapoints.
#   NCC  = number of cycles to use.
#--------------------------------------------------------------------------------
def crop_data(data,NCC):
    if NCC < 0: return data,[]
    # Reduce noise by only using 1 point of out PPS.
    sign = 0
    PPS = 2
    for i in range(10):
        cycles = []
        samples = range(0,data.shape[0],PPS)
        strain_diff = np.sign(np.diff(data[samples,0]))
        for i,s in enumerate(strain_diff):
            if s!=sign:
                sign = s
                cycles.append(PPS*i)
        cycles.append(data.shape[0])
        NC = (len(cycles)-1.0)/2.0
        print 'Found cycles:             ',NC
        if NC>20 or NC!=int(NC):
            print '\tWARNING: Failed to find cycles with PPS =',PPS
            PPS = PPS+1
        else: break
    else:
        utility.print_error('Failed to find cycles.',False)
        print '\n'
        return data,[]

    print 'Number of cycles cropped: ',NC-NCC
    for i in range(len(cycles)-1,1,-2):
        if i==2.0*(NC-NCC):
            indices = range(0,cycles[i]+PPS)
            data = np.delete(data,indices,axis=0)
            break
        else:
            indices = range(cycles[i-1]-PPS,cycles[i]+PPS)
            data = np.delete(data,indices,axis=0)
    return data,cycles


#--------------------------------------------------------------------------------
# Calculate the strain and stress from force and area.
#   data = experimental datapoints.
#   A    = applied-force area.
#--------------------------------------------------------------------------------
def calc_strain_stress(data,A):
    print 'Cross-sectional area:     ',A
    min_idx = np.argmin(np.abs(data),axis=0)    # Index of minimum force value.
    z = data[min_idx[1],0]                      # Length at minimum force.
    print 'Minimum (abs val) force:  ',data[min_idx[1],1]
    print 'Length at minimum force:  ',z
    data[:,0] = (data[:,0] - z)/z       # Calculate (engineering) strain.
    #  data[:,0] = np.log(data[:,0]/z)  # Calculate (logarithmic) strain.
    data[:,1] = data[:,1]/A             # Calculate nominal stress.
    return data


#--------------------------------------------------------------------------------
# Plot the force or stress versus time.
#   ax     = axes from figure to plot on.
#   dataO  = untouched experimental datapoints.
#   data   = cropped experimental datapoints.
#   plabel = label for plot.
#--------------------------------------------------------------------------------
def plot_vs_time(ax,dataO,data,plabel):
    ax[0].plot(dataO[:,1],label=plabel)
    ax[0].set_title('Untouched Experimental Data')
    ax[0].set_xlabel('Sample #')
    ax[0].set_ylabel('Force (N)')
    ax[1].plot(data[:,1],label=plabel)
    ax[1].set_title('Cleaned Experimental Data')
    ax[1].set_xlabel('Sample #')
    ax[1].set_ylabel(r'Nominal Stress $\left({}^N\!/{}_{mm^2}\right)$')
    pyplot.tight_layout()


#--------------------------------------------------------------------------------
# Plot the force vs length and stress vs strain.
#   ax     = axes from figure to plot on.
#   cycles = bounds for each cycle, [] for no coloring.
#   dataO  = untouched experimental datapoints.
#   data   = cropped experimental datapoints.
#   plabel = label for plot.
#--------------------------------------------------------------------------------
def plot_vs_length(ax,cycles,dataO,data,plabel):
    if len(cycles)>0:
        idx = 1
        colors = iter(pyplot.cm.rainbow(np.linspace(0,1,(len(cycles)-1.0)/2.0,endpoint=True)))
        for i in range(0,len(cycles)-1,2):
            I = range(cycles[i],cycles[i+2])
            C = next(colors)
            W = 4.0 - 3.0*i/(len(cycles)-3.0)
            ax[0].plot(dataO[I,0],dataO[I,1],linewidth=W,color=C)
            idx = idx+1
    else:
        ax[0].plot(dataO[:,0],dataO[:,1])
    ax[0].set_title('Force-Length Curves')
    ax[0].set_xlabel('Length (mm)')
    ax[0].set_ylabel('Force (N)')
    ax[1].plot(data[:,0],data[:,1],label=plabel)
    ax[1].set_title('Stress-Strain Curves')
    ax[1].set_xlabel('Engineering Strain')
    ax[1].set_ylabel(r'Nominal Stress $\left({}^N\!/{}_{mm^2}\right)$')
    x0,x1 = ax[0].get_xlim()
    y0,y1 = ax[0].get_ylim()
    ax[0].set_aspect(abs(x1-x0)/abs(y1-y0))
    x0,x1 = ax[1].get_xlim()
    y0,y1 = ax[1].get_ylim()
    ax[1].set_aspect(abs(x1-x0)/abs(y1-y0))
    pyplot.tight_layout()


#--------------------------------------------------------------------------------
# Main.
#--------------------------------------------------------------------------------
if __name__ == "__main__":
    # Handle user input.
    # TODO - add option to set point-width to reduce bias.
    # TODO - write some of the screen output into the file header.
    parser = argparse.ArgumentParser(description="Calculate stress and strain given force and length for a"
            "given dataset. Length in testing direction is read from the datafile name if possible.",
            epilog="Example: calc_stressstrain.py --crop_cycles=2 --area=100 test1.dat test2.dat")
    parser.add_argument("-f","--cmbimg",metavar="FNAME",default="combined",help="Filename for combined plots (default \'combined\').")
    parser.add_argument("-t","--title",help="Title for combined plots.")
    parser.add_argument("-c","--crop_cycles",type=int,default=1,help="Number of cycles to use, considering Mullin's Effect (default 1). Disable with -1.")
    parser.add_argument("-a","--area",type=float,help="Area over which force was applied. Default is to read from filenames.")
    parser.add_argument("--format",choices=['png','eps'],help="Image format, default is png.")
    parser.add_argument("datafiles",nargs='+',help="Dataset(s) to process.")
    args = parser.parse_args()
    if args.format: fmt = args.format
    else:           fmt = 'png'
    pyplot.rc('mathtext',default='regular') # Don't use italics for mathmode.
    #  pyplot.rc('font',size=16)
    pyplot.rc('axes',grid=True)
    pyplot.rc('figure',dpi=300)
    pyplot.rc('savefig',dpi=300)

    # Remove any datafiles that have already been processed.
    datafiles = []
    for dfile in args.datafiles:
        if not path.splitext(dfile)[0].endswith('--CLEAN'): datafiles.append(dfile)

    # Try to sort the data.
    try:
        from natsort import natsort
        datafiles = natsort(datafiles)
    except:
        print 'WARNING: no natsort module found, sorting not available.'

    # Prepare the multiplots.
    if len(datafiles)>1:
        cm = pyplot.get_cmap('gist_rainbow')
        colors = [cm(float(i)/len(datafiles)) for i in range(len(datafiles))]
        fig1,ax1 = pyplot.subplots(3,1,figsize=(15,11))
        fig2,ax2 = pyplot.subplots(1,3,figsize=(15,6))
        if args.title:
            fig1.suptitle(args.title,fontweight="bold")
            fig2.suptitle(args.title,fontweight="bold")
        ax1[2].axis('off')
        ax2[2].axis('off')
        ax1[0].set_color_cycle(colors)
        ax1[1].set_color_cycle(colors)
        ax2[0].set_color_cycle(colors)
        ax2[1].set_color_cycle(colors)
    # Iterate over datafiles.
    for datafile in datafiles:
        # Read in the given dataset.
        print '--------------------------------------------------------------------'
        print ' Loading datafile',datafile
        print '  (assuming column0=time, column1=length, column2=force, additional columns ignored)'
        print '--------------------------------------------------------------------'
        try:
            dataO = np.loadtxt(datafile,usecols=(1,2),dtype=np.longdouble)
        except:
            utility.print_error('Failed to read file.',False)
            print '\n'
            continue
        print 'Number of points imported:',dataO.shape[0]
        fname = path.splitext(datafile)
        data = np.zeros(dataO.shape)
        data[:,0] = dataO[:,0]
        data[:,1] = dataO[:,1]

        # Calculate/determine the area over which the force is applied.
        if args.area: A = args.area
        else:
            try:
                size = fname[0].split('--')[2]
                dims = size.split('x')
                if len(dims)>2: raise Exception('Too many dimensions.')
                A = float(dims[0])*float(dims[1])
            except:
                utility.print_error('Failed to read dimensions from filename, use --area option.',True)

        # Remove datapoints with negative forces (compression).
        data = data[data[:,1]>0.1]
        # Clean and prepare the data.
        data,cycles = crop_data(data,args.crop_cycles)
        if not args.crop_cycles<0 and (len(cycles)==0 or len(data)==0):
            utility.print_error('Failed to extract data!',False)
            continue
        data = calc_strain_stress(data,A)

        # Plot the force or stress versus time.
        print 'Plotting data...'
        if len(datafiles)>1:
            pyplot.figure(fig1.number)
            plot_vs_time(ax1,dataO,data,fname[0])
        fig,ax = pyplot.subplots(2)
        plot_vs_time(ax,dataO,data,'none')
        pyplot.savefig(fname[0]+'--time.'+fmt)
        pyplot.close()

        # Sort the data.
        data = data[data[:,0].argsort()]

        # Plot the force vs length and stress vs strain.
        if len(datafiles)>1:
            pyplot.figure(fig2.number)
            plot_vs_length(ax2,[],dataO,data,fname[0])
        fig,ax = pyplot.subplots(1,2)
        plot_vs_length(ax,cycles,dataO,data,'none')
        pyplot.savefig(fname[0]+'--stress_strain.'+fmt)
        pyplot.close()

        # Write the new dataset to disk.
        outfile = fname[0]+'--CLEAN'+fname[1]
        head = 'Nominal stress (col0) and engineering strain (col1) data from '+datafile
        data = np.fliplr(data)              # Swap columns, now [stress,strain].
        np.savetxt(outfile,data,delimiter=',',header=head)
        print 'Number of points exported:',data.shape[0]
        print 'Data written to disk at:  ',outfile
        print
    if len(datafiles)>1:
        pyplot.figure(fig1.number)
        L = ax1[1].legend(ncol=2,loc=9,bbox_to_anchor=(0.5,-0.2),frameon=False,framealpha=0)
        for legobj in L.legendHandles: legobj.set_linewidth(3.0)
        pyplot.savefig(args.cmbimg+'--time.'+fmt)
        pyplot.close()
        pyplot.figure(fig2.number)
        L = ax2[1].legend(loc=2,bbox_to_anchor=(1,1),frameon=False,framealpha=0)
        for legobj in L.legendHandles: legobj.set_linewidth(3.0)
        pyplot.savefig(args.cmbimg+'--stress_strain.'+fmt)
        pyplot.close()
