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

import argparse, warnings
import os.path as path
import numpy as np
import utility
from matplotlib import pyplot,axes


#--------------------------------------------------------------------------------
# Calculate the time from the beginning of the step for each step.
#   data   = experimental datapoints (time,length,force).
#   moduli = output moduli (time,shear,bulk)
#--------------------------------------------------------------------------------
def calc_steptime(data,moduli):
    # Step ends when length changes again after 10 constant datapoints.
    tstart = data[0,0]
    slen = data[0,1]
    nc = 0
    stepS = [0]
    stepE = []
    # Find the start and end of rise of each step.
    for i in range(data.shape[0]):
        if np.abs(data[i,1]-slen)<0.0001:
            nc = nc + 1
            if nc==10:
                # stepE.append(i-10)
                stepE.append(i)     # TODO - this crops 10 points at the beginning.
                tstart = data[i,0]
        else:
            if nc>10:
                # tstart = data[i,0]
                stepS.append(i)
            slen = data[i,1]
            nc = 0
        data[i,0] = data[i,0] - tstart
        moduli[i,0] = moduli[i,0] - tstart
    stepS.append(data.shape[0])

    # Remove rising data from moduli.
    for i in reversed(range(len(stepE))):
        moduli = np.delete(moduli,range(stepS[i]+1,stepE[i]),axis=0)
        moduli[stepS[i],:] = [np.nan,np.nan,np.nan]
    return data,moduli,stepS,stepE


#--------------------------------------------------------------------------------
# Calculate the moduli from force, area, length, and Poisson's ratio.
#   data       = experimental datapoints (time,length,force).
#   moduli     = output moduli (time,shear,bulk)
#   planarBOOL = True if planar data.
#   A          = applied-force area.
#   nu         = Poisson's ratio.
#--------------------------------------------------------------------------------
def calc_moduli(data,moduli,planarBOOL,A,nu):
    print 'Cross-sectional area:     ',A
    min_idx = np.argmin(np.abs(data),axis=0)    # Index of minimum force value.
    z = data[min_idx[2],1]-0.001                # Length at minimum force, -delta avoids /0
    print 'Minimum (abs val) force:  ',data[min_idx[2],2]
    print 'Length at minimum force:  ',z
    data[:,1] = (data[:,1] - z)/z                   # Engineering strain.
    # data[:,1] = np.log(data[:,1]/z)               # Logarithmic strain.
    data[:,2] = data[:,2]/A                         # Nominal stress.
    E = data[:,2]/data[:,1]                         # Young's moduli (E).
    if planarBOOL: E = E*(1.0-nu*nu)
    moduli[:,1] = E/(2.0*(1.0+nu))                  # Shear moduli (G)
    if nu<0.5: moduli[:,2] = E/(3.0*(1.0-2.0*nu))   # Bulk moduli (G)
    return data,moduli


#--------------------------------------------------------------------------------
# Plot the force or stress versus time.
#   ax     = axes from figure to plot on.
#   data   = experimental datapoints (time,strain,stress).
#   moduli = output moduli (time,shear,bulk)
#   stepS  = list of starting rows for steps.
#   plabel = label for plot.
#--------------------------------------------------------------------------------
def plot_vs_time(ax,data,moduli,stepS,plabel):
    # Plot strain and stress.
    for i in range(len(stepS)-1):
        ax[0].plot(data[stepS[i]:stepS[i+1],0],data[stepS[i]:stepS[i+1],1],'b-',label=plabel)
    for i in range(len(stepS)-1):
        ax[1].plot(data[stepS[i]:stepS[i+1],0],data[stepS[i]:stepS[i+1],2],'r-',label=plabel)
    ax[0].set_xscale('log')
    ax[0].set_xlabel('Time (s)')
    ax[0].set_ylabel('Engineering Strain', color='b')
    ax[1].set_xscale('log')
    ax[1].set_xlabel('Time (s)')
    ax[1].set_ylabel(r'Nominal Stress $\left({}^N\!/{}_{mm^2}\right)$', color='r')
    for t in ax[0].get_yticklabels(): t.set_color('b')
    for t in ax[1].get_yticklabels(): t.set_color('r')
    # Plot shear moduli (bulk is the same, just scaled).
    ax[2].plot(moduli[:,0],moduli[:,1],label=plabel)
    ax[2].set_xscale('log')
    ax[2].set_xlabel('Time (s)')
    ax[2].set_ylabel(r'Shear Moduli $\left({}^N\!/{}_{mm^2}\right)$')
    pyplot.tight_layout()


#--------------------------------------------------------------------------------
# Plot the the curves as stress-strain.
#   ax     = axes from figure to plot on.
#   data   = experimental datapoints (time,strain,stress).
#   stepS  = list of starting rows for steps (where strain starts to rise).
#   stepE  = list of ending rows for steps (where strain stops rising).
#--------------------------------------------------------------------------------
def plot_stress_strain(ax,data,stepS,stepE):
    strain_inst = data[stepE,1]
    stress_inst = data[stepE,2]
    eq = np.concatenate((stepS[:-1],[stepS[-1]-1]))
    strain_eq = data[eq,1]
    stress_eq = data[eq,2]
    # ax.plot(strain_inst,stress_inst,'r.-',label='Immediate')
    # ax.plot(strain_eq,stress_eq,'b.-',label='Relaxed')
    ax.plot(data[:,1],data[:,2],'g',label='Time-Series')
    ax.set_xlabel('Engineering Strain')
    ax.set_ylabel(r'Nominal Stress $\left({}^N\!/{}_{mm^2}\right)$')
    ax.legend(loc='best',frameon=False,framealpha=0)
    pyplot.tight_layout()


#--------------------------------------------------------------------------------
# Main.
#--------------------------------------------------------------------------------
if __name__ == "__main__":
    # Handle user input.
    parser = argparse.ArgumentParser(description="Calculate the shear and bulk relaxation moduli versus time for a given dataset.",
            epilog="Example: calc_viscomoduli.py 100 test1.dat test2.dat")
    parser.add_argument("-f","--cmbimg",metavar="FNAME",default="combined",help="Filename for combined plots (default \'combined\').")
    parser.add_argument("-t","--title",help="Title for combined plots.")
    parser.add_argument("--format",choices=['png','eps'],help="Image format, default is png.")
    parser.add_argument("--planar",action='store_true',help="Data is from a planar test, rather than a uniaxial test.")
    parser.add_argument("area",type=float,help="Area over which force was applied.")
    parser.add_argument("nu",type=float,help="Poisson's ratio for this material.")
    parser.add_argument("datafiles",nargs='+',help="Dataset(s) to process.")
    args = parser.parse_args()
    if args.format: fmt = '.'+args.format
    else:           fmt = '.png'
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
        print 'WARNING - multiplots are broken.'    # TODO
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
            dataO = np.loadtxt(datafile,usecols=(0,1,2),dtype=np.longdouble)
        except:
            utility.print_error('Failed to read file.',False)
            print '\n'
            continue
        print 'Number of points imported:',dataO.shape[0]
        fname = path.splitext(datafile)
        data = np.zeros(dataO.shape)
        data[:,0] = dataO[:,0]      # Time (s).
        data[:,1] = dataO[:,1]      # Total length (mm).
        data[:,2] = dataO[:,2]      # Force (N).
        moduli = np.zeros(data.shape)
        moduli[:,0] = data[:,0]

        # Calculate the moduli.
        data,moduli = calc_moduli(data,moduli,args.planar,args.area,args.nu)
        # Break up the time into step-time.
        data,moduli,stepS,stepE = calc_steptime(data,moduli)

        # Plot the stress, strain, and moduli versus time.
        print 'Plotting data...'
        if len(datafiles)>1:
            pyplot.figure(fig1.number)
            plot_vs_time(ax1,data,moduli,stepS,fname[0])
        ax = []
        fig = pyplot.figure(figsize=(15,5))
        ax.append(fig.add_subplot(141))
        ax.append(ax[0].twinx())
        ax.append(fig.add_subplot(142))
        ax.append(fig.add_subplot(1,4,(3,4)))
        plot_vs_time(ax,data,moduli,stepS,'none')
        plot_stress_strain(ax[-1],data,stepS,stepE)
        # Catch a warning about less_equal on NaN.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pyplot.savefig(fname[0]+fmt)
        pyplot.close()

        # Write the new dataset to disk.
        outfile = fname[0]+'--CLEAN'+fname[1]
        head = 'Time (s), Shear Moduli (N/mm2), Bulk Moduli (N/mm2) from '+datafile
        np.savetxt(outfile,moduli,delimiter=',',header=head)
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
