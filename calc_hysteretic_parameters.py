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

import argparse, time, sys, os
import os.path as path
import numpy as np
import scipy.optimize as opt
import utility as U
import run_tests as RT
from matplotlib import pyplot,axes


#--------------------------------------------------------------------------------
# Minimization function.
#--------------------------------------------------------------------------------
def fit_error(data,cmd,A,L,*popt):
    print '\n********************************************************************'
    print '** Testing with parameters: ',popt
    print '********************************************************************'
    # Get the simulation results.
    if not hasattr(fit_error, "counter"):
        fit_error.counter = 0
    vdata = run_viscotest(cmd,A,L,fit_error.counter,*popt)
    fit_error.counter = fit_error.counter + 1
    if vdata==[]:
        print '\n********************************************************************'
        print ' Abaqus failed, we return an average error of 1000.0'
        print '********************************************************************'
        return 1000.0
    # Find the comparable experimental results.
    edata = np.zeros(vdata.shape)
    for i in range(vdata.shape[0]):
        vstrain = vdata[i,1]
        for j in range(data.shape[0]):
            if data[j,1] >= vdata[i,1]:
                edata[i,:] = data[j,:]
                break
    # Make the comparison.
    err = np.linalg.norm(edata[:,0]-vdata[:,0])
    # Plot results.
    fig,ax = pyplot.subplots()
    ax.plot(data[:,1],data[:,0],label='Experimental Data')
    ax.plot(vdata[:,1],vdata[:,0],'o-',label='Simulation Data')
    ax.set_title(str(popt)+'\nitr='+str(fit_error.counter)+', err='+str(err))
    ax.set_xlabel('Engineering Strain')
    ax.set_ylabel('Nominal Stress')
    ax.grid()
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.set_aspect(abs(xlim[1]-xlim[0])/abs(ylim[1]-ylim[0])) # Square axis.
    ax.legend(loc='best',frameon=False,framealpha=0)
    pyplot.tight_layout()
    pyplot.savefig('stress_strain--'+str(fit_error.counter)+'.png')
    pyplot.savefig('stress_strain--latest.png')
    pyplot.close()
    print '\n********************************************************************'
    print ' Calculated sum-of-squares error:',err
    print '********************************************************************'
    return err

#--------------------------------------------------------------------------------
# Create and run an Abaqus viscoelastic test.
#   cmd     = Abaqus command to create test.
#   A       = cross-sectional area.
#   L       = sample length in testing direction.
#   i       = iteration number.
#   args    = S,A,m,C visco parameters.
#--------------------------------------------------------------------------------
def run_viscotest(cmd,A,L,i,*args):
    PROCS = 4
    cmd = cmd+map(str,args)
    U.run_cmd(cmd)
    if RT.run_abq('test_visco.inp','.',PROCS,str(i)):
        if RT.run_post('test_visco.inp','.',str(i)):
            results = np.loadtxt(path.join('test_visco-'+str(i),'data.rpt'),skiprows=3)
            results = results[:,1:]
            results[:,0] = results[:,0] / A         # Calculate nominal stress.
            results[:,1] = results[:,1] / (0.5*L)   # Calculate engineering strain (with symmetry).
            return results
    return []

#--------------------------------------------------------------------------------
# Main.
#--------------------------------------------------------------------------------
if __name__ == "__main__":
    tinit = time.time()
    # Handle user input.
    parser = argparse.ArgumentParser(description="Fit the given strain-rate datasets to hysteresis parameters.",
                                     epilog="Example: calc_hysteretic_parameters.py")
    parser.add_argument("-m","--method",choices=['l-bfgs-b','slsqp'],default='l-bfgs-b',help="Method to use for fitting.")
    # parser.add_argument("--format",choices=['png','eps'],help="Image format, default is png.")
    parser.add_argument("dirname",help="Name of directory to create with all the tests in it.")
    parser.add_argument("dims",nargs=3,help="Depth, width, and length of sample (mm). Length will be divided by 2 for symmetry.")
    parser.add_argument('rate',type=float,help='Strain rate (mm/s).')
    parser.add_argument("datafile",help="Dataset to fit, consisting of (stress, strain) columns.")
    parser.add_argument("matfile",help="File containing hyperelastic material properties.")
    args = parser.parse_args()
    # 'suppress' disables scientific notation for small numbers.
    np.set_printoptions(precision=4,linewidth=130,suppress=True)
    # np.seterr(all='raise')
    pyplot.rc('savefig',dpi=300)
    pyplot.rc('font',size=8)
    pyplot.rc('mathtext',default='regular') # Don't use italics for mathmode.

    # Read in the given datasets.
    print '--------------------------------------------------------------------'
    print ' Importing dataset...'
    data = np.loadtxt(args.datafile,comments='#',delimiter=',',dtype=np.longdouble)
    print '  Imported',data.shape[0],'datapoints.'
    # TODO - smooth dataset?

    # Prep variables.
    matfile = path.abspath(args.matfile)
    max_strain = max(data[:,1])
    A = float(args.dims[0]) * float(args.dims[1])
    L = float(args.dims[2])     # Divided by two later, in abq_hysttest.py
    time = L*max_strain / args.rate

    # Change directory.
    if path.exists(args.dirname): U.print_error('Output directory exists.',True)
    os.makedirs(args.dirname)
    rootdir = os.getcwd()
    os.chdir(args.dirname)
    script = U.find_file('python/abq_hysttest.py')

    # Calculate optimal parameters.
    print '---------------------------------------------------------'
    print ' Calculating parameters...'
    cmd0 = ['abaqus','cae','noGUI='+script,'--',matfile,
            args.dims[0],args.dims[1],args.dims[2],
            str(max_strain),str(time)]
    maxitr = 1000
    popt = [1.6, 0.5, 4.0, -1.0]
    min_fun = lambda args: fit_error(data,cmd0,A,L,*args)
    OR = opt.minimize(min_fun,popt,method=args.method,options={'maxiter':maxitr,'eps':0.05},
            bounds=((0,None),(0,None),(0,None),(-1,0)))
    if not OR.success: print '\t\tWARNING: '+OR.message
    popt = OR.x
    # S = calculate_rsquared(data,poisson,popt)
    # print '\t\tRsquared: ',S


    print '\n\n--------------------------------------------------------------------'
    print ' Results for fit to strain-rate data.'
    print '--------------------------------------------------------------------'
    np.set_printoptions(suppress=False)
    # print 'Rsquared:',R2
    print 'Hysteresis parameters:'
    print 'S: ',popt[0]
    print 'A: ',popt[1]
    print 'm: ',popt[2]
    print 'C: ',popt[3]
    os.chdir(rootdir)
    print 'TOTAL TIME ELAPSED: ',U.time_elapsed(tinit)
