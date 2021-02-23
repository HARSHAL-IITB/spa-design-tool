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
import argparse
import os.path as path
import numpy as np


#--------------------------------------------------------------------------------
# Main.
#--------------------------------------------------------------------------------
if __name__ == "__main__":
    # Handle user input.
    scriptname = path.basename(__file__)
    parser = argparse.ArgumentParser(
            description='Create a PCA-style plot showing scatter of optimization results.',
            epilog='Example: '+scriptname+' scatterplot.png 4 Displacement 8 Force 7 records4.csv --percent_err 6 5 70 6 --xlim 0 100 --ylim 0 100 --aspect_equal')
    parser.add_argument('img',help='Output image to create.')
    parser.add_argument('nerr',type=int,help='Number of error columns.')
    parser.add_argument('xvar',help='Variable on x-axis.')
    parser.add_argument('xcol',type=int,help='Column in dataset for x-axis.')
    parser.add_argument('yvar',help='Variable on y-axis.')
    parser.add_argument('ycol',type=int,help='Column in dataset for y-axis.')
    parser.add_argument('csv',nargs='+',help='Dataset(s) to plot, typically records.csv files from optimize_geometric_parameters.py runs.')
    parser.add_argument('--xlim',nargs=2,type=float,metavar=('MIN','MAX'),help='Optional min and max for x-axis.')
    parser.add_argument('--ylim',nargs=2,type=float,metavar=('MIN','MAX'),help='Optional min and max for y-axis.')
    parser.add_argument('--legend',action='store_true',help='Include a legend in the plot.')
    parser.add_argument('--aspect_equal',action='store_true',help='Set an euqal aspect-ratio for plot.')
    parser.add_argument('--percent_err',nargs='+',type=float,
            help='Convert R2 values to percent errors. Values here are the \'goal\' values for each error column.')
    args = parser.parse_args()
    np.set_printoptions(precision=4,linewidth=130,suppress=True) # Disable scientific notation for small numbers.
    pyplot.rc('mathtext',default='regular') # Don't use italics for mathmode.
    # pyplot.rc('savefig',dpi=300)
    # pyplot.rc('font',size=8)

    # Read in the given datasets.
    data = dict()
    minimizer_itr = 35      # Recently the optimize_geometric_parameters.py script uses 35. Before it used 40.
    num_basins_total = 0
    for fname in args.csv:
        fdata = np.loadtxt(fname,comments='#',delimiter=',')
        data[fname] = fdata
        num_basins_total = num_basins_total + int(np.ceil(max(fdata[:,0]/minimizer_itr)))
        print 'Imported',fdata.shape[0],'rows from',fname

    # Convert R2 values to real-world values, if requested.
    print '----------------------------------------------------------------------------------'
    if args.percent_err:
        if args.nerr != len(args.percent_err):
            print 'ERROR: number of arguments to --percent_err must be equal to the NERR parameter.'
            exit(1)
        print 'Errors are calculated as: err = 100.0 * abs(sim-goal)/goal == percent error'
        for fname in data.keys():
            for col,goal in enumerate(args.percent_err):
                SStot = 2.0 * ((args.percent_err[col]/2.0)**2.0)
                dcol = data[fname].shape[1]-args.nerr+col
                data[fname][:,dcol] = 100.0 * np.sqrt(SStot*data[fname][:,dcol]) / args.percent_err[col]
    else:
        print 'Errors are calculated as: err = sum((sim-goal)^2) / sum((goal-goal_avg)^2) == R2'
    print 'Changing between percent error and R2 error changes the order of results slightly,'
    print 'but typically the top geometries remain the same to the first few digits.'
    print '----------------------------------------------------------------------------------'

    # Prepare the plots.
    fig,axis = pyplot.subplots(1,1)
    # cm = pyplot.get_cmap('gist_rainbow')
    # cm = pyplot.get_cmap('gist_ncar')
    # cm = pyplot.get_cmap('jet')
    cm = pyplot.get_cmap('hsv')

    # Given the axes limits, determine how many colors we need to display the points.
    xmax = args.xlim[1] if args.xlim else 1e100
    ymax = args.ylim[1] if args.ylim else 1e100
    if args.xlim or args.ylim:
        num_basins_total = 0
        for fdata in data.itervalues():
            num_basins = int(np.ceil(max(fdata[:,0]/minimizer_itr)))
            for b in range(num_basins):
                bmin = b*minimizer_itr
                bmax = b*minimizer_itr+minimizer_itr
                basin_pts = (bmin<=fdata[:,0]) & (fdata[:,0]<bmax)
                if sum(basin_pts)==0: continue
                if min(fdata[basin_pts,args.xcol])<xmax or min(fdata[basin_pts,args.ycol])<ymax:
                    num_basins_total = num_basins_total+1
    colors = [cm(float(i)/num_basins_total) for i in range(num_basins_total)]

    # Create scatterplot.
    cidx = 0
    for fname,fdata in data.iteritems():
        num_basins = int(np.ceil(max(fdata[:,0]/minimizer_itr)))
        for b in range(num_basins):
            bmin = b*minimizer_itr
            bmax = b*minimizer_itr+minimizer_itr
            basin_pts = (bmin<=fdata[:,0]) & (fdata[:,0]<bmax)
            print '  Points in basin '+str(b+1)+' ('+fname+'):',sum(basin_pts)
            if sum(basin_pts)==0: continue
            if min(fdata[basin_pts,args.xcol])<xmax or min(fdata[basin_pts,args.ycol])<ymax:
                axis.scatter(fdata[basin_pts,args.xcol],fdata[basin_pts,args.ycol],
                             color=colors[cidx],alpha=0.5,s=40,
                             label='Basin '+str(b+1)+' ('+fname+')')
                cidx = cidx + 1
            else: print '    Not plotting basin '+str(b+1)+', no points within given limits.'
    if args.percent_err:
        axis.set_xlabel('Absolute Percent Difference in '+args.xvar+' (%)')
        axis.set_ylabel('Absolute Percent Difference in '+args.yvar+' (%)')
    else:
        axis.set_xlabel(r'$R^2$ Difference in '+args.xvar)
        axis.set_ylabel(r'$R^2$ Difference in '+args.yvar)
    if args.xlim: axis.set_xlim(args.xlim)
    if args.ylim: axis.set_ylim(args.ylim)
    if args.aspect_equal: axis.set_aspect('equal')
    axis.grid()

    if args.legend:
        # axis.legend(loc='best',frameon=False,framealpha=0)
        lgd = axis.legend(loc=2,bbox_to_anchor=(1,1),frameon=False,framealpha=0)
        pyplot.savefig(args.img,bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
        pyplot.tight_layout()
        pyplot.savefig(args.img)

    # Calculate the total error at each point.
    # Sort the array by error value, and print the best values.
    print 'Best result from each basin:'
    for fname,fdata in data.iteritems():
        num_basins = int(np.ceil(max(fdata[:,0]/minimizer_itr)))
        basin_best = np.zeros((num_basins,fdata.shape[1]))
        basin_best_errs = []
        err = np.sum(fdata[:,-args.nerr:],axis=1)
        for b in range(num_basins):
            bmin = b*minimizer_itr
            bmax = b*minimizer_itr+minimizer_itr
            basin_pts = (bmin<=fdata[:,0]) & (fdata[:,0]<bmax)
            if sum(basin_pts)==0:
                basin_best_errs.append(1.e5)
            else:
                basin_data = fdata[basin_pts,:]
                basin_errs = err[basin_pts]
                basin_data = basin_data[np.argsort(basin_errs),:]
                basin_best[b,:] = basin_data[0,:]
                basin_best_errs.append(min(basin_errs))
        basin_best = basin_best[np.argsort(basin_best_errs),:]
        basin_best_errs.sort()  # This sort needs to come after argsort, or we lose the indexing.
        print '  '+fname+', total error for each shown basin:'
        print '    '+str(basin_best_errs)
        print basin_best

    # pyplot.show()     # Interactive exploration of plot data.
    pyplot.close()
