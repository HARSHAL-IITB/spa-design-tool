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
import os.path as path
import model_utility as MU
import calc_hyperelastic_parameters as CHP
import argparse

#--------------------------------------------------------------------------------
# Main.
#--------------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description="Create a plot from hyperelastic data and/or fits.",
            epilog="Example: ./create_hyperelastic_plots.py --xlim 0 3 plot.eps ogden3--UniaxialPlanar.mat E30--uni11.dat none E30--planar6b.dat")
    parser.add_argument("img",help="Output image to create (eps or png).")
    parser.add_argument("matfile",help="Material file defining model fit.")
    parser.add_argument("uniaxial",nargs='?',default='none',help="Optional uniaxial dataset to plot (default 'none').")
    parser.add_argument("biaxial",nargs='?',default='none',help="Optional biaxial dataset to plot (default 'none').")
    parser.add_argument("planar",nargs='?',default='none',help="Optional planar dataset to plot (default 'none').")
    parser.add_argument("--discard",type=float,help="Discard datapoints below a certain strain.")
    parser.add_argument("--xlim",nargs=2,type=float,default=[0.0,0.0],help="Min,Max values on xscale (AUTO if undefined).")
    parser.add_argument("--ylim",nargs=2,type=float,default=[0.0,0.0],help="Min,Max values on yscale (AUTO if undefined).")
    parser.add_argument("--ylimerr",nargs=2,type=float,default=[0.0,0.0],help="Min,Max values on yscale for error plots (AUTO if undefined).")
    parser.add_argument("--paper",action='store_true',help="Create plots designed for the paper, rather than for general use.")
    args = parser.parse_args()

    # Mathematical rendering can be enabled by enclosing the symbols in $'s
    pyplot.rc('savefig',dpi=300)
    pyplot.rc('font',size=(14 if args.paper else 8))
    pyplot.rc('mathtext',default='regular') # Don't use italics for mathmode.

    # Read in the material file.
    mat = MU.read_matfile(args.matfile)
    # Read in the data files.
    data = dict()
    MU.import_dataset(args.uniaxial,data,'u')
    MU.import_dataset(args.biaxial,data,'b')
    MU.import_dataset(args.planar,data,'p')
    # Discard datapoints if requested.
    if args.discard:
        for key in 'ubp':
            if not key in data: continue
            data[key] = data[key][data[key][:,0]>=args.discard,:]
            data[key] = np.vstack(([0.0,0.0],data[key]))
    # Get the model.
    mname = (mat['model']+str(mat['order'])) if 'order' in mat else mat['model']
    model = MU.get_model(mname)

    # Calculate the Rsquared quality of fit.
    R2 = CHP.calculate_rsquared(data,model,mat['params'],'ubp')
    CHP.create_plots(model,data,mat['params'],R2,mat['D'],args.img,mat['matname'],args.xlim,args.ylim,args.ylimerr,args.paper)
