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

# Fit data to a curve.
import scipy.optimize as opt
import numpy as np
from numpy.linalg import *
from models import *


#--------------------------------------------------------------------------------
# Calculate Rsquared coefficient of determination.
#--------------------------------------------------------------------------------
def rsquared(ydata,yfit):
    # Sum of the squares of the residuals.
    SSres = np.sum(np.power(ydata-yfit,2))
    # Sum of the squares.
    yavg = np.average(ydata)
    SStot = norm(ydata-yavg)
    return 1.0 - (SSres/SStot)


#--------------------------------------------------------------------------------
# Fit data to a curve with one function.
#--------------------------------------------------------------------------------
def curve_fit1(F,xdata,ydata,method,p0,maxitr):
    min_fun = lambda args: norm(F(xdata,*args)-ydata)
    #  OR = opt.minimize(min_fun,p0,method=method,options={'maxiter':maxitr,'disp':True})
    OR = opt.minimize(min_fun,p0,method=method,options={'maxiter':maxitr})
    if not OR.success: print '\t\tWARNING: '+OR.message
    return OR.x


#--------------------------------------------------------------------------------
# Fit data to a curve using basinhopping.
#--------------------------------------------------------------------------------
def curve_fit1_basinhopping(F,xdata,ydata,p0,nitr):
    min_fun = lambda args: norm(F(xdata,*args)-ydata)
    OR = opt.basinhopping(min_fun,p0,niter=nitr)
    return OR.x


#--------------------------------------------------------------------------------
# Fit data to a curve with two functions.
#--------------------------------------------------------------------------------
def curve_fit2(F1,xdata1,ydata1,F2,xdata2,ydata2,method,p0,maxitr):
    min_fun = lambda args: norm(F1(xdata1,*args)-ydata1) + \
                           norm(F2(xdata2,*args)-ydata2)
    #  OR = opt.minimize(min_fun,p0,method=method,options={'maxiter':maxitr,'disp':True})
    OR = opt.minimize(min_fun,p0,method=method,options={'maxiter':maxitr})
    if not OR.success: print '\t\tWARNING: '+OR.message
    return OR.x


#--------------------------------------------------------------------------------
# Fit data to a curve with three functions.
#--------------------------------------------------------------------------------
def curve_fit3(F1,xdata1,ydata1,F2,xdata2,ydata2,F3,xdata3,ydata3,method,p0,maxitr):
    min_fun = lambda args: norm(F1(xdata1,*args)-ydata1) + \
                           norm(F2(xdata2,*args)-ydata2) + \
                           norm(F3(xdata3,*args)-ydata3)
    #  OR = opt.minimize(min_fun,p0,method=method,options={'maxiter':maxitr,'disp':True})
    OR = opt.minimize(min_fun,p0,method=method,options={'maxiter':maxitr})
    if not OR.success: print '\t\tWARNING: '+OR.message
    return OR.x

