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

# Use this script to extract displacements or blocked forces from ABQ results.
# Note that you probably never want to run this script directly;
# instead prefer to use the run_tests.py wrapper script.
# Run using: abaqus cae noGUI=abq_extract_data.py -- --help
#            Output will be written to the abaqus.rpy file.
# Note that if you're editing this file, you may want to run the following command in ABQ
# in order to get more useful history in the rpy file:
#   session.journalOptions.setValues(replayGeometry=COORDINATE,recoverGeometry=COORDINATE)

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
import sys, pprint, os
sys.path.append('/usr/lib/python2.7/')
import argparse
import inspect, os
sys.path.append(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
import model_utility as MU


#--------------------------------------------------------------------------------
# Main.
#--------------------------------------------------------------------------------
if __name__ == "__main__":
    # Handle user input. We can't do any optional arguments, because we'd have
    # to find them in the sys.argv, which includes a bunch of Abaqus junk.
    # Also, Abaqus won't let us pass "--xxx=xxx" style arguments anyway.
    NARGS = 2
    run_cmd = 'abaqus cae noGUI=abq_create_geom.py -- '
    if '--help' in sys.argv: sys.argv = [run_cmd,'--help']
    else: sys.argv = [run_cmd]+sys.argv[len(sys.argv)-NARGS:]
    parser = argparse.ArgumentParser(
            description='Extract data from an Abaqus simulation.',
            epilog='Example: '+run_cmd+'job.odb data.rpt')
    parser.add_argument('odb',help='Abaqus results database.')
    parser.add_argument('rpt',help='Filename for extracted data.')
    args = parser.parse_args()
    if not args.odb.endswith('.odb'):
        print 'ERROR: incorrect inputs, need an ODB file. Found arguments:'
        print sys.argv
        exit(1)

    # Open the ODB file.
    o1 = session.openOdb(name=args.odb)
    session.viewports['Viewport: 1'].setValues(displayedObject=o1)
    odb = session.odbs[args.odb]

    # Determine the type of problem.
    try:
        t = odb.rootAssembly.instances['SPA'].nodeSets['SYMM_SHELL']
        actuator = 'bnd'
        scale_factor = 2.0
    except:
        try:
            t = odb.rootAssembly.instances['SPA'].nodeSets['SYMM_T']
            actuator = 'bubble'
            scale_factor = 2.0
        except:
            actuator = 'lin'
            scale_factor = 4.0
    print 'Actuator =',actuator
    try:
        t = odb.rootAssembly.nodeSets['FPT']
        test = 'F'
    except:
        try:
            t = odb.rootAssembly.instances['SPA'].nodeSets['UPT']
            # t = odb.rootAssembly.nodeSets['UPT']
            test = 'U'
        except:
            try:
                t = odb.rootAssembly.nodeSets['EPT']
                test = 'visco'
            except:
                test = 'unknown'
    print 'Test =',test

    # Create XY Data.
    if test=='F':
        session.xyDataListFromField(odb=odb, outputPosition=ELEMENT_FACE,
                variable=(('P', ELEMENT_FACE), ), elementSets=('SPA.PRESSURE', ))
        x0 = session.xyDataObjects[session.xyDataObjects.keys()[0]]
        session.XYDataFromHistory(name='RF1', odb=odb, steps=('Step-1', ),
            outputVariableName='Reaction force: RF1 PI: rootAssembly Node 2 in NSET FPT')
        session.XYDataFromHistory(name='RF2', odb=odb, steps=('Step-1', ),
            outputVariableName='Reaction force: RF2 PI: rootAssembly Node 2 in NSET FPT')
        session.XYDataFromHistory(name='RF3', odb=odb, steps=('Step-1', ),
            outputVariableName='Reaction force: RF3 PI: rootAssembly Node 2 in NSET FPT')
        session.XYDataFromHistory(name='ZMIN', odb=odb, steps=('Step-1', ),
            outputVariableName='Coordinates: COOR3 PI: rootAssembly Node 2 in NSET FPT')
        session.XYDataFromHistory(name='ZMAX', odb=odb, steps=('Step-1', ),
            outputVariableName='Coordinates: COOR3 PI: rootAssembly Node 1 in NSET LPT')
        x1 = session.xyDataObjects['RF1']
        x2 = session.xyDataObjects['RF2']
        x3 = session.xyDataObjects['RF3']
        x4 = session.xyDataObjects['ZMIN']
        x5 = session.xyDataObjects['ZMAX']
        # Scale data, due to symmetry.
        x1S = []
        x2S = []
        x3S = []
        for i,(x,y) in enumerate(x1): x1S.append((x,scale_factor*y))
        for i,(x,y) in enumerate(x2): x2S.append((x,scale_factor*y))
        for i,(x,y) in enumerate(x3): x3S.append((x,scale_factor*y))
        session.xyDataObjects['RF1'].setValues(data=x1S)
        session.xyDataObjects['RF2'].setValues(data=x2S)
        session.xyDataObjects['RF3'].setValues(data=x3S)
        data = (x0, x1, x2, x3, x4, x5)
    elif test=='U':
        session.xyDataListFromField(odb=odb, outputPosition=ELEMENT_FACE,
                variable=(('P', ELEMENT_FACE), ), elementSets=('SPA.PRESSURE', ))
        x0 = session.xyDataObjects[session.xyDataObjects.keys()[0]]
        for i in range(1,100):
            try:
                if actuator=='lin' or actuator=='bnd':
                    ovn1='Spatial displacement: U1 PI: SPA Node '+str(i)+' in NSET UPT'
                    ovn2='Spatial displacement: U2 PI: SPA Node '+str(i)+' in NSET UPT'
                    ovn3='Spatial displacement: U3 PI: SPA Node '+str(i)+' in NSET UPT'
                    ovn4='Coordinates: COOR3 PI: SPA Node '+str(i)+' in NSET UPT'
                    ovn5='Coordinates: COOR3 PI: rootAssembly Node 1 in NSET LPT'
                    session.XYDataFromHistory(name='ZMIN', odb=odb, steps=('Step-1', ), outputVariableName=ovn4)
                    session.XYDataFromHistory(name='ZMAX', odb=odb, steps=('Step-1', ), outputVariableName=ovn5)
                elif actuator=='bubble':
                    ovn1='Spatial displacement: U1 at Node '+str(i)+' in NSET UPT'
                    ovn2='Spatial displacement: U2 at Node '+str(i)+' in NSET UPT'
                    ovn3='Spatial displacement: U3 at Node '+str(i)+' in NSET UPT'
                session.XYDataFromHistory(name='U1', odb=odb, steps=('Step-1', ), outputVariableName=ovn1)
                session.XYDataFromHistory(name='U2', odb=odb, steps=('Step-1', ), outputVariableName=ovn2)
                session.XYDataFromHistory(name='U3', odb=odb, steps=('Step-1', ), outputVariableName=ovn3)
                break
            except: continue
        x1 = session.xyDataObjects['U1']
        x2 = session.xyDataObjects['U2']
        x3 = session.xyDataObjects['U3']
        if actuator=='lin' or actuator=='bnd':
            x4 = session.xyDataObjects['ZMIN']
            x5 = session.xyDataObjects['ZMAX']
            data = (x0, x1, x2, x3, x4, x5)
        elif actuator=='bubble':
            data = (x0, x1, x2, x3)
    elif test=='visco':
        session.XYDataFromHistory(name='RF', odb=odb, steps=('Step-1', ),
            outputVariableName='Reaction force: RF3 PI: rootAssembly Node 1 in NSET EPT')
        session.XYDataFromHistory(name='U3', odb=odb, steps=('Step-1', ),
            outputVariableName='Spatial displacement: U3 PI: rootAssembly Node 1 in NSET EPT')
        x0 = session.xyDataObjects['RF']
        x1 = session.xyDataObjects['U3']
        data = (x0, x1)
    else:
        print 'ERROR: unknown test type.'
        sys.exit(1)

    # Write the data to disk.
    print 'Writing data to disk at',args.rpt
    session.xyReportOptions.setValues(numDigits=8, numberFormat=AUTOMATIC)
    session.writeXYReport(fileName=args.rpt, appendMode=OFF, xyData=data)

    # Try to extract energies.
    try:
        session.XYDataFromHistory(name='IE', odb=odb, steps=('Step-1', ),
                outputVariableName='Internal energy: ALLIE for Whole Model')
        session.XYDataFromHistory(name='KE', odb=odb, steps=( 'Step-1', ),
                outputVariableName='Kinetic energy: ALLKE for Whole Model')
        ie = session.xyDataObjects['IE']
        ke = session.xyDataObjects['KE']
        data = (ie, ke)
        fname = os.path.splitext(args.rpt)
        fname = fname[0]+'-energy'+fname[1]
        print 'Writing energy data to disk at',fname
        session.writeXYReport(fileName=fname, appendMode=OFF, xyData=data)
    except:
        print 'No energy data found.'

    # Try to extract strain data.
    # NOTE: there may be a difference here between the numbers extracted and the numbers shown
    # in Abaqus contour plots. For the contour plots Abaqus seems to be recalculating the values
    # at the nodes, whereas the values are initially calculated at integration points. If you
    # use the 'probe values' tool in the Abaqus GUI, you'll find the same values as extracted here.
    try:
        time = []
        vals = []
        for step in odb.steps.values():
            for frame in step.frames:
                if frame.fieldOutputs.has_key('NE'):
                    maxval = 0.0
                    for val in frame.fieldOutputs['NE'].values:
                        maxval = max(val.maxPrincipal,maxval)
                    time.append(frame.frameValue)
                    vals.append(maxval)
        fname = os.path.splitext(args.rpt)
        fname = fname[0]+'-strain.csv'
        print 'Writing strain data to disk at',fname
        with open(fname,'w') as of:
            of.write('# simulation time, nominal (engineering) strain\n')
            for i in range(len(time)):
                of.write(str(time[i])+','+str(vals[i])+'\n')
    except:
        print 'No strain data found.'


    o1.close() # Close the output database before exiting.
