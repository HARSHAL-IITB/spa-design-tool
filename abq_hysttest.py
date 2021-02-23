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

# Use this script to create Abaqus geometry for a viscoelasticity test.
# Note that you probably never want to run this script directly;
# instead prefer to use the optimize_hyperelastic_properties.py wrapper script.
# Run using: abaqus cae noGUI=abq_viscotest.py -- --help
#            Output will be written to the abaqus.rpy file.
# Note that if you're editing this file, you may want to run the following command in ABQ
# in order to get more useful history in the rpy file:
#   session.journalOptions.setValues(replayGeometry=COORDINATE,recoverGeometry=COORDINATE)

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
import sys, pprint, os, shutil
sys.path.append('/usr/lib/python2.7/')
import argparse
import inspect, os
sys.path.append(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
import model_utility as MU


#--------------------------------------------------------------------------------
# Return the model given a string model name.
#--------------------------------------------------------------------------------
def get_abq_model(mname):
    mn = mname.lower()
    if mn=='ab':    return ARRUDA_BOYCE
    if mn=='mr':    return MOONEY_RIVLIN
    if mn=='neoh':  return NEO_HOOKE
    if mn=='ogden': return OGDEN
    if mn=='poly':  return POLYNOMIAL
    if mn=='rpoly': return REDUCED_POLYNOMIAL
    if mn=='vdw':   return VAN_DER_WAALS
    if mn=='yeoh':  return YEOH
    print 'ERROR: invalid model name \''+mname+'\''
    exit(1)

#--------------------------------------------------------------------------------
# Main.
#--------------------------------------------------------------------------------
if __name__ == "__main__":
    # Handle user input. We can't do any optional arguments, because we'd have
    # to find them in the sys.argv, which includes a bunch of Abaqus junk.
    # Also, Abaqus won't let us pass "--xxx=xxx" style arguments anyway.
    NARGS = 10
    run_cmd = 'abaqus cae noGUI=abq_viscotest.py -- '
    if '--help' in sys.argv: sys.argv = [run_cmd,'--help']
    else: sys.argv = [run_cmd]+sys.argv[len(sys.argv)-NARGS:]
    parser = argparse.ArgumentParser(
            description='Create and run an Abaqus test to calibrate hysteresis parameters.')
    parser.add_argument('matfile',help='File containing hyperelastic material properties.')
    parser.add_argument('D',type=float,help='Depth of sample.')
    parser.add_argument('W',type=float,help='Width of sample.')
    parser.add_argument('L',type=float,help='Length of sample (will be divided by 2 for symmetry).')
    parser.add_argument('strain',type=float,help='Maximum strain (not percent) to be applied.')
    parser.add_argument('time',type=float,help='Time period over which loading takes place (s).')
    parser.add_argument('S',type=float,help='Viscoelastic stress scaling factor.')
    parser.add_argument('A',type=float,help='Viscoelastic creep parameter.')
    parser.add_argument('m',type=float,help='Viscoelastic effective stress exponent.')
    parser.add_argument('C',type=float,help='Viscoelastic creep strain exponent.')
    args = parser.parse_args()

    #--------------------------------------------------------------------------------
    # Read the material file.
    #--------------------------------------------------------------------------------
    mat = MU.read_matfile(args.matfile)

    #--------------------------------------------------------------------------------
    # Calculate dimensions, assuming mirror-symmetry in z.
    #   depth  = Y dimension.
    #   width  = X dimension.
    #   length = Z dimension.
    #--------------------------------------------------------------------------------
    sample_D = args.D           # Depth of sample.
    sample_W = args.W           # Width of sample.
    sample_L = args.L / 2.0     # length of sample.
    displacement = args.strain*sample_L
    time = args.time


    #--------------------------------------------------------------------------------
    # Open the viewport.
    #--------------------------------------------------------------------------------
    session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=300, height=200)
    executeOnCaeStartup()
    Mdb()

    #--------------------------------------------------------------------------------
    # Create the part.
    #--------------------------------------------------------------------------------
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',sheetSize=200.0)
    s.setPrimaryObject(option=STANDALONE)

    # Extrude a solid rectangle for body (assuming mirror-symmetry in x).
    #
    #     STRAIN BC
    #     +-------+
    #    /       /|
    #   +-------+ |
    #   |       | |         Z (length)
    #   |       | |         ^
    #   |       | |         |     Y (depth)
    #   |       | +         |   /
    #   |       |/          | /
    #   +-------+           O-------> X (width)
    #   SYMMETRY
    #
    s.rectangle(point1=(0.0, 0.0), point2=(sample_W, sample_D))
    p = mdb.models['Model-1'].Part(name='SPA', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['SPA']
    p.BaseSolidExtrude(sketch=s, depth=sample_L)
    s.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']

    #--------------------------------------------------------------------------------
    # Define some variables.
    #--------------------------------------------------------------------------------
    d = p.datums
    c = p.cells
    f = p.faces

    #--------------------------------------------------------------------------------
    # Create sets and surfaces.
    # The global coordinate system has the origin at the bottom left of the symmetry
    # face, with body length defined in Z.
    #--------------------------------------------------------------------------------
    p.Set(cells=c, name='ALL')
    # Symmetry boundary condition.
    faces = f.findAt(((0.1, 0.1, 0.0),),)
    p.Set(faces=faces, name='SYMM_Z')
    # Strain boundary condition.
    faces = f.findAt(((0.1, 0.1, sample_L),),)
    p.Set(faces=faces, name='STRAIN_BC')

    #--------------------------------------------------------------------------------
    # Mesh the geometry.
    #--------------------------------------------------------------------------------
    elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD,
        kinematicSplit=AVERAGE_STRAIN, hourglassControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)
    pickedRegions =(c, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))
    # p.seedPart(size=sample_D/4.0, deviationFactor=0.1, minSizeFactor=0.1)
    # p.seedPart(size=sample_D, deviationFactor=0.1, minSizeFactor=0.1)
    p.seedPart(size=1000.0, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    #--------------------------------------------------------------------------------
    # Create and assign the material.
    #--------------------------------------------------------------------------------
    if mat['density']==0.0: print 'WARNING: material density equals 0.0'
    if mat['matname']=='???': print 'WARNING: material name is given as',mat['matname']
    mdb.models['Model-1'].Material(name=mat['matname'])
    mdb.models['Model-1'].materials[mat['matname']].Density(table=((mat['density'], ), ))
    allparams = mat['params']+mat['D']
    if 'order' in mat:
        mdb.models['Model-1'].materials[mat['matname']].Hyperelastic(
            materialType=ISOTROPIC, testData=OFF, volumetricResponse=VOLUMETRIC_DATA,
            # moduliTimeScale=INSTANTANEOUS, type=get_abq_model(mat['model']),
            moduliTimeScale=LONG_TERM, type=get_abq_model(mat['model']),
            n=mat['order'], table=(allparams, ))
    else:
        mdb.models['Model-1'].materials[mat['matname']].Hyperelastic(
            materialType=ISOTROPIC, testData=OFF, volumetricResponse=VOLUMETRIC_DATA,
            # moduliTimeScale=INSTANTANEOUS, type=get_abq_model(mat['model']),
            moduliTimeScale=LONG_TERM, type=get_abq_model(mat['model']),
            table=(allparams, ))
    mdb.models['Model-1'].materials[mat['matname']].hyperelastic.Hysteresis(
        table=((args.S, args.A, args.m, args.C), ))
    mdb.models['Model-1'].HomogeneousSolidSection(name='SPA',material=mat['matname'],thickness=None)
    region = p.sets['ALL']
    p.SectionAssignment(region=region, sectionName='SPA', offset=0.0,
        offsetType=MIDDLE_SURFACE, offsetField='',thicknessAssignment=FROM_SECTION)


    #--------------------------------------------------------------------------------
    # Create the assembly.
    #--------------------------------------------------------------------------------
    a = mdb.models['Model-1'].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    a.Instance(name='SPA', part=p, dependent=ON)

    #--------------------------------------------------------------------------------
    # Create a step.
    #--------------------------------------------------------------------------------
    # Defaults: initialInc=1.0, maxNumInc=100, minInc=1e-5, maxInc=1.0
    mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial', timePeriod=time,
            nlgeom=ON, initialInc=0.02*time, maxNumInc=1000, minInc=1e-06, maxInc=0.02*time)
            # , extrapolation=PARABOLIC)
    # mdb.models['Model-1'].steps['Step-1'].setValues(useLongTermSolution=True) # Only applies to static steps.
    # mdb.models['Model-1'].ViscoStep(name='Step-1', previous='Initial', cetol=0.1000, amplitude=RAMP,
            # timePeriod=time, nlgeom=ON,initialInc=0.1, maxNumInc=500, minInc=1e-06, maxInc=0.1)
    # Enable line search, see manual section 7.2.2
    mdb.models['Model-1'].steps['Step-1'].control.setValues(allowPropagation=OFF,
            resetDefaultValues=OFF, lineSearch=(5.0, 1.0, 0.0001, 0.25, 0.1))

    #--------------------------------------------------------------------------------
    # Create a reference point for calculating the length of the actuator.
    #--------------------------------------------------------------------------------
    # Create a reference point, and then make a set out of it.
    a.ReferencePoint(point=(0.5*sample_W, 0.5*sample_D, sample_L))
    rp = a.referencePoints[a.referencePoints.keys()[0]]
    a.Set(referencePoints=(rp,), name='EPT')
    # Create a coupling constraint, so that the point follows the average movement
    # of the measured surface.
    region1 = a.sets['EPT']
    region2 = a.instances['SPA'].sets['STRAIN_BC']
    mdb.models['Model-1'].Coupling(name='MEASUREMENT', controlPoint=region1,
        surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=DISTRIBUTING,
        weightingMethod=UNIFORM, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON,
        ur2=ON, ur3=ON)
    # Create a history output request for the displacement at the RP.
    regionDef=a.sets['EPT']
    mdb.models['Model-1'].HistoryOutputRequest(name='MEASUREMENT',
        createStepName='Step-1', region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE,
        variables=('U3', 'RF3', 'TF3'))

    #--------------------------------------------------------------------------------
    # Create the loads and boundary conditions.
    #--------------------------------------------------------------------------------
    region = a.instances['SPA'].sets['SYMM_Z']
    mdb.models['Model-1'].ZsymmBC(name='SYMM_Z',createStepName='Step-1',region=region,localCsys=None)
    # region = a.instances['SPA'].sets['STRAIN_BC']
    region = a.sets['EPT']
    mdb.models['Model-1'].DisplacementBC(name='STRAIN_BC', createStepName='Step-1',
        region=region, u1=UNSET, u2=UNSET, u3=displacement, ur1=UNSET, ur2=UNSET, ur3=UNSET,
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)


    #--------------------------------------------------------------------------------
    # Adjust the output requests.
    #--------------------------------------------------------------------------------
    mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
        'S', 'PE', 'PEEQ', 'PEMAG', 'EE', 'IE', 'NE', 'LE', 'U', 'RF', 'TF',
        'CF', 'P', 'CSTRESS', 'CDISP', 'CFORCE'))



    #================================================================================
    # Create the job and save the file.
    # Note: number of CPUs isn't written to the .inp file, it's only used if you
    #       run this job from the GUI.
    #================================================================================
    mdb.Job(name='test_visco', model='Model-1', description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE, nodalOutputPrecision=FULL, echoPrint=OFF,
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
        scratch='', multiprocessingMode=DEFAULT, numCpus=8, numDomains=8,
        numGPUs=0)
    mdb.jobs['test_visco'].writeInput(consistencyChecking=OFF)
    mdb.saveAs(pathName='test_visco.cae')
    # Open .inp file to add diagnostics message.
    with open('test_visco.tmp','w') as of, open('test_visco.inp','r') as f:
        for line in f:
            line = line.rstrip()    # Remove newline/carriage-return character.
            if line=='*End Step':
                of.write('**\n')
                of.write('** Disable error message about Poisson\'s ratio (see bug #3).\n')
                of.write('**\n')
                of.write('*DIAGNOSTICS, NONHYBRID INCOMPRESSIBLE=WARNING\n')
            of.write(line+'\n')
    shutil.move('test_visco.tmp','test_visco.inp')

