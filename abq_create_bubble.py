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

# Use this script to create Abaqus geometry for a SPA.
# Note that you probably never want to run this script directly;
# instead prefer to use the create_geom.py wrapper script.
# Run using: abaqus cae noGUI=abq_create_geom.py -- --help
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
    NARGS = 9
    run_cmd = 'abaqus cae noGUI=abq_create_bubble.py -- '
    if '--help' in sys.argv: sys.argv = [run_cmd,'--help']
    else: sys.argv = [run_cmd]+sys.argv[len(sys.argv)-NARGS:]
    parser = argparse.ArgumentParser(
            description='Create Abaqus geometry for a bubble-style SPA.',
            epilog='Example: '+run_cmd+'test.cae U 2.0 0.05 8.0 3.0 2.0 matfile.mat')
    parser.add_argument('cae',help='Output Abaqus database file to create.')
    parser.add_argument('test',choices=['U','F'],help='Type of test to run.')
    parser.add_argument('mesh_size',type=float,help='Approximate size of mesh elements.')
    parser.add_argument('pressure',type=float,help='Pressure in chamber, probably in N/mm2. 0.1=100kPa. If not converging, try really small (0.01) values here.')
    parser.add_argument('time',type=float,help='Time to end simulation (s).')
    parser.add_argument('diameter',type=float,help='Diameter of bubble (will be divided by two for symmetry).')
    parser.add_argument('thickness',type=float,help='Wall thickness (will be divided by two for symmetry).')
    parser.add_argument('inlet',type=float,help='Width of inlet tunnel (will be divided by two for symmetry).')
    parser.add_argument('matfile',help='File containing material properties.')
    args = parser.parse_args()
    if not args.cae.endswith('.cae'):
        print 'WARNING: incorrect inputs. Found arguments:'
        print sys.argv
        exit(1)
    #--------------------------------------------------------------------------------
    # Read the material file.
    #--------------------------------------------------------------------------------
    mat = MU.read_matfile(args.matfile)

    #--------------------------------------------------------------------------------
    # Calculate dimensions, assuming mirror-symmetry in X and Z.
    #   width = X
    #   height = Y
    #   thickness = Z
    #--------------------------------------------------------------------------------
    mesh_size = args.mesh_size
    diameter = args.diameter
    radius = args.diameter / 2.0
    thickness = args.thickness / 2.0
    inlet = args.inlet / 2.0

    buff = radius
    inlet_buff = 5.0*inlet
    body_W = radius + buff
    body_H = diameter + 2.0*buff + inlet_buff
    time = 200


    #--------------------------------------------------------------------------------
    # Open the viewport.
    #--------------------------------------------------------------------------------
    session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=300, height=200)
    session.viewports['Viewport: 1'].maximize()
    executeOnCaeStartup()
    Mdb()

    #--------------------------------------------------------------------------------
    # Create the part.
    #--------------------------------------------------------------------------------
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',sheetSize=20.0)
    s.setPrimaryObject(option=STANDALONE)

    # Extrude the shape for the glued sheet.
    #   O = origin
    #   H = body_H
    #   W = body_W
    #
    #        +------W-----+
    #        |            |
    #        +--__        |
    #             \       |
    #              |      H
    #             /       |
    #          +-         |
    #          |          |
    #          |          |
    #        O +----------+
    #
    s.ArcByCenterEnds(center=(0.0, body_H-buff-radius), point1=(0.0, body_H-buff), point2=(0.0, body_H-diameter-buff), direction=CLOCKWISE)
    s.Line(point1=(inlet, body_H-buff-radius), point2=(inlet, 0.0))
    s.Line(point1=(inlet, 0.0), point2=(body_W, 0.0))
    s.Line(point1=(body_W, 0.0), point2=(body_W, body_H))
    s.Line(point1=(body_W, body_H), point2=(0.0, body_H))
    s.Line(point1=(0.0, body_H), point2=(0.0, body_H-buff))
    s.autoTrimCurve(curve1=s.geometry.findAt((inlet, 0.1)), point1=(inlet, body_H-buff-radius-0.1))
    s.autoTrimCurve(curve1=s.geometry.findAt((0.0, body_H-buff)), point1=(0.0, body_H-diameter-buff))
    p = mdb.models['Model-1'].Part(name='BUBBLE', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['BUBBLE']
    p.BaseSolidExtrude(sketch=s, depth=thickness)
    s.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts['BUBBLE']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']


    # Extrude the shape for the free sheet.
    #   O = origin
    #
    #        +--__
    #        |    \
    #        |     |
    #        |    /
    #        | +-
    #        | |
    #        | |
    #        O-+
    #
    f, e = p.faces, p.edges
    t = p.MakeSketchTransform(sketchPlane=f.findAt(coordinates=(inlet+0.1, 0.1, thickness)),
            sketchUpEdge=e.findAt(coordinates=(body_W, 0.1, thickness)),
            sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, thickness))
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=23.82, gridSpacing=0.59, transform=t)
    s.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    s.ArcByCenterEnds(center=(0.0, body_H-buff-radius), point1=(0.0, body_H-buff), point2=(0.0, body_H-diameter-buff), direction=CLOCKWISE)
    s.Line(point1=(inlet, body_H-buff-radius), point2=(inlet, 0.0))
    s.Line(point1=(inlet, 0.0), point2=(0.0, 0.0))
    s.Line(point1=(0.0, 0.0), point2=(0.0, body_H-buff))
    s.autoTrimCurve(curve1=s.geometry.findAt((inlet, 0.1)), point1=(inlet, body_H-buff-radius-0.1))
    s.autoTrimCurve(curve1=s.geometry.findAt((0.0, body_H-buff)), point1=(0.0, body_H-diameter-buff))
    p.SolidExtrude(sketchPlane=f.findAt(coordinates=(inlet+0.1, 0.1, thickness)),
        sketchUpEdge=e.findAt(coordinates=(body_W, 0.1, thickness)), sketchPlaneSide=SIDE1,
        sketchOrientation=RIGHT, sketch=s, depth=thickness, flipExtrudeDirection=ON,
        keepInternalBoundaries=ON)
    s.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']

    #--------------------------------------------------------------------------------
    # Define some variables.
    #--------------------------------------------------------------------------------
    d = p.datums
    c = p.cells
    f = p.faces
    e = p.edges
    v = p.vertices

    #--------------------------------------------------------------------------------
    # Create sets and surfaces. See diagram above for global origin and directions.
    #--------------------------------------------------------------------------------
    p.Set(cells=c, name='ALL')
    # Symmetry boundary condition in thickness.
    faces = f.findAt(((inlet+0.1, 0.1, 0.0),),)
    p.Set(faces=faces, name='SYMM_T')
    # Symmetry boundary condition in width.
    faces = f.findAt(((0.0, 0.1, 0.1),),
                     ((0.0, body_H-0.1, 0.1),),)
    p.Set(faces=faces, name='SYMM_W')
    # Pressure set, for measuring purposes.
    faces = f.findAt(((0.1, 0.1, 0.0),),)
    p.Set(faces=faces, name='PRESSURE')
    # Set for measuring forces in blocked-force tests.
    if args.test=='F':
        faces = f.findAt(((0.1,       0.1, thickness),),
                         ((inlet+0.1, 0.1, thickness),),)
        p.Surface(side1Faces=faces, name='BLOCKED')

    # Fixed boundary condition.
    faces = f.findAt(((body_W-0.1, 0.0,    0.1),),  # Bottom edge.
                     ((body_W,     0.1,    0.1),),  # Right edge.
                     ((0.1,        body_H, 0.1),),) # Top edge.
    p.Set(faces=faces, name='FIXED')

    # Pressure boundary condition, in chambers and connecting tubes.
    faces = f.findAt(((0.1, 0.1, 0.0),),)
    p.Surface(side1Faces=faces, name='PRESSURE')

    # Create a partition so we have a measuring point for displacements.
    if args.test=='U':
        pickedCells = c.findAt(((0.1,0.1,0.1),))
        p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=body_H-buff-radius)
        p.PartitionCellByDatumPlane(datumPlane=d.values()[-1], cells=pickedCells)


    #--------------------------------------------------------------------------------
    # Mesh the geometry.
    #--------------------------------------------------------------------------------
    # Quadratic elements. These provide a much faster rate of convergence compared to linear elements.
    elemType1 = mesh.ElemType(elemCode=C3D20RH, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=C3D15, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD)
    # Linear elements.
    # elemType1 = mesh.ElemType(elemCode=C3D8RH, elemLibrary=STANDARD, kinematicSplit=AVERAGE_STRAIN, hourglassControl=DEFAULT)
    # elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
    # elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)
    region = p.sets['ALL']
    p.setElementType(regions=region, elemTypes=(elemType1, elemType2, elemType3))
    p.seedPart(size=args.mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    #--------------------------------------------------------------------------------
    # Set the output frequency.
    #--------------------------------------------------------------------------------
    num_nodes = len(p.nodes)
    output_frequency = 1

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
            type=get_abq_model(mat['model']), n=mat['order'], table=(allparams, ),
            moduliTimeScale=LONG_TERM)
            # moduliTimeScale=INSTANTANEOUS)
    else:
        mdb.models['Model-1'].materials[mat['matname']].Hyperelastic(
            materialType=ISOTROPIC, testData=OFF, volumetricResponse=VOLUMETRIC_DATA,
            type=get_abq_model(mat['model']), table=(allparams, ),
            moduliTimeScale=LONG_TERM)
            # moduliTimeScale=INSTANTANEOUS)
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
    inc = 0.1*time
    # Defaults: initialInc=1.0, maxNumInc=100, minInc=1e-5, maxInc=1.0
    # mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial', timePeriod=time,
            # nlgeom=ON, initialInc=inc, maxNumInc=1000, minInc=1e-06, maxInc=inc,
            # solutionTechnique=QUASI_NEWTON,
            # extrapolation=PARABOLIC)
    # Enable stabilization.
    # mdb.models['Model-1'].steps['Step-1'].setValues(stabilizationMagnitude=0.0002,
        # stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
        # continueDampingFactors=False, adaptiveDampingRatio=0.05)
    # mdb.models['Model-1'].ExplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=time)
    # For bending actuators with membrane elements we use an implicit dynamics step.
    # Unloaded membrane elements can rapidly become unstable in a static analysis
    # because they have no bending stiffness.
    # For linear actuators, they also let us converge at higher pressures.
    mdb.models['Model-1'].ImplicitDynamicsStep(name='Step-1', previous='Initial',
        timePeriod=time, application=QUASI_STATIC, initialInc=inc, minInc=1e-06, maxInc=inc,
        maxNumInc=10000, nohaf=OFF, amplitude=RAMP, alpha=DEFAULT, initialConditions=OFF, nlgeom=ON)
        # solutionTechnique=QUASI_NEWTON,
    # mdb.models['Model-1'].materials[mat['matname']].Damping(alpha=0.002*time, beta=0.0)
    # mdb.models['Model-1'].materials[mat['matname']].Damping(alpha=0.4, beta=0.0)
    # Enable line search for Quasi-Newton algorithm, see manual section 7.2.2
    # mdb.models['Model-1'].steps['Step-1'].control.setValues(allowPropagation=OFF,
            # resetDefaultValues=OFF, lineSearch=(5.0, 1.0, 0.0001, 0.25, 0.1))

    #--------------------------------------------------------------------------------
    # Create the loads and boundary conditions.
    #--------------------------------------------------------------------------------
    region = a.instances['SPA'].sets['SYMM_T']
    mdb.models['Model-1'].ZsymmBC(name='SYMM_T',createStepName='Step-1',region=region,localCsys=None)
    region = a.instances['SPA'].sets['SYMM_W']
    mdb.models['Model-1'].XsymmBC(name='SYMM_W',createStepName='Step-1',region=region,localCsys=None)
    region = a.instances['SPA'].sets['FIXED']
    mdb.models['Model-1'].DisplacementBC(name='FIXED', createStepName='Initial',
        region=region, u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET,
        amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
    region = a.instances['SPA'].surfaces['PRESSURE']
    mdb.models['Model-1'].Pressure(name='PRESSURE', createStepName='Step-1',
        region=region, distributionType=UNIFORM, field='', magnitude=args.pressure,
        amplitude=UNSET)

    # If this is blocked-force testing, set additional BCs.
    # if args.test=='F':
        # region = a.instances['SPA'].sets['BC_CUPS']
        # mdb.models['Model-1'].DisplacementBC(name='CUPS', createStepName='Initial',
            # region=region, u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET,
            # amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)


    #--------------------------------------------------------------------------------
    # Adjust the output requests.
    #--------------------------------------------------------------------------------
    mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(frequency=output_frequency,variables=(
        'S',                          # Stress components and invariants.
        # 'PE', 'PEEQ', 'PEMAG',        # Plastic strain, Equivalent plastic strain, Plastic strain magnitude.
        'EE', 'IE', 'NE', 'LE',       # Elastic strain, Inelastic strain, Nominal strain, Logarithmic strain.
        'U', 'V', 'A',                # Displacement, Velocity, Acceleration.
        'RF', 'CF', 'P',              # Reaction forces and moments, Concentrated forces and moments, Pressure loads.
        'CSTRESS', 'CDISP', 'CFORCE', # Contact stresses, Contact displacements, Contact forces.
        'ENER',                       # All energy magnitudes.
        # 'EVOL',                       # Element volume.
        ))
    regionDef = a.allInstances['SPA'].sets['PRESSURE']
    mdb.models['Model-1'].FieldOutputRequest(name='PRESSURE', createStepName='Step-1',
        variables=('P', ), region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE, frequency=output_frequency)
    # H-Output-1 is the default history output request with the default variables.
    mdb.models['Model-1'].historyOutputRequests['H-Output-1'].setValues(frequency=output_frequency)





    #================================================================================
    # For a displacement test, add a nodeset to measure displacements.
    #================================================================================
    if args.test=='U':
        verts = v.findAt(((0.0, body_H-buff-radius, thickness), ))
        p.Set(vertices=verts, name='UPT')
        # Create a history output request for the displacement at the RP.
        regionDef = a.allInstances['SPA'].sets['UPT']
        mdb.models['Model-1'].HistoryOutputRequest(name='DISPLACEMENT', frequency=output_frequency,
            createStepName='Step-1', region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE,
            variables=('U1', 'U2', 'U3', 'COOR3'))


    #================================================================================
    # If the test is a blocked force test, add blocked force BCs.
    #================================================================================
    elif args.test=='F':
        # Create the analytic rigid surface part.
        s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=20.0)
        s.setPrimaryObject(option=STANDALONE)
        s.Line(point1=(0.0, 0.0), point2=(body_W, 0.0))
        s.HorizontalConstraint(entity=s.geometry[2], addUndoState=False)
        pr = mdb.models['Model-1'].Part(name='WALL', dimensionality=THREE_D, type=ANALYTIC_RIGID_SURFACE)
        pr = mdb.models['Model-1'].parts['WALL']
        pr.AnalyticRigidSurfExtrude(sketch=s, depth=body_H)  # Depth is for display only (infinite depth).
        s.unsetPrimaryObject()
        pr = mdb.models['Model-1'].parts['WALL']
        del mdb.models['Model-1'].sketches['__profile__']
        # Add an instance of the rigid part, and move into place.
        a.Instance(name='WALL', part=pr, dependent=ON)
        a.rotate(instanceList=('WALL', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(1.0, 0.0, 0.0), angle=90.0)
        a.translate(instanceList=('WALL', ), vector=(0.0, 0.5*body_H, thickness))
        # Create the reference points for rigid body constraints.
        a.ReferencePoint(point=(0.0, 0.5*body_H, thickness))
        a.features.changeKey(fromName='RP-1', toName='FORCE_RP')
        rp = a.referencePoints[a.referencePoints.keys()[0]]
        a.Set(referencePoints=(rp,), name='FPT')
        # Lock the reference point in place.
        region = a.sets['FPT']
        mdb.models['Model-1'].EncastreBC(name='WALL_RP', createStepName='Step-1', region=region, localCsys=None)
        # Constrain the rigid body "wall" to the reference point, which is fixed.
        s = a.instances['WALL'].faces
        faces = s.findAt(((0.1, 0.1, thickness), ))
        region1=regionToolset.Region(side1Faces=faces)
        region2 = a.sets['FPT']
        mdb.models['Model-1'].RigidBody(name='FIXWALL', refPointRegion=region2, surfaceRegion=region1)
        # Define the interaction properties between the SPA and the wall.
        mdb.models['Model-1'].ContactProperty('SPA_WALL')
        mdb.models['Model-1'].interactionProperties['SPA_WALL'].NormalBehavior(
            pressureOverclosure=HARD, allowSeparation=ON, constraintEnforcementMethod=DEFAULT)
        # Define which surfaces and methods are used to enforce interaction.
        s = a.instances['WALL'].faces
        faces = s.findAt(((0.1, 0.1, thickness), ))
        region1=regionToolset.Region(side2Faces=faces)
        region2=a.instances['SPA'].surfaces['BLOCKED']
        mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='SPA_WALL',
            createStepName='Step-1', master=region1, slave=region2, sliding=FINITE,
            enforcement=NODE_TO_SURFACE, thickness=OFF, interactionProperty='SPA_WALL',
            surfaceSmoothing=NONE, adjustMethod=NONE, smooth=0.2,
            # surfaceSmoothing=NONE, adjustMethod=OVERCLOSED, smooth=0.2, tied=OFF,
            initialClearance=OMIT, datumAxis=None, clearanceRegion=None)
        # Create a history output request for the blocked force at the RP.
        regionDef=a.sets['FPT']
        mdb.models['Model-1'].HistoryOutputRequest(name='BLOCKED_FORCE', frequency=output_frequency,
            createStepName='Step-1', region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE,
            variables=('RF1', 'RF2', 'RF3', 'COOR3'))


    #================================================================================
    # Create the job and save the file.
    # Note: number of CPUs isn't written to the .inp file, it's only used if you
    #       run this job from the GUI.
    #================================================================================
    jname = os.path.split(os.path.splitext(args.cae)[0])
    mdb.Job(name=jname[1], model='Model-1', description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE, nodalOutputPrecision=FULL, echoPrint=OFF,
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
        scratch='', multiprocessingMode=DEFAULT, numCpus=8, numDomains=8,
        numGPUs=0)
    if jname[0]!='':
        print 'WARNING: Abaqus can\'t write a .inp file anywhere other than the current directory.'
        print 'Abaqus will try to write it in the current directory, and then we will move it.'
        if os.path.isfile(jname[1]+'.inp'):
            print 'ERROR: file exists!'
            exit(1)
    mdb.jobs[jname[1]].writeInput(consistencyChecking=OFF)
    if jname[0]!='':
        path1 = jname[1]+'.inp'
        path2 = os.path.join(jname[0],path1)
        shutil.move(path1, path2)
        print 'Successfully moved the file from',path1,'to',path2
    mdb.saveAs(pathName=args.cae)
