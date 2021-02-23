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

# Wrapper script for the abq_create_geom.py command.
import argparse, sys, os, time, random, shutil
import os.path as path
from datetime import datetime
import utility
ABAQUS='abaqus'                # Assume Abaqus is on the path.


#--------------------------------------------------------------------------------
# Main.
#--------------------------------------------------------------------------------
if __name__ == '__main__':
    tinit = time.time()
    #------------------------------------------------------------------------------------
    # Handle user input.
    #------------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='Create a geometry and mesh, ready to run in Abaqus.',
            epilog='Example: create_geom.py actuator lin U ogden3.mat outfile.cae 0.05 200 4 --mesh_size 2.0 --chamber 8 8 2 --wall 7')
    subparsers = parser.add_subparsers(dest='cmd',help='Choose a geometry type.')

    parser_ACT = subparsers.add_parser('actuator',help='Create a linear or a bending actuator. '
            'Note that all widths will be divided by two for symmetry, and linear actuator heights will also be divided by two for symmetry.')
    parser_ACT.add_argument('actuator',choices=['lin','bnd'],help='Type of actuator to create.')
    parser_ACT.add_argument('test',choices=['U','F'],help='Type of test to run.')
    parser_ACT.add_argument('matfile',help='File containing material properties.')
    parser_ACT.add_argument('cae',help='Output Abaqus database file to create.')
    parser_ACT.add_argument('press',help='Pressure in chamber, in N/mm2 (0.1N/mm2=100kPa).')
    parser_ACT.add_argument('time',help='Pseudo-time to run simulation (unitless unless using viscoelasticity, try 200).')
    parser_ACT.add_argument('num_chambers',help='Number of chambers.')
    parser_ACT.add_argument('--visco',default='none',help='Optional viscoelastic material definition.')
    parser_ACT.add_argument('--randomize',nargs=2,metavar=('PARAMS','RANGE'),help='Add a random value to certain parameters. '
            'PARAMS are comma separated, eg "ch,cw,cl,ih,iw,w" or "all", RANGE defines width of random spread.')
    parser_ACT.add_argument('--mesh_size','-m',default='2.0',metavar='S',help='Approximate size of mesh elements (default 2.0).')
    parser_ACT.add_argument('--chamber',type=float,nargs=3,default=[6.0,6.0,2.0],metavar=('H','W','L'),
            help='Height, width, and length of a chamber (default \'6.0 6.0 2.0\').')
    parser_ACT.add_argument('--inlet',type=float,nargs=2,default=[2.0,2.0],metavar=('H','W'),
            help='Height and width of inlet tunnel (default \'2.0 2.0\').')
    parser_ACT.add_argument('--wall',type=float,default=3.0,help='Wall thickness (default 3.0).')
    parser_ACT.add_argument('--mass_scaling',default='1.0',help='Mass scaling factor (default 1.0, ie disabled).')
    parser_ACT.add_argument('--ale',action='store_true',default=False,help='Use an Arbitrary-Lagrangian-Eularian mesh.')
    parser_ACT.add_argument('--doubleIC',action='store_true',default=False,help='Doubles the interchamber spacing (2x wall).')
    parser_ACT.add_argument('--dist_to_force',default='0.0',help='Distance to the blocked force plate in blocked-force tests (default 0.0), in mm for linF or degrees for bndF.')
    parser_ACT.add_argument('--maxnuminc',default='100000',help='Maximum number of increments (default 100000).')

    parser_BUB = subparsers.add_parser('bubble',help='Create a bubble-style actuator. All dimensions will be divided by two for symmetry.')
    parser_BUB.add_argument('test',choices=['U','F'],help='Type of test to run.')
    parser_BUB.add_argument('matfile',help='File containing material properties.')
    parser_BUB.add_argument('cae',help='Output Abaqus database file to create.')
    parser_BUB.add_argument('press',help='Pressure in chamber, in N/mm2 (0.1=100kPa).')
    parser_BUB.add_argument('time',help='Pseudo-time to run simulation (unitless unless using viscocity, try 200).')
    parser_BUB.add_argument('--diameter',type=float,default=8.0,help='Bubble diameter (default 8.0).')
    parser_BUB.add_argument('--thickness',type=float,default=3.0,help='Film thickness (default 3.0).')
    parser_BUB.add_argument('--inlet',type=float,default=2.0,help='Inlet width (default 2.0).')
    parser_BUB.add_argument('--mesh_size','-m',default='2.0',metavar='S',help='Approximate size of mesh elements (default 2.0).')

    args = parser.parse_args()

    #------------------------------------------------------------------------------------
    # Prep the arguments.
    #------------------------------------------------------------------------------------
    # Create the directory, if necessary.
    caedir = path.split(args.cae)[0]
    if caedir=='': caedir = os.getcwd()
    elif not os.path.exists(caedir): os.makedirs(caedir)
    # If the cae file already exists, remove it because Abaqus won't necessarily overwrite it.
    if path.isfile(args.cae): os.remove(args.cae)
    # Change directory, since Abaqus is limited in where it can write files.
    if   args.cmd=='actuator': acg_file = utility.find_file('python/abq_create_geom.py')
    elif args.cmd=='bubble':   acg_file = utility.find_file('python/abq_create_bubble.py')
    rootdir = os.getcwd()
    os.chdir(caedir)
    acg_file = path.normpath(path.join(rootdir,acg_file))
    cae_file = path.split(args.cae)[1]
    inp_file = path.splitext(cae_file)[0]+'.inp'
    mat_file = path.normpath(path.join(rootdir,args.matfile))
    if args.cmd=='actuator': vis_file = path.normpath(path.join(rootdir,args.visco)) if args.visco!='none' else 'none'

    #------------------------------------------------------------------------------------
    # Linear and bending actuators.
    #------------------------------------------------------------------------------------
    if args.cmd=='actuator':
        # Determine the dimensions.
        if args.randomize:
            random.seed()
            if args.randomize[0]=='all': args.randomize[0]='ih,iw,ch,cw,cl,w'
            rp = args.randomize[0].split(',')
            scale = float(args.randomize[1])
        else: rp = []
        NC = int(args.num_chambers)
        IH = args.inlet[0]   + (random.uniform(-0.5,0.5)*scale if 'ih' in rp else 0.0)
        IW = args.inlet[1]   + (random.uniform(-0.5,0.5)*scale if 'iw' in rp else 0.0)
        CH = args.chamber[0] + (random.uniform(-0.5,0.5)*scale if 'ch' in rp else 0.0)
        CW = args.chamber[1] + (random.uniform(-0.5,0.5)*scale if 'cw' in rp else 0.0)
        CL = args.chamber[2] + (random.uniform(-0.5,0.5)*scale if 'cl' in rp else 0.0)
        W  = args.wall       + (random.uniform(-0.5,0.5)*scale if 'w'  in rp else 0.0)
        TL = NC*CL + 2.0*W + (NC-1.0)*W*(2.0 if args.doubleIC else 1.0)

        # Create the geometry.
        print 'Creating Abaqus geometry with dimensions:'
        print '\tInlet height:   ',IH,('\t\t(randomized from '+str(args.inlet[0])+')' if 'ih' in rp else '')
        print '\tInlet width:    ',IW,('\t\t(randomized from '+str(args.inlet[1])+')' if 'iw' in rp else '')
        print '\tChamber height: ',CH,('\t\t(randomized from '+str(args.chamber[0])+')' if 'ch' in rp else '')
        print '\tChamber width:  ',CW,('\t\t(randomized from '+str(args.chamber[1])+')' if 'cw' in rp else '')
        print '\tChamber length: ',CL,('\t\t(randomized from '+str(args.chamber[2])+')' if 'cl' in rp else '')
        print '\tWall thickness: ',W,('\t\t(randomized from '+str(args.wall)+')' if 'w' in rp else '')
        print '\tChamber spacing:',('2x_wall' if args.doubleIC else '1x_wall')
        print '\tTOTAL LENGTH:   ',TL
        ale = 'ale' if args.ale else 'noale'
        cspace = '2wall' if args.doubleIC else '1wall'
        acg_args = [cae_file, ale, args.actuator, args.test, args.mesh_size, args.press, args.time, args.num_chambers,
                    cspace, str(IH), str(IW), str(CH), str(CW), str(CL), str(W), args.mass_scaling, mat_file, vis_file,
                    args.dist_to_force, args.maxnuminc]
        cmd = [ABAQUS,'cae','noGUI='+acg_file,'--']+acg_args
    #------------------------------------------------------------------------------------
    # Bubble actuators.
    #------------------------------------------------------------------------------------
    if args.cmd=='bubble':
        acg_args = [cae_file, args.test, args.mesh_size, args.press, args.time, str(args.diameter),
                    str(args.thickness), str(args.inlet), mat_file]
        cmd = [ABAQUS,'cae','noGUI='+acg_file,'--']+acg_args

    #------------------------------------------------------------------------------------
    # Attempt to run Abaqus to create the geometry.
    #------------------------------------------------------------------------------------
    try:
        utility.run_cmd_screen(cmd)
    except Exception:
        utility.print_error(sys.exc_info()[0],False)
        print 'Check above to see if Abaqus printed an error. Otherwise, check the Abaqus log files.'
        print 'Tried to run command:'
        print ' '.join(cmd)
        print 'In directory:',caedir
        print
    else:
        # Cleanup.
        os.remove('abaqus.rpy')
        os.remove(path.splitext(cae_file)[0]+'.jnl')
        # Read .inp file to determine number of elements and nodes.
        with open(inp_file+'.tmp','w') as of, open(inp_file,'r') as f:
            countN = False
            countE = False
            nn = 0
            ne = 0
            for line in f:
                line = line.rstrip()    # Remove newline/carriage-return character.
                if line.startswith('** Generated by: Abaqus/CAE'):
                    of.write('** Created: '+str(datetime.now())+'\n')
                    if args.cmd=='actuator':
                        of.write('**   Actuator type:  '+args.actuator+args.test+'\n')
                        of.write('**   Mesh size:      '+args.mesh_size+'\n')
                        of.write('**   Pressure:       '+args.press+'\n')
                        of.write('**   Num Chambers:   '+args.num_chambers+'\n')
                        of.write('**   Inlet height:   '+str(IH)+('\t\t(randomized from '+str(args.inlet[0])+')\n' if 'ih' in rp else '\n'))
                        of.write('**   Inlet width:    '+str(IW)+('\t\t(randomized from '+str(args.inlet[1])+')\n' if 'iw' in rp else '\n'))
                        of.write('**   Chamber height: '+str(CH)+('\t\t(randomized from '+str(args.chamber[0])+')\n' if 'ch' in rp else '\n'))
                        of.write('**   Chamber width:  '+str(CW)+('\t\t(randomized from '+str(args.chamber[1])+')\n' if 'cw' in rp else '\n'))
                        of.write('**   Chamber length: '+str(CL)+('\t\t(randomized from '+str(args.chamber[2])+')\n' if 'cl' in rp else '\n'))
                        of.write('**   Wall thickness: '+str(W)+ ('\t\t(randomized from '+str(args.wall)+')\n' if 'w' in rp else '\n'))
                        of.write('**   Chmaber spacing:'+('2x_wall' if args.doubleIC else '1x_wall')+'\n')
                        of.write('**   TOTAL LENGTH:   '+str(TL)+'\n')
                # if line=='*End Step':
                    # of.write('**\n')
                    # of.write('** Disable error message about Poisson\'s ratio (see bug #3).\n')
                    # of.write('**\n')
                    # of.write('*DIAGNOSTICS, NONHYBRID INCOMPRESSIBLE=WARNING\n')
                of.write(line+'\n')
                if line.startswith('*'):
                    countN=False
                    countE=False
                if   line=='*Node': countN=True
                elif line.startswith('*Element, type='): countE=True
                elif countN: nn = nn+1
                elif countE: ne = ne+1
        shutil.move(inp_file+'.tmp',inp_file)
        print 'Created','{:,}'.format(ne),'elements and','{:,}'.format(nn),'nodes.'
    finally:
        os.chdir(rootdir)
        print 'TIME ELAPSED: ',utility.time_elapsed(tinit)
