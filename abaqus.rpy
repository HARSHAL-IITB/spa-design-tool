# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 6.14-1 replay file
# Internal Version: 2014_06_05-00.11.02 134264
# Run by admin on Fri May 08 13:32:39 2015
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(1.76302, 1.76389), width=259.517, 
    height=174.978)
session.viewports['Viewport: 1'].makeCurrent()
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
execfile('abq_create_geom.py', __main__.__dict__)
#: usage: abaqus cae noGUI=abq_create_geom.py --  [-h]
#:                                                cae {ale,noale} {lin,bnd} {U,F}
#:                                                mesh_size pressure time
#:                                                num_chambers {1wall,2wall}
#:                                                inlet_H inlet_W chamber_H
#:                                                chamber_W chamber_L wall
#:                                                mass_scaling matfile visfile
#:                                                dist_to_force maxnuminc
#: 
#: Create Abaqus geometry for a SPA.
#: 
#: positional arguments:
#:   cae            Output Abaqus database file to create.
#:   {ale,noale}    Whether or not to create an ALE mesh.
#:   {lin,bnd}      Type of actuator to create.
#:   {U,F}          Type of test to run.
#:   mesh_size      Approximate size of mesh elements.
#:   pressure       Pressure in chamber, probably in N/mm2. 0.1=100kPa. If not
#:                  converging, try really small (0.01) values here.
#:   time           Time to end simulation (s).
#:   num_chambers   Number of chambers.
#:   {1wall,2wall}  Spacing between chambers.
#:   inlet_H        Height of inlet tunnel (divided by two for linear actuators).
#:   inlet_W        Width of inlet tunnel (will be divided by two for symmetry).
#:   chamber_H      Height of chamber (divided by two for linear actuators).
#:   chamber_W      Width of chamber (will be divided by two for symmetry).
#:   chamber_L      Length of chamber.
#:   wall           Wall thickness.
#:   mass_scaling   Mass scaling factor (1.0 to disable).
#:   matfile        File containing material properties.
#:   visfile        File containing viscoelastic material properties (or 'none'.
#:   dist_to_force  Distance to blocked force plate in blocked-force tests.
#:   maxnuminc      Maximum number of increments.
#: 
#: optional arguments:
#:   -h, --help     show this help message and exit
#: 
#: Example: abaqus cae noGUI=abq_create_geom.py -- problem1.cae ale lin U 2.0
#: 0.05 200 4 1wall 2.0 2.0 6.0 6.0 2.0 3.0 1.0 ecoflex30.mat visfile.mat 0.0
#* Exit code: 0
#* File "abq_create_geom.py", line 71, in <module>
#*     args = parser.parse_args()
#* File "C:\SIMULIA\Abaqus\6.14-1\tools\SMApy\python2.7\lib\argparse.py", line 
#* 1688, in parse_args
#*     args, argv = self.parse_known_args(args, namespace)
#* File "C:\SIMULIA\Abaqus\6.14-1\tools\SMApy\python2.7\lib\argparse.py", line 
#* 1720, in parse_known_args
#*     namespace, args = self._parse_known_args(args, namespace)
#* File "C:\SIMULIA\Abaqus\6.14-1\tools\SMApy\python2.7\lib\argparse.py", line 
#* 1926, in _parse_known_args
#*     start_index = consume_optional(start_index)
#* File "C:\SIMULIA\Abaqus\6.14-1\tools\SMApy\python2.7\lib\argparse.py", line 
#* 1866, in consume_optional
#*     take_action(action, args, option_string)
#* File "C:\SIMULIA\Abaqus\6.14-1\tools\SMApy\python2.7\lib\argparse.py", line 
#* 1794, in take_action
#*     action(self, namespace, argument_values, option_string)
#* File "C:\SIMULIA\Abaqus\6.14-1\tools\SMApy\python2.7\lib\argparse.py", line 
#* 995, in __call__
#*     parser.exit()
#* File "C:\SIMULIA\Abaqus\6.14-1\tools\SMApy\python2.7\lib\argparse.py", line 
#* 2335, in exit
#*     _sys.exit(status)
