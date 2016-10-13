__author__ = 'adherron'
__doc__ = '''\n
vtk_reader is a script designed to convert a large number of VTK files output by Shawn Coleman's fix_time_ave_saed \n
int hdf5 files that can be read by the DiffractionSim.py class. It prompts the user to enter a series of directories and \n
searches those directories for VTK files. It prompts the user to enter the zone and Atom Count information for each \n
of those files first and then processes them all afterward, creating a new file of the same name as the VTK file, but \n
with the ending .h5 instead of .vtk\n
CAUTION: this will overwrite any other existing h5 file of the same name without warning!\n
'''

from numpy import *
from DiffractionSim import *
import os

paths = {}
print "vtk_reader asks for directories and converts all VTK files in those directories to hdf5 files\n"
new_path = raw_input("Please enter a directory path containing VTK files:\n")
if len(new_path) <= 0:
    print "ERROR: No path specified!"
    raise IOError
paths[new_path] = {}
while True:
    new_path = raw_input('\nEnter another directory path or hit ENTER to continue:\n')
    if len(new_path) <= 0:
        break
    else:
        paths[new_path] = {}
print "After Reading Each File Name, Type Its Zone x y z and press" \
      " ENTER then type the Atom Count for that structure and press ENTER:"
for path in paths:
    print "\nMoved to: ", path
    for vtk_file in os.listdir(path):
        if vtk_file[-3::] == 'vtk':
            while True:
                try:
                    zone = []
                    print "\n", vtk_file
                    z = raw_input("Zone x y z: ")
                    for string in z.split(' '):
                        zone.append(float(string))
                    if zone == [1., 1., 1.]:
                        zone = [3., 0., 2.14]
                    if zone == [1., -1., 1.]:
                        zone = [0., 3., 2.14]
                    if zone == [1., 1., 0.]:
                        zone = [1., 0., 0.]
                    if zone == [1., -1., 0.]:
                        zone = [0., 1., 0.]
                    paths[path][vtk_file] = {}
                    paths[path][vtk_file]['zone'] = array(zone)
                    break
                except ValueError:
                    print "Try Again - not correct x y z format for zone assignment!"
            paths[path][vtk_file]['n_atoms'] = input("Atom Count: ")
for path in paths:
    print "Now Entering: ", path
    os.chdir(path)
    for vf in paths[path]:
        f_name = vf
        h_name = f_name[0:-3] + 'h5'
        zone = paths[path][vf]['zone']
        print "Processing", f_name, "with zone: ", zone, " . . ."

        A = DiffractionSim()
        A.read_data(f_name, zone)
        A.n_atoms = paths[path][vf]['n_atoms']
        # A.map2saed_screen(pixel_dim=512)
        A.save(h_name)
