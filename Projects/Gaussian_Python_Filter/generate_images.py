__author__ = 'adherron'
__version__ = '1.0.1'
__doc__ = '''
generate_images is designed to run after vtk_reader has processed a directory full of VTK files. It is a script that \n
prompts the user for a series of directories, opens all the hdf5 files in that directory, and generates images from \n
their mapped screen_i data at a desired exposure level It saves these images as .png files of the same base name as \n
the .h5 file they are created from. CAUTION! It will overwrite any other images of the same name without warning!
'''
from numpy import *
from DiffractionSim import *
import os

paths = {}
new_path = raw_input("Enter a directory path containing saved DiffractionSim.py (hdf5) files:\n")
if len(new_path) <= 0:
    print "ERROR: No path specified!"
    raise IOError
paths[new_path] = []
while True:
    new_path = raw_input('\nEnter another directory path or hit ENTER to continue:\n')
    if len(new_path) <= 0:
        break
    else:
        paths[new_path] = []
exposure_c = input("Please enter the desired exposure level [float > 0]: ")
# shell_thickness = input('Please enter the desired Ewald Shell Thickness (inverse Angstroms) [float > 0]: ')
shell_thickness = 0.04
for path in paths:
    print "\nMoved to: ", path
    for h5_file in os.listdir(path):
        if h5_file[-2::] == 'h5':
            paths[path].append(h5_file)
for path in paths:
    print "Now Entering: ", path
    os.chdir(path)
    for f_name in paths[path]:
        print "Processing", f_name, ' . . . '
        s_name = f_name[0:-2] + 'png'
        A = DiffractionSim()
        A.read_data(f_name)
        if (A.zone[0] == 3) and (A.zone[1] == 0) and (A.zone[2] == 2):
            A.zone[2] = 2.14
        elif (A.zone[0] == 0) and (A.zone[1] == 3) and (A.zone[2] == 2):
            A.zone[2] = 2.14
        A.map2saed_screen(pixel_dim=512, shell_thickness=shell_thickness)
        A.image_screen(A.screen_i, exposure_c, 2, 2, plot=False, file_name=s_name)
