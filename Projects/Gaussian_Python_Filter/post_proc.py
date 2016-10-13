## This script will read a vtk, txt, or h5 file and create an image like those we created before - use h5 if possible
## It will also create an h5 version of any *.vtk file so that the next time you want to generate an image, it is much
## faster. Simply enter the path of the *.h5 file below for f_path on subsequent runs of this script

## edit the script first then run it by entering this:
##                   python post_proc.py

## If HDF5 gives you an version error when you try to run DiffractionSim scripts like this, run this first:
##                   export HDF5_DISABLE_VERSION_CHECK=2

from DiffractionSim import *


f_path = 'C:\Users\Adherron\Documents\Co_Ni_Al PROJECT\saed_data\New_SET3\Al29Co36Ni35_110MC900_16x16x1.h5'    # string with path to the desired vtk or hdf5 file (don't use '~/')
i_path = 'C:\Users\Adherron\Documents\Co_Ni_Al PROJECT\saed_data\New_SET3\Al29Co36Ni35_110MC900_16x16x1.png'    # string with path to the final image (eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, tiff formats accepted)
zone = [1, 0, 0]    # zone (in cartesian coordinates) ie: NiCoAl [111] is [3.0, 0.0, 2.14]
shell_thickness = 0.04    # total thickness of the post-process Ewald Shell Slice ( float <= 2 * Shawn Coleman's dR_Ewald in LAMMPS for best results)
n_atoms = 3072000    # integer number of atoms in structure (if unknown, enter -1 (scales from max intensity instead))
exposure_c = 6    # exposure level (float > 0)
pixel_dim = 512    # How many pixels across is the desired image?


A = DiffractionSim()

A.read_data(f_path, zone)

A.n_atoms = n_atoms

if (f_path[-3:] == 'vtk') or (f_path[-3:] == 'txt'):
    s_path = f_path[0:-3] + 'h5'
    A.save(s_path)

elif f_path[-2:] == 'h5':
    s_path = f_path

else:
    print "ERROR: only vtk, h5, and txt files accepted!"
    raise ValueError

A.map2saed_screen(zone=zone, pixel_dim=pixel_dim, shell_thickness=shell_thickness)

A.image_screen(A.screen_i, exposure_c, 2, 2, plot=True, file_name=i_path)


