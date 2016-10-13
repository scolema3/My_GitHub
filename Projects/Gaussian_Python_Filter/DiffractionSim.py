__author__ = 'Adam D. Herron'
__version__ = 1.1
__edited_ = '6 AUG 2015'
__doc__ =''' DiffractionSim - toolset to perform simulated SAED and Kikuchi Diffraction\n
        :argument A = DiffractionSim( [wavelength] ) - see documentation for __init__ method\n
        :return Instance of the DiffractionSim Class\n

        Python/System requirements:
        All system requirements are met by using Anaconda, a free scientific computing python distrubution. If you don't
        want to use Anaconda, the individual requirement are as follows:

        Python 2.7 (developed in 2.7.9)
        HDF5 (required for save() and load() functionality)

        Numpy module
        Scipy module
        Tkinter module
        tkFileDialog module
        matplotlib module
        h5py module

        Basics of how to use this module:\n
        * To get started, import the module into python: import DiffractionSim \n
        * Create an instance of the DiffractionSim class by calling [var] = DiffractionSim( [wavelength] ) \n
        * Use the methods of the DiffractionSim class to manipulate that object. For instance: \n
            - read in VTK or HDF5 data using [var].read_data()
            - filter data in various ways using [var].threshold_data()
            - search for peaks in 3D data using [var].find_peaks3D()
            - map 3D data (or slice of 3D data) to 2D screen using [var].map2saed_screen()
            - view the 2D screen using [var].image_screen()
            - after 2D mapping, search for 2D peaks using [var].find_peaks2D()
            - after 2D peak search view the 2D peaks using view_peaks2D()
            - after finding 3D peak data, call [var].view_peaks3D() to plot a 3D figure
            - save a copy of [var] to manipulate/image/etc later using [var].save()
            - load an old copy of [var] using [var].read_data() or [var].load (read_data() calls load() if HDF5)
        '''

import Tkinter
import tkFileDialog
from numpy import *
import matplotlib
import h5py as h5
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy
import scipy.spatial
import scipy.ndimage as ndimage


rad2degrees = 180 / 3.1415927
degrees2rad = 3.1415927 / 180

def norm(x):
    ''' calculates the Euclidean norm of a 1x3 vector.
    :param x: 1x3 vector of floats or integers
    :return: the Euclidean norm of the vector |x|
    '''
    if len(x) == 3:
        return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])
    else:
        print "Error: norm not defined for lists not 1x3"

def rotate3D(points, axes, angles):
    ''' Function that accepts an (n x 3) numpy array and a series of rotation calls - returns that array rotated\n
    :argument pointsB = rotate3d( pointsA, axes, angles )
    :parameter [points]: (n x 3) numpy array which is to be rotated\n
    :parameter [axes]: (1 x m) list or tuple of integers (0 ,1, or 2) for the (x, y, and z) axies respectively. The
    order specified here is the order in which the rotations will occur.\n
    :parameter [angles]: (1 X m) list or tuple of float angle values (in radians) which correspond to the amount that
    the point array should be rotated about the corresponding axis. Length of [axes] and [angles] must be the same.\n
    :returns [ndarray]: (n x 3) array the same size as [points] which represents the rotated point cloud.\n
    - rotate3D constructs a 3D rotation matrix based on the axes and angles specified, does left-side matrix
    multiplication on the transpose of a copy of the point cloud and then re-transposes the point cloud and returns it.\n
    '''
    n_rot = len(axes)
    if n_rot == 0 or len(angles) != n_rot:
        print "ERROR: number of axes and number of angles must be the same"
        return 0
    r = array([[1,0,0],[0,1,0],[0,0,1]])
    for i in range(0,n_rot):
        c = cos(angles[i])
        s = sin(angles[i])
        if axes[i] == 0:
            r = dot(array([[1,0,0],[0,c,-s],[0,s,c]]), r)
            print "--- Rotated %f Degrees About the x Axis ---" % (angles[i]*rad2degrees)
        elif axes[i] == 1:
            r = dot(array([[c,0,s],[0,1,0],[-s,0,c]]), r)
            print "--- Rotated %f Degrees About the y Axis ---" % (angles[i]*rad2degrees)
        elif axes[i] == 2:
            r = dot(array([[c,-s,0],[s,c,0],[0,0,1]]), r)
            print "--- Rotated %f Degrees About the z Axis ---" % (angles[i]*rad2degrees)
        else:
            print "ERROR: axis must be 0,1, or 2 for 3D rotation"
    rotated_points = dot(r, copy(points).transpose()).transpose()
    return rotated_points

class DiffractionSim:
    __doc__ = 'Class for handling and visualizing diffraction simulation data'
    def __init__(self, wavelength=0.0251):   # Construct a DiffractionSim Object
        ''' Creates an instance of the DiffractionSim class\n
        :argument A = DiffractionSim( [wavelength] )\n
        :parameter [wavelength]: characteristic wavelength of the simulated rediation in the system. It is used to
        determine the radius of the Ewald Sphere. It is a public member and can be changed at any time\n
        '''
        self.__version__ = __version__
        self.__doc__ = __doc__
        if wavelength == 0:
            print "ERROR: wavelength of incident radiation cannot be 0 - bad things happens!"
            raise ValueError
            return 0
        self.wavelength = wavelength
        self.reciprocal_points = []
        self.reciprocal_spacing = -1
        self.intensity = []
        self.file_header = {}
        self.n_k = -1
        self.k_mag = []
        self.zone = array([1, 0, 0])
        self.screen_distance = -1
        self.dpm = -1
        self.p_sweep = -1
        self.n_p = -1
        self.a_sweep = -1
        self.n_a = -1
        self.search_radius = -1
        self.integral_radius = -1
        self.peaks = []
        self.screen_x = []
        self.screen_y = []
        self.screen_z = []
        self.screen_i = []
        self.n_atoms = -1
        self.complex = False
        print "Created DiffractionSim Object with characteristic wavelength", wavelength, "Angstroms"

    def read_data(self, file_name='FileDialogGUI', zone=None, n_atoms=None):
        '''Reads K-points and intensity values from data file or loads a previous DiffractionSim state from an hdf5 file\n
        At this point, read_data() only accepts vtk files output from Shawn Coleman's LAMMPS code fix_time_ave_saed,
        txt files output by Eric Homer's intSAED LAMMPS code, and HDF5 files created using the save() method\n
        :argument read_data( [file_name] )\n
        :parameter [file_name]: name of or path to desired data file - the default value opens a TK file dialog window\n
        :parameter [zone]: zone associated with the file to read in\n
        :parameter [n_atoms]: number of atoms associated with this simulation\n
        :returns [None]: stores data in self.reciprocal_points, self.intensity, self.n_k, and self.k_mag\n
        '''
        root = Tkinter.Tk()
        root.withdraw()
        if file_name == 'FileDialogGUI':
            print "--- Opened File Dialog Gui ---"
            file_name = tkFileDialog.askopenfilename()
            if file_name == '':
                print "No file specified --- returned None"
                return None
        fin = open(file_name, 'r')
        if (file_name[-3:] == 'vtk'):
            if fin:
                print "Opened", file_name
                print "Reading in intensity values . . ."
                threshold = 0
                head = []
                for i in range(0,4):
                    head.append(fin.readline())
                d = fin.readline().split()
                aspect_ratio = fin.readline().split()
                self.file_header['head'] = head
                self.file_header['dimensions'] = d
                self.file_header['aspect_ratio'] = aspect_ratio
                self.file_header['origin'] = fin.readline().split()
                self.file_header['point_data'] = int(fin.readline().split()[1])
                self.file_header['scalars'] = fin.readline().split()
                self.file_header['lookup_table'] = fin.readline().split()
                intensity = loadtxt(file_name, dtype=float, skiprows=10, usecols=[0])
                pos = empty((int(d[1]) * int(d[2]) * int(d[3]), 3),float)
                n = 0

                ## Shawn's Way ##
                # for k in range(1,int(d[3])+1):
                #     for j in range(1,int(d[2])+1):
                #         for i in range(1,int(d[1])+1):
                #             pos[n,:] = [i, j, k]
                #             n += 1
                # pos = pos - dot(ones([len(pos), 1]),
                #                array([[(int(d[1])+1)/2.0, (int(d[2])+1)/2.0, (int(d[3])+1)/2.0]]))
                # aspect_ratio = array([float(aspect_ratio[1]), float(aspect_ratio[2]), float(aspect_ratio[3])])
                # for i in range(0,3):
                #     pos[:,i] = pos[:,i] * aspect_ratio[i]

                ## Try new way ##
                origin = self.file_header['origin']
                origin = [float(origin[1]), float(origin[2]), float(origin[3])]
                aspect_ratio = array([float(aspect_ratio[1]), float(aspect_ratio[2]), float(aspect_ratio[3])])
                for k in range(0,int(d[3])):
                    for j in range(0,int(d[2])):
                        for i in range(0,int(d[1])):
                            pos[n,:] = [origin[0]+i*aspect_ratio[0], origin[1]+j*aspect_ratio[1], origin[2]+k*aspect_ratio[2]]
                            n += 1

                test = intensity > threshold
                # intensiry = intensity[test]
                intensity = intensity[test] * (aspect_ratio[0] * aspect_ratio[1] * aspect_ratio [2] / 3e-07)
                self.intensity = intensity
                self.reciprocal_points = pos[test, :]
                self.n_k = len(self.reciprocal_points)
                print "n: ", n
                print "file_n: ", self.file_header['point_data']
                self.k_mag = sqrt(self.reciprocal_points[:,0]**2 + \
                                  self.reciprocal_points[:,1]**2 + self.reciprocal_points[:,2]**2)
                print "--- Finished reading vtk file ---"
                self.zone = zone
                self.n_atoms = n_atoms
                self.reciprocal_spacing = aspect_ratio
            elif not fin:
                print file_name, " did not open correctly"
        elif (file_name[-2:] == 'h5'):
            if fin:
                self.load(file_name)
                if zone is not None:
                    self.zone = zone
                if n_atoms is not None:
                    self.n_atoms = n_atoms
            else:
                print file_name, " did not open correctly"
        elif (file_name[-3:] == 'txt'):
            if fin:
                print "Opened", file_name
                print "Reading in data . . ."
                header = []
                header.append(fin.readline())
                header.append(fin.readline())
                header.append(fin.readline().split())
                header.append(fin.readline().split())
                if header[2][-1][-2] == '5':
                    self.complex = False
                    self.n_k = int(header[3][1])
                    self.reciprocal_points = loadtxt(file_name, skiprows=4, usecols=[1, 2, 3])
                    theta = loadtxt(file_name, skiprows=4, usecols=[4])
                    self.k_mag = 2 * sin(theta) / self.wavelength
                    self.intensity = loadtxt(file_name, skiprows=4, usecols=[5])
                elif header[2][-1][-2] == '6':
                    self.complex = True
                    self.n_k = int(header[3][1])
                    self.reciprocal_points = loadtxt(file_name, skiprows=4, usecols=[1, 2, 3])
                    theta = loadtxt(file_name, skiprows=4, usecols=[4])
                    self.k_mag = 2 * sin(theta) / self.wavelength
                    intensity = loadtxt(file_name, skiprows=4, usecols=[5, 6])
                    self.intensity = 1j * intensity[:,1] + intensity[:,0]
                else:
                    print "Unable to read file - incorrect header format for txt file type"
                print "--- Finished Reading File ---"
            else:
                print file_name, " did not open correctly"
        else:
            print "ERROR: read_data() only accepts txt, vtk and specific hdf5 files at present"
            exit(0)

    def threshold_data(self, threshold_magnitude=-1, k_min=-1, k_max=1.1, threshold=None):
        '''Selects points from self.recipocal_points via threshold criterion and returns an (n x 5) numpy array\n
        :argument k = threshold_data([threshold_magnitude=5], [k_min=0.05], [k_max=1.1], [threshold])\n
        :parameter [threshold magnitude]: orders of magnitude of intensity below which all points are ignored\n
        If a negative value is passed as threshold magnitude, the threshold will be zero.\n
        :parameter [threshold]: overides the threshold_magnitude parameter by letting the user specify the
        actaul threshold value itself\n
        :parameter [k_min]: distance from origin below which all points are ignored\n
        :parameter [k_max]: distance from origin above which all points are ignored\n
        :returns [k]: one (n x 5) numpy array which contains n rows: [reciprocal_points index, h, k, l, intensity]\n
        '''
        if threshold is not None:
            t = threshold
        else:
            intensity_max = max(self.intensity)
            if threshold_magnitude > 0:
                t = intensity_max * (10**(-threshold_magnitude))
            elif threshold_magnitude <= 0:  # negative values clipped off by default
                t = 0
        test = (self.intensity > t) * (self.k_mag >= k_min) * (self.k_mag <= k_max)

        k = ones([sum(test),5])
        k[:,0] = range(0,sum(test))
        k[:,1:4] = self.reciprocal_points[test,:]
        if self.complex:
            k = k.astype(complex_)
        k[:,4] = self.intensity[test]
        return k

    def _kikuchi_screen(self, zone=array([1,0,0]), screen_distance=3,
                       dpm=27.16244, p_sweep=3.1416, a_sweep=3.1416):
        ''' Generate 3 arrays of pixel location indicies that describe a view screen
        :argument screen_vectors([self], [zone=[1,0,0]], [screen_distance=3], [dpm=256], [p_sweep=pi], [a_sweep=pi])\n
        :parameter [zone]: direction of zone axis pointing toward the center of the viewscreen\n
        :parameter [screen_distance]: distance ( in mm ) from origin to center of view screen\n
        :parameter [dpm]: approximate resolutions of the viewscreen ( in dots-per-millimeter squared )\n
        :parameter [p_sweep]: total polar (horizontal) angle capture of view screen\n
        :parameter [a_sweep]: total azimuthal (vertical) angle capture of view screen\n
        :returns None: generates self.screenx/y/z (m x n) arrays of pixel location coordinates --> bmp style layout\n
        '''
        print "Generating screen vectors. . ."
        if type(zone) == list:
            zone = array(zone)
        if type(zone) == ndarray and len(zone == 3):
            self.zone = zone
            self.screen_distance = screen_distance
            self.dpm = dpm
            self.p_sweep = p_sweep
            self.a_sweep = a_sweep
            ilambda = 1 / self.wavelength
            n_p = int(round(dpm * p_sweep * screen_distance))
            n_a = int(round(dpm * a_sweep * screen_distance))
            self.n_p = n_p
            self.n_a = n_a
            delta_p = p_sweep/n_p
            delta_a = a_sweep/n_a
            unit_zone = zone/norm(zone)
            c_angle = [arccos(unit_zone[0]), arcsin(unit_zone[2])]
            if unit_zone[1] < 0:
                c_angle[0] = -c_angle[0]
            p_angles = linspace(c_angle[0]-(p_sweep/2), c_angle[0]+(p_sweep/2), num=n_p)
            a_angles = linspace(c_angle[1]+(a_sweep/2), c_angle[1]-(a_sweep/2), num=n_a)
            x = ilambda * outer(cos(a_angles), cos(p_angles))
            y = ilambda * outer(cos(a_angles), sin(p_angles))
            z = ilambda * outer(sin(a_angles), ones(len(p_angles)))
            print "--- Finished generating screen vectors ---", n_a, "x", n_p, "bitmap array"
            self.screen_x, self.screen_y, self.screen_z, self.screen_i = x, y, z, full((len(x), len(x[0])), 0.0, float32)
        else:
            print "ERROR: zone must be a list or ndarray of length 3"

    def _saed_screen(self, k_max=1.1, pixel_dim=256):
        ''' Generates a screen to map SAED pattern data onto\n
        :argument self._saed_screen(k_max, pixel_dim)\n
        :param [k_max]: width/height of square view screen (in reciprocal space)\n
        :param [pixel_dim]: number of pixels wide/hgih that the screen is\n
        :return: None - stores pixel locations/intensities in self.screen_x/y/z/i\n
        '''
        values = linspace(-k_max, k_max, pixel_dim, dtype=float)
        oneslist = ones((pixel_dim, 1), dtype=float)
        self.screen_y, self.screen_x = meshgrid(values, values)
        self.screen_z = zeros((pixel_dim, pixel_dim))
        self.screen_i = copy(self.screen_z)

    def image_screen(self, screen=None, exposure_c=2, sigma=2, blur_steps=2, plot=True, file_name=None):
        ''' Smooths and plots (or saves) data stored in an nxm image-like array (defaults to self.screen_i)\n
        :param [screen]: nxm image-like array to show/save - defaults to self.screen_\n
        :param [exposure_c]: exposure constant (1 currently tunes to a good value for Ni-Co-Aluminum)\n
        :param [sigma]: sigma for the Gaussian Convolution\n
        :param [blur_steps]: number of sequential convolutions to perform on the image\n
        :param [plot]: True opens a matplotlib window and plots the image for the user to view\n
        :param [file_name]: if specified, an image is saved to disk under this name (in current working directory) or path\n
        :return: None
        '''
        if screen is None:
            screen2 = real((self.screen_i))
        else:
            screen2 = real(screen)
        for i in range(0,blur_steps):
            screen2 = ndimage.filters.gaussian_filter(screen2, sigma=sigma)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        vmin = 0
        if self.n_atoms != -1:
            vmax = 1000 / (log(self.n_atoms) * exposure_c)
        else:
            print "WARNING: self.n_atoms (number of atoms in simulation) is not specified - exposure not calibrated!"
            vmax = amax(screen2) / exposure_c
        if file_name is not None:
            matplotlib.image.imsave(file_name, screen2, vmin=vmin, vmax=vmax, cmap='gray')
            print "--- Saved screen image --- ", file_name
        if plot:
            image = ax.imshow(screen2.astype(float32))
            image.set_cmap('gray')
            # image.set_cmap('hot')
            # image.set_cmap('bone')
            image.set_clim(vmin=vmin, vmax=vmax)
            # plt.colorbar(image)
            print "--- Displayed screen image ---"
            plt.show()

    def find_peaks3D(self, search_radius=0.09, integral_radius=0.09, threshold_magnitude=-1,
                     threshold=None, k_min=0, k_max=1.1):
        """ method to locate intensity peaks and calculate a numerical integral intensity for each peak\n
        :argument self.find_peaks([search_radius], [integral_radius])\n
        :parameter [search_radius]: distance from each point that the algorithm will search for maxima\n
        :parameter [integral_radius]: radius from maxima for inclusion of neighbor inetensity in integral calculation\n
        :parameter [k_min]: distance from origin to filter out data (removes center peak)\n
        :parameter [k_max]: distance from origin to filter out data (removes outside fringe)\n
        :returns None: stores data in self.peaks as an m x 5 array with peak data:
        m rows of:  [reciprocal_point index, h, k, l, integral_intensity]\n
        """
        print "Searching for intensity peaks . . ."
        self.search_radius = search_radius
        self.integral_radius = integral_radius
        if len(self.reciprocal_points) == 0:
            print "ERROR: DiffractionSim object has no k points loaded - load them using read_data() method"
            return False
        else:
            test = (self.k_mag >= k_min) * (self.k_mag <= k_max)
            n_k = sum(test)
            k = empty((n_k,5),float_)
            k[:,0] = array(range(0,n_k))
            k[:,1:4] = self.reciprocal_points[test,:]
            k[:,4] = absolute(self.intensity[test])
            pool = k[:,4] > threshold
            pool_size = copy(sum(pool))
            maxima = zeros(n_k,dtype=bool)
            srq = search_radius**2
            irq = integral_radius**2
            rt3 = sqrt(3)
            while sum(pool) > 0:
                poi = ((k[pool,0])[k[pool,4] == max(k[pool,4])])[0]
                rt3_test = abs(k[:,1] - full((len(k[:,1])), k[poi,1])) <= rt3
                rt3_test *= abs(k[:,2] - full((len(k[:,2])), k[poi,2])) <= rt3
                rt3_test *= abs(k[:,3] - full((len(k[:,3])), k[poi,3])) <= rt3
                dist_sq = ones((sum(rt3_test),3))
                dist_sq[:,0] = k[rt3_test,1] - dist_sq[:,0]*k[poi,1]
                dist_sq[:,1] = k[rt3_test,2] - dist_sq[:,1]*k[poi,2]
                dist_sq[:,2] = k[rt3_test,3] - dist_sq[:,2]*k[poi,3]
                dist_sq = dist_sq[:,0]**2 + dist_sq[:,1]**2 +dist_sq[:,2]**2
                neighbors = (k[rt3_test,:])[dist_sq < srq,:]
                for index in neighbors[:,0]:
                    if pool[index]: pool[index] = False
                max_neighbors = neighbors[neighbors[:,4] == max(neighbors[:,4]), :]
                for index in max_neighbors[:,0]:
                    if not maxima[index]:
                        if index != poi:
                            pool[index] = True
                        elif index == poi:
                            maxima[index] = True
                # print "%d points left of %d total" % (sum(pool), pool_size)
            local_maxima = k[maxima,:]
            for index in local_maxima[:, 0]:
                dist_sq = ones((n_k,3), dtype=float)
                dist_sq[:,0] = k[:,1] - dist_sq[:,0]*k[index,1]
                dist_sq[:,1] = k[:,2] - dist_sq[:,1]*k[index,2]
                dist_sq[:,2] = k[:,3] - dist_sq[:,2]*k[index,3]
                dist_sq = dist_sq[:,0]**2 + dist_sq[:,1]**2 +dist_sq[:,2]**2
                neighbors = dist_sq < irq
                local_maxima[local_maxima[:,0]==index,4] = sum(k[neighbors,4])
            print "--- {n_peaks} Peaks Located ---".format(n_peaks=len(local_maxima[:,0]))
            self.peaks = empty((len(local_maxima[:,0]), 5), float)
            self.peaks[:,:] = local_maxima

    def find_peaks2D(self, sigma=2, blur_steps=4):
        ''' method to find peaks in the 2D self.screen_i after gaussian blurring\n
        argument: \t\t peaks = self.find_peaks2D(sigma, blur_steps)\n
        :param [sigma]: float, standard deviation used in the gaussian blur\n
        :param [blur_steps]: integer >= 0, number of times the image is seqentially blurred before the maximum filter is applied\n
        :return [peaks]: a boolean array the same size as self.screen_i which is true at local maxima\n
        '''
        screen = copy(self.screen_i)
        if blur_steps < 0 or type(blur_steps) is not int:
            print "ERROR: blur_steps must be a non-negative integer!"
            raise TypeError
        for i in range(0,blur_steps):
            screen = ndimage.filters.gaussian_filter(screen, sigma, mode='nearest')
        neighborhood = array([[0,1,1,1,0],[1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1],[0,1,1,1,0]])
        # neighborhood = array([[0,1,0],[1,1,1],[0,1,0]])
        # neighborhood = array([[0,0,0],[0,1,0],[0,0,0]])
        local_max = ndimage.maximum_filter(screen, footprint=neighborhood) == screen
        background = screen <= 0
        eroded_background = ndimage.binary_erosion(background, structure=neighborhood, border_value=1)
        peaks = local_max - eroded_background
        return peaks

    def view_peaks3D(self, zone=None):
        ''' method for plotting intensity peaks in 3D from reciprocal_point data\n
        argument: self.view_peaks3D(zone)
        :param [zone]: the direction of the view axis pointing from the origin toward the user's viewpoint.
        :return None: opens a plot for the user to view and manipulate
        '''
        if len(self.peaks) <= 0:
            print "ERROR: No peaks to view! Call self.find_peaks3D before calling self.view_peaks3D"
            raise AttributeError

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        s = log(self.peaks[:, 4])
        s = 800 * s / max(s)
        s = s - full(s.size, (min(s) - 10))

        ax.scatter(self.peaks[:,1], self.peaks[:,2], self.peaks[:,3], c=s, s=s, marker='o', edgecolors='none', alpha=0.75)
        ax.set_xlabel('X Axis')
        ax.set_ylabel('Y Axis')
        ax.set_zlabel('Z Axis')
        ax.set_xlim([-1.1, 1.1])
        ax.set_ylim([-1.1, 1.1])
        ax.set_zlim([-1.1, 1.1])
        ax.set_aspect('equal', 'box')

        if zone is None:
            zone2 = self.zone
        else:
            zone2 = zone
        if type(zone2) == list:
            zone2 = array(zone2)
        if zone2 is None or (zone2[0] == 0 and zone2[1] == 0 and zone2[2] == 0):
            uz = array([1, 0, 0])
        else:
            uz = zone2 / norm(zone2)

        ax.plot([-uz[0], uz[0]], [-uz[1], uz[1]], [-uz[2], uz[2]], color='black', lw=1)

        if uz[2] == 0:
            elev = 0
        else:
            elev = arcsin(uz[2]) * (180 / 3.1415927)
        if uz[1] == 0:
            azim = 0
        else:
            azim = arccos(uz[0] / cos(elev * (3.1415927 / 180)))  * (180 / 3.1415927)
            if arcsin(uz[1]) < 0:
                azim = -azim
        print "Polar Angle: {p_angle} degrees\nAzimuthal Angle: {a_angle} degrees".format(p_angle=azim, a_angle=elev)
        ax.view_init(elev=elev, azim=azim)

        plt.show()

    def view_peaks2D(self, screen=None, peaks=None, exposure_c=2, sigma=2, blur_steps=2, plot=True, file_name=None):
        ''' method to view only those regions of a 2D image which are deemed 'peaks'
        :param screen: n x n numpy array containing a 2D image (defaults to self.screen_i)
        :param peaks: n x n boolean array containing peak locations in image (defaults to a new peak search)
        :param exposure_c: exposure level of final image
        :param sigma: standard deviation of Gaussian blurring
        :param blur_steps: interger number of Gaussian blurrings
        :param plot: Boolean that plots the final image for user if True
        :param file_name: if PATH (as a string) is specified, final image is saved to PATH
        :returns: None - saves and/or plots image
        '''
        if peaks is None:
            peaks2 = self.find_peaks2D(sigma)
        else:
            peaks2 = peaks

        self.image_screen(peaks2 * self.screen_i, exposure_c, sigma, blur_steps, plot, file_name)

    def save(self, file_name='FileDialogGUI'):
        ''' Saves the current state of the DiffractionSim object to an hdf5 file for later use - load it using load()\n
        :argument self.save( [file_name] )\n
        :parameter [ file_name ]: If a file name is passed, this name is used and same-name files are overwirtten\n
        Otherwise, a TK file dialog box is opened for the user\n
        '''
        root = Tkinter.Tk()
        root.withdraw()
        if file_name == 'FileDialogGUI':
            print "--- Opened File Dialog GUI ---"
            file_name = tkFileDialog.asksaveasfilename()
            if len(file_name) <= 0:
                return 0
            if file_name[-3:] != '.h5':
                file_name = file_name + '.h5'
        file = h5.File(file_name, 'w')
        print "Saving data to ", file_name, ". . ."
        file.create_group('DiffractionSim_Data')
        for i in dir(self):
            if not callable(getattr(self, i)) and not i[0:2] == '__' and not i[0:4] == 'file':
                file['/DiffractionSim_Data'].create_dataset(name=i, data=getattr(self, i))
        file.close()
        print "--- Finished Save ---"

    def load(self, file_name='FileDialogGUI'):
        '''Loads a previously saved version of the DiffractionSim object saved using self.save().\n
        :argument self.load( [file_name] )\n
        :parameter [file_name]: If a file name is passed, this name is used. Otherwise, a TK file dialog box is opened
        for the user\n
        '''
        root = Tkinter.Tk()
        root.withdraw()
        if file_name == 'FileDialogGUI':
            file_name = tkFileDialog.askopenfilename()
        file = h5.File(file_name, 'r')
        print "Loading data from ", file_name, ". . ."
        if file.__contains__('/DiffractionSim_Data'):
            for i in file['/DiffractionSim_Data']:
                setattr(self, i, file['/DiffractionSim_Data/' + i][...])
                print "Loaded DiffractionSim_Data/" + i
            print "--- Finished loading data ---"
            file.close()
            if len(self.intensity) == 0 and len(self.reciprocal_points) == 0:
                print "WARNING: The data file contianed no reciprocal points or intensity data"
        else:
            print "ERROR: File does not contain DiffractionSim_Data"

    def map2saed_screen(self,k=None, zone=None, threshold_magnitude=-1, k_min=0, k_max=1.1, pixel_dim=256,
                        shell_thickness=None, threshold=None):
        ''' method that maps 3D reciprocal_point data onto a 2D saed_screen for image manipulation and viewing\n
        :param [k]: points to map onto the screen (nx5) array of row [index, h, k, l, intensity]\n
        :param [zone]: zone axis to view - defaults to self.zone set or loaded during self.read_data()\n
        :param [threshold_magnitude]: number of orders of magnitude of intensity below the maximu below which all data
        is filtered out\n
        :param [k_min]: minimum distance from origin below which all data is filtered out\n
        :param [k_max]: maximum distance from origin above which all data is filtered out\n
        :param [pixel_dim]: number of pixels high/wide that the square screen is\n
        :param [shell_thickness]: if specified, this limits the points mapped to only those within a swath of this width
        surrounding the ewald sphere\n
        :param [threshold]: can be used instead of threshold_magnitude - points below this intensity are filtered out\n
        :return [None]: mapped points are stored in self.screen_i, a pixel_dim x pixel_dim array of floats
        '''
        print "Mapping SAED data to 2-D screen . . ."
        self._saed_screen(k_max, pixel_dim)
        if k is None:
            k2 = self.threshold_data(threshold_magnitude,k_min,k_max,threshold)
        else:
            k2 = k
        if zone is None:
            zone2 = self.zone
        else:
            zone2 = zone
        if type(zone2) == list:
            zone2 = array(zone2)
        if zone2 is None or (zone2[0] == 0 and zone2[1] == 0 and zone2[2] == 0):
            print "ERROR: Zone is None or [0,0,0]!"
            uz = array([1, 0, 0])
        else:
            uz = zone2 / norm(zone2)
            self.zone = zone2
        if shell_thickness is not None:
            ec = uz / self.wavelength
            n_k2 = len(k2[:,0])
            # edsq = k2[:,1:4] - dot(ones((len(k2[:,0]),3), float), array([[ec[0], 0, 0],[0,ec[1],0],[0,0,ec[2]]]))
            edsq = k2[:,1:4] - concatenate((full((n_k2,1),ec[0]), full((n_k2,1),ec[1]), full((n_k2,1),ec[2])), axis=1)
            edsq = edsq[:,0]**2 + edsq[:,1]**2 + edsq[:,2]**2
            t = (edsq >= ((1/self.wavelength) - (shell_thickness / 2))**2) *\
                   (edsq <= ((1/self.wavelength) + (shell_thickness / 2))**2)
            k2 = k2[t,:]
        if uz[2] == 0:
            elev = 0
        else:
            elev = arcsin(uz[2])
        if uz[0] == 0:
            polar = ((3.1415927 / 2)* sign(uz[1]))
        else:
            polar = arctan(uz[1]/uz[0])
        azim = (3.1415927 / 2) - elev
        k2[:,1:4] = rotate3D(k2[:,1:4],[2,1],[-polar, -azim])
        print "Please Wait . . . "
        sqrt2res = sqrt(2) * 2 * k_max / (pixel_dim - 1)
        screen_shape = self.screen_i.shape
        x = ravel(self.screen_x)
        y = ravel(self.screen_y)
        screen_i = ravel(self.screen_i)
        # tree = scipy.spatial.KDTree(zip(x,y))
        tree = scipy.spatial.cKDTree(zip(x,y))
        if self.complex:
            screen_i = screen_i.astype(complex_)
            d, index = tree.query(real(k2[:,1:3]))
        else:
            d, index = tree.query(k2[:,1:3])
        for i in range(0,len(index)):
            screen_i[index[i]] += k2[i,4]
        self.screen_i = reshape(screen_i,screen_shape)
        print "--- Finished Mapping ---"
    def hr_transform(self, screen=None, shift=0):
        ''' Method that tests the ability of this class to generate high-resolution images by performing a 2D fast
        fourier transform and then a phase shift on a complex 2D image\n
        :param [screen]: an n x n complex numpy array that contains the image to be transformed (defaults to self.screen_i)
        :param [shift]: phase shift angle (in radians)
        :returns [screen_2]: transformed n x n complex numpy array
        '''
        if screen is None:
            screen_2 = copy(self.screen_i).astype(complex_)
        else:
            screen_2 = copy(screen).astype(complex_)
        screen_2 = fft.fft2(screen_2)
        screen_2 *= exp(1j * shift)
        return screen_2
