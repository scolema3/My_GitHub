#this script will act as a cask for the conversion of GB_Rots.m
import numpy as np
import scipy.spatial as sp
import re
import math as mt
def grain_rotation(var_dict):

	# form an index lattice to search through for allowable periodic bounds

index = np.indices((2*var_dict['maxmiller']+1,2*var_dict['maxmiller']+1,2*var_dict['maxmiller']+1)).T.reshape(-1,3)
index = np.fliplr(index)
lat_ind = np.dot(index,var_dict['lattice1'])

lat_ind = np.delete(lat_ind,0,0)

#check for uniqueness, remove repeated directions

norm = np.empty((len(lat_ind),1))

for i in range(0,len(lat_ind)):
	norm[i] = np.linalg.norm(lat_ind[i,:])

norm_ind = np.divide(lat_ind,norm)
check_unique = np.full((len(norm_ind),1),True,dtype=bool)

for i in range(0,len(norm_ind)):
	for j in range(0,len(norm_ind)):
		if i == j:
			continue
		elif check_unique[i] == False or check_unique[j] == False:
			continue
		elif np.array_equal(norm_ind[i],norm_ind[j]) == True:
			if norm[i] > norm[j]:
				check_unique[i] = False
			else:
				check_unique[j] = False

check_unique = check_unique.reshape(len(check_unique))
lat_ind_unique = lat_ind[check_unique,:]
norm_ind_unique =norm_ind[check_unique,:]

#find all points that lie within a plane perpendicular to the gb axis

axlat = np.dot(var_dict['axis'],var_dict['lattice1'])
axlat_norm = axlat/np.linalg.norm(axlat)
cos_off = np.full((len(norm_ind_unique),),(1-np.cos(mt.radians(var_dict['deg_tol']))))
plane_true = np.less_equal(np.dot(norm_ind_unique,axlat_norm),cos_off) 

lat_ind_plane = lat_ind_unique[plane_true,:]
norm_ind_plane = norm_ind_unique[plane_true,:]

if len(lat_ind_plane) == 0:
	sys.exit('No points found on plane, try increasing deg_tol.')
elif var_dict['verbose'] == True:
	print 'found %i points within %f degrees of requested plane' % (len(lat_ind_plane),var_dict['deg_tol'])

#pair orthogonal points to form a basis
perp_vec = np.empty((len(norm_ind_plane),1))

#for i in range(0,len(norm_ind_plane)):
perp_vec = np.cross(norm_ind_plane,axlat_norm)
lat_ind_perp = np.cross(lat_ind_plane,axlat_norm)

norm_tree = sp.cKDTree(norm_ind_plane)  # finish implimenting test condition for when a perpendicular vector is within the tolerance angle of a basis vector
perp_query = norm_tree.query_ball_point(perp_vec,np.cos(mt.radians(var_dict['deg_tol'])))
perp_match = np.empty((len(perp_query),1),dtype = bool)

for i in range(0,len(perp_query)):
	if len(perp_query[i]) > 0:
		perp_match[i] = True

perp_match = perp_match.reshape(len(perp_match),)
perp_vec = perp_vec[perp_match,:]
lat_ind_perp = lat_ind_perp[perp_match,:]
lat_ind_comp = lat_ind_plane[perp_match,:]
norm_ind_comp = norm_ind_plane[perp_match,:]
ax_ind = np.full((len(lat_ind_perp),3),axlat)

#check for symmetry before defining allowable

allowable = {'base_x':lat_ind_perp,'base_y':lat_ind_comp,'base_z':ax_ind}

var_dict['lat_or'] = allowable
