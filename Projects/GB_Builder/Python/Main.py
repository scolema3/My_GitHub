# this script will act as the master for GB_builder_python, Variable control and function flow will be handled here

#read variables from a config file
import grain_rot
import numpy as np
import scipy.spatial as sp
import re
import math as mt
var_dict = {}

f = open('example_config.txt','r')
datastr = f.read().strip().split(';\n')
f.close()

for i in range(0,len(datastr)):
	data_entry = str(datastr[i])
	data_entry = re.split('\n| = | =\n',data_entry)
	while '#' in data_entry[0] or data_entry[0]=='':
		data_entry.remove(data_entry[0])
	data_value = ','
	data_value = data_value.join(data_entry[1:])
	var_dict[data_entry[0]] = data_value

#convert stringed variables to proper data format

if len(var_dict['lattice2']) == 0:
	var_dict['lattice2'] = var_dict['lattice1']

arrays = {'lattice1':3,'lattice2':3,'basis':4,'orientation':3}
vectors = ['axis','normslabdim','translate','angle_range','write_style']
bools = ['symmetry','verbose','fully_periodic','stoich','forcestoich','summary','archive']
floats = ['masses','sym_tol','dis_tol','deg_tol','strain_tol','maxarea','overlap_tol','vacuum','gbregion','pbc_overlap_tol']

var_dict['maxmiller']=int(var_dict['maxmiller'])

for key,value in arrays.items():
	tmp = np.fromstring(var_dict[key],sep=',')
	tmp = tmp.reshape(len(tmp)/value,value)
	var_dict[key] = tmp

for i in range(len(bools)):
	var_dict[bools[i]] = var_dict[bools[i]]=='True'

for i in range(len(floats)):
	var_dict[floats[i]] = float(var_dict[floats[i]])

for i in range(len(vectors)):
	tmp = np.fromstring(var_dict[vectors[i]],sep=',')
	var_dict[vectors[i]] = tmp

#if var_dict['ngbs'] == 'NULL':
#	var_dict['ngbs'] = len(var_dict['orientation'])/3
#else:
#	var_dict['ngbs'] = int(var_dict['ngbs'])

var_check = grain_rot.grain_rotation(var_dict)
