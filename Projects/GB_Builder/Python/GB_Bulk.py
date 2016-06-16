#this script will act as a replacement for GB_Bulk.m
import numpy as np

def GB_Bulk(GBorientations,Write_Style,AtomStyle,Archive,Dir_Base,Verbose):

Lattice1=GBorientations[1]['Info']['Lattice'][1]
Basis1=GBorientations[1]['Info']['Basis'][1]
Species1=GBorientations[1]['Info']['Species'][1]
Masses1=GBorientations[1]['Info']['Masses'][1]
Axis1=GBorientations[1]['Info']['Axis'][1]
Direction1=GBorientations[1]['Info']['Direction'][1]
MaxMiller1=GBorientations[1]['Info']['MaxMiller'][1]
Dis_Tol1=GBorientations[1]['Info']['Dis_Tol'][1]
Deg_Tol1=GBorientations[1]['Info']['Deg_Tol'][1]
Symmetry1=GBorientations[1]['Info']['Symmetry'][1]
Sym_Tol1=GBorientations[1]['Info']['Sym_Tol'][1]

Lattice2=GBorientations[2]['Info']['Lattice'][2]
Basis2=GBorientations[2]['Info']['Basis'][2]
Species2=GBorientations[2]['Info']['Species'][2]
Masses2=GBorientations[2]['Info']['Masses'][2]
Axis2=GBorientations[2]['Info']['Axis'][2]
Direction2=GBorientations[2]['Info']['Direction'][2]
MaxMiller2=GBorientations[2]['Info']['MaxMiller'][2]
Dis_Tol2=GBorientations[2]['Info']['Dis_Tol'][2]
Deg_Tol2=GBorientations[2]['Info']['Deg_Tol'][2]
Symmetry2=GBorientations[2]['Info']['Symmetry'][2]
Sym_Tol2=GBorientations[2]['Info']['Sym_Tol'][2]

#find minimum bulk volume
MinVol1=1*10**10
MinVol2=1*10**10

for i in range(0,len(GBorientations)):
	if GBorientations[i]['Lat'][2]['volume'] < MinVol1:
		MinVol1 = GBorientations[i]['Lat'][2]['volume']
		Orient = GBorientations[i]['Lat'][2]
		Lx1 = GBorientations[i]['Lat'][2]['norms'][1]
		Lx1 = GBorientations[i]['Lat'][2]['norms'][2]
		Lx1 = GBorientations[i]['Lat'][2]['norms'][3]

	if GBorientations[i]['Lat'][2]['volume'] < MinVol2:
		MinVol1 = GBorientations[i]['Lat'][2]['volume']
		Orient = GBorientations[i]['Lat'][2]
		Lx1 = GBorientations[i]['Lat'][2]['norms'][1]
		Lx1 = GBorientations[i]['Lat'][2]['norms'][2]
		Lx1 = GBorientations[i]['Lat'][2]['norms'][3]

ComLat = np.linalge.norm(Lattice1-Lattice2) < (1*10**(-9))
ComBasL = len(Basis1) == len(Basis2)
ComSpecL = len(Species1) == len(Species2)
ComMassL = len(Masses1) == len(Masses2)
ComBasis = 0
ComSpec = 0
ComMass = 0

if ComBasL:
	ComBasis = np.linalg.norm(np.asarray(Basis1)-np.asarray(Basis2))<(1*10**(-9))
if ComSpecL:
	ComSpec = Species1 == Species2
if ComMassL:
	ComMass = np.linalg.norm(np.asarray(Masses1)-np.asarray(Masses2))<(1*10**(-9))
if ComLat and ComBasis and ComSpec and ComMass:
	if MinVol1<=MinVol2:
		AtomData = GB_FillRegion(Lattice1,Basis1,Orient1,Lx1,Ly1,Lz1)
		Corners = np.array([0,Lx1,0,Ly1,0,Lz1])
	else:
		AtomData = GB_FillRegion(Lattice2,Basis2,Orient2,Lx2,Ly2,Lz2)
		Corners = np.array([0,Lx2,0,Ly1,0,Lz2]) #possible typo with Ly1
	AtomData = GB_WrapPBC(AtomData,Corners[1::2],.5)
	Header = {'Head1':'# Bulk: 1','Head2':'','Head3':'# Bulk: 1'}
	GB_WriteFiles('Bulk',Write_Style,AtomData,Species1,Masses1,Corners,Header,AtomStyle,Archive,Dir_Base,Verbose)
else:
	AtomData1 = GB_FillRegion(Lattice1,Basis1,Orient1,Lx1,Ly1,Lz1)
	Corners1 = np.array([0,Lx1,0,Ly1,0,Lz1]) 
	AtomData1 = GB_WrapPBC(AtomData1,Corners1[1::2],.5)
	Header1 = {'Head1':'# Bulk: 2\n# Lat 1','Head2':'','Head3':'# Bulk: 2| Lat 1'}
	GB_WriteFiles('Bulk',Write_Style,AtomData1,Species1,Masses1,Corners1,Header1,AtomStyle,Archive,Dir_Base,Verbose)
	
	AtomData2 = GB_FillRegion(Lattice2,Basis2,Orient2,Lx2,Ly2,Lz2)
	Corners2 = np.array([0,Lx2,0,Ly2,0,Lz2]) 
	AtomData2 = GB_WrapPBC(AtomData1,Corners1[1::2],.5)
	Header2 = {'Head1':'# Bulk: 2\n# Lat 1','Head2':'','Head3':'# Bulk: 2| Lat 1'}
	GB_WriteFiles('Bulk',Write_Style,AtomData2,Species2,Masses2,Corners2,Header2,AtomStyle,Archive,Dir_Base,Verbose)
