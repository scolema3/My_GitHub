import numpy
import math
import os
import uuid
import mdcs


user='admin'
pswd='admin'
host='http://127.0.0.1:8000/'
template_name='record-calculation-grain-boundary'

user='scolema3'
pswd='MDCS1pswd'
host='https://mdcs1.nist.gov'
template_name='record-calculation-grain-boundary'


# Sigma 3 = 1,1,0
# Sigma 9,11 = 1, -1, 0
# Sigma 5, 13 = 0, 0, 1
sigmaval=0
tilt=numpy.array([0,0,1])
stgbval=001

recordtemplate='/Users/Coleman/Research/MDCS/Example_Data/Curate/template-record-tschopp-Al-data.xml'
infofile='/Users/Coleman/Research/MDCS/Example_Data/Curate/lammps15May2015_Al/info.stgb'+'001'+'.txt'
datadir='/Users/Coleman/Research/MDCS/Example_Data/Curate/lammps15May2015_Al'


#Read in orientations from info file
g = open(infofile,'r')
g.readline()
g.readline()
g.readline()
g.readline()
g.readline()
g.readline()
for i in range(0,1000):

    line = g.readline()
    data=line.split()

    if len(data)<1:
	break

    index=data[0]
    angle=data[1]
    Al_GBE=data[2]
    Cu_GBE=data[3]

    ori1z=tilt
    ori2y=numpy.array([float(data[4]), float(data[5]), float(data[6])])
    ori2z=tilt
    ori1y=numpy.array([float(data[7]), float(data[8]), float(data[9])])
    ori2x=numpy.cross(ori2y,ori2z)
    ori1x=numpy.cross(ori1y,ori1z)

    xsize1=numpy.linalg.norm(ori1x)
    xsize2=numpy.linalg.norm(ori2x)
    if xsize1<=xsize2:
        xsize=xsize1
    else:
        xsize=xsize2

    zsize1=numpy.linalg.norm(ori1z)
    zsize2=numpy.linalg.norm(ori2z)
    if zsize1<=zsize2:
        zsize=zsize1
    else:
        zsize=zsize2

    x=numpy.array([1,0,0])
    y=numpy.array([0,1,0])
    z=numpy.array([0,0,1])

    R1=numpy.matrix([[float(numpy.dot(x,ori1x)),float(numpy.dot(y,ori1x)),float(numpy.dot(z,ori1x))],[float(numpy.dot(x,ori1y)),float(numpy.dot(y,ori1y)),float(numpy.dot(z,ori1y))],[float(numpy.dot(x,ori1z)),float(numpy.dot(y,ori1z)),float(numpy.dot(z,ori1z))]])
    R2=numpy.matrix([[float(numpy.dot(x,ori2x)),float(numpy.dot(y,ori2x)),float(numpy.dot(z,ori2x))],[float(numpy.dot(x,ori2y)),float(numpy.dot(y,ori2y)),float(numpy.dot(z,ori2y))],[float(numpy.dot(x,ori2z)),float(numpy.dot(y,ori2z)),float(numpy.dot(z,ori2z))]])


    # Read in data file information
    datafile=datadir+'/lammps15May15.Al.stgb'+'001'+'_'+repr(int(index))+'.dat'
    d = open(datafile,'r')
    d.readline()
    d.readline()
    line2=d.readline()
    data2=line2.split()
    natoms=int(data2[0])
    d.readline()
    d.readline()
    line2=d.readline()
    data2=line2.split()
    lx=float(data2[1])-float(data2[0])
    comx=float(data2[0])+lx/2
    line2=d.readline()
    data2=line2.split()
    com1y=float(data2[0])/2
    com2y=float(data2[1])/2
    ly=float(data2[1])-float(data2[0])
    line2=d.readline()
    data2=line2.split()
    lz=float(data2[1])-float(data2[0])
    comz=float(data2[0])+lz/2
    area=lz*lx
    d.close()


    UIDGrain1=uuid.uuid4()
    UIDGrain2=uuid.uuid4()
    UIDRecord=uuid.uuid4()
    UIDGB=uuid.uuid4()
    GBname='Al.stgb'+'001'+'_'+repr(int(index))

    print 'uploading blob'
    url = mdcs.blob.upload(datafile,host,user,pswd)
    print url


    replacements = {'HERE_CALCULATION_KEY_HERE':str(UIDRecord),
                    'HERE_GB_KEY_HERE':str(UIDGB),
                    'HERE_GB_NAME_HERE':GBname,
                    'HERE_XDIM_HERE':str(lx),
                    'HERE_YDIM_HERE':str(ly),
                    'HERE_ZDIM_HERE':str(lz),
                    'HERE_N_ATOMS_TOTAL_HERE':str(natoms),
                    'HERE_GB_AREA_HERE':str(area),
                    'HERE_GRAIN1_KEY_HERE':str(UIDGrain1),
                    'HERE_ORIENT1_X1_HERE':str(int(R1[0,0])),
                    'HERE_ORIENT1_X2_HERE':str(int(R1[0,1])),
                    'HERE_ORIENT1_X3_HERE':str(int(R1[0,2])),
                    'HERE_ORIENT1_Y1_HERE':str(int(R1[1,0])),
                    'HERE_ORIENT1_Y2_HERE':str(int(R1[1,1])),
                    'HERE_ORIENT1_Y3_HERE':str(int(R1[1,2])),
                    'HERE_ORIENT1_Z1_HERE':str(int(R1[2,0])),
                    'HERE_ORIENT1_Z2_HERE':str(int(R1[2,1])),
                    'HERE_ORIENT1_Z3_HERE':str(int(R1[2,2])),
                    'HERE_COM_X_HERE':str(comx),
                    'HERE_COM1_Y_HERE':str(com1y),
                    'HERE_COM2_Y_HERE':str(com2y),
                    'HERE_COM_Z_HERE':str(comz),
                    'HERE_GRAIN2_KEY_HERE':str(UIDGrain2),
                    'HERE_ORIENT2_X1_HERE':str(int(R2[0,0])),
                    'HERE_ORIENT2_X2_HERE':str(int(R2[0,1])),
                    'HERE_ORIENT2_X3_HERE':str(int(R2[0,2])),
                    'HERE_ORIENT2_Y1_HERE':str(int(R2[1,0])),
                    'HERE_ORIENT2_Y2_HERE':str(int(R2[1,1])),
                    'HERE_ORIENT2_Y3_HERE':str(int(R2[1,2])),
                    'HERE_ORIENT2_Z1_HERE':str(int(R2[2,0])),
                    'HERE_ORIENT2_Z2_HERE':str(int(R2[2,1])),
                    'HERE_ORIENT2_Z3_HERE':str(int(R2[2,2])),
                    'HERE_MISORIENTATION_ANGLE_HERE':str(angle),
                    'HERE_SIGMA_HERE':str(sigmaval),
                    'HERE_GBE_HERE':str(Al_GBE),
                    'HERE_GBSTRUCTURE_BLOB_HERE':str(url)
                    }

    recordfile=datadir+'/lammps15May15.Al.stgb'+'001'+'_'+repr(int(index))+'.xml'
    recordname='lammps15May15.Al.stgb'+'001'+'_'+repr(int(index))

    with open(recordtemplate) as infile, open(recordfile, 'w') as outfile:
        for line in infile:
            for src, target in replacements.iteritems():
                line = line.replace(src, target)
            outfile.write(line)

    print 'uploading record'
    response = mdcs.curate_as(recordfile,recordname,host,user,pswd,template_title=template_name)
    print "Response:",response


g.close()
