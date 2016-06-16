#!/usr/bin/env python

import os
import commands
import numpy as np
from math import acos,cos,sin,pi,sqrt

def header_data(fn,xy,xz,yz,n,xl,xh,yl,yh,zl,zh,U):
    with open(fn,'w') as f:
        f.write('#%s structure: DFT minimized energy %g\n'%(fn[:-5],U))
        f.write('%d atoms\n2 atom types\n\n'%n)
        f.write('%g %g xlo xhi\n'%(xl,xh))
        f.write('%g %g ylo yhi\n'%(yl,yh))
        f.write('%g %g zlo zhi\n'%(zl,zh))
        f.write('%g %g %g xy xz yz\n\nAtoms\n\n'%(xy,xz,yz))

def rotate_pos(abc,V,A,B,C,pos):
#    return np.dot(abc,\
#                  np.array([float(pos[1]),float(pos[2]),float(pos[3])]))/V
    return np.dot(abc,\
           np.dot(np.array([np.cross(B,C),np.cross(C,A),np.cross(A,B)]),\
                  np.array([float(pos[1]),float(pos[2]),float(pos[3])])))/V

path='./'
dirs=os.listdir(path)
elem={'C':2, 'B':1}
lx=0.;ly=0.;lz=0.
al=0.;be=0.;ga=0.
xy=0.;xz=0.;yz=0.
U=0.; Ud={}

status, files = commands.getstatusoutput("find ./ -name \"*out*\" | sort ")
files=files.split('\n')
files=[ f for f in files if 'b1' in f and 'ata' not in f ]

for f in files:
    fn=f.replace('/outdump','').replace('./','').replace('-','')
    fn=''.join([ fs if fs == 'p' else ( fs if fs == 'e' else fs.upper() )\
                                                             for fs in fn ])
    L=open(f, 'r').readlines()
    data={}; coords=0; qM=0; qH=0; Ud[fn]=0
    for l in L:
        if 'Vector a' in l:
            lx=float(l.split()[4])
            if 'CELL_REF' in l:
                A=np.array(l.split()[4:7],dtype=float)
        elif 'Vector b' in l:
            ly=float(l.split()[5])
            xy=float(l.split()[4])
            if 'CELL_REF' in l:
                B=np.array(l.split()[4:7],dtype=float)
        elif 'Vector c' in l:
            lz=float(l.split()[6])
            xz=float(l.split()[4])
            yz=float(l.split()[5])
            if 'CELL_REF' in l:
                C=np.array(l.split()[4:7],dtype=float)
        elif 'alpha' in l:
            al=float(l.split()[-1])/180.*pi
        elif 'beta' in l:
            be=float(l.split()[-1])/180.*pi
        elif 'gamma' in l:
            ga=float(l.split()[-1])/180.*pi
        elif 'ATOMIC COORDINATES' in l:
            coords=1
        elif coords:
            ls=l.split()
            try:
                data[int(ls[0])]=\
                    [elem[ls[2]],0,0,float(ls[4]),float(ls[5]),float(ls[6])]
            except:
                if len(data)==120:
                    coords=0
        elif 'Mulliken Pop' in l:
            qM=1
        elif qM:
            ls=l.split()
            try:
                data[int(ls[0])][1]=float(ls[-1])
            except:
                if 'Total charge' in l:
                    qM=0
        elif 'Hirshfeld Ch' in l:
            qH=1
        elif qH:
            ls=l.split()
            try:
                data[int(ls[0])][2]=float(ls[-1])
            except:
                if 'Total Charge' in l:
                    qH=0
        elif 'ENERGY' in l:
            U=float(l.split()[-1])
            Ud[fn]=U

    a=lx;b=sqrt(ly*ly+xy*xy);c=sqrt(lz*lz+xz*xz+yz*yz)
    abc=np.array([[lx,0.,0.],[b*cos(ga),b*sin(ga),0.],[c*cos(be),0.,0.]])
    abc[2][1]=(b*c*cos(al)-abc[1][0]*abc[2][0])/abc[1][1]
    abc[2][2]=sqrt(c*c-abc[2][1]**2.-abc[2][0]**2.)
    abc=np.array([[lx,0,0],[xy,ly,0],[xz,yz,lz]])
    abcT=abc.T
    #C=[c,0.,0.]
    #A=[a*sin(be),0.,a*cos(be)]
    #de=acos((cos(ga)-cos(al)*cos(be))/sin(al)/sin(be))
    #B=[b*sin(al)*cos(de),b*sin(al)*sin(de),b*cos(al)]
    #A=[1.,0.,0.];B=[0.,1.,0.];C=[0.,0.,1.]
    V=a*b*c*sqrt(1-(cos(al))**2.-(cos(ga))**2.-(cos(be))**2.\
                       +2.*cos(al)*cos(be)*cos(ga))
    header_data('%sMull.data'%fn,xy,xz,yz,len(data),0,lx,0,ly,0,lz,U)
    with open('%sMull.data'%fn,'a') as fout:
        for (i, dat) in sorted(data.items()):
            atype=dat[0]
            pos=rotate_pos(abcT,V,A,B,C,dat[-4:])
            q=dat[1]
            fout.write('%d %d %.4f %g %g %g\n'%\
                       (i,atype,q,pos[0],pos[1],pos[2]))
        fout.write('\nMasses\n\n')
        fout.write('1 10.811000\n2 12.010700')

    header_data('%sHirsh.data'%fn,xy,xz,yz,len(data),0,lx,0,ly,0,lz,U)
    with open('%sHirsh.data'%fn,'a') as fout:
        for (i, dat) in sorted(data.items()):
            atype=dat[0]
            pos=rotate_pos(abcT,V,A,B,C,dat[-4:])
            q=dat[2]
            fout.write('%d %d %.4f %g %g %g\n'%\
                       (i,atype,q,pos[0],pos[1],pos[2]))
        fout.write('\nMasses\n\n')
        fout.write('1 10.811000\n2 12.010700')

