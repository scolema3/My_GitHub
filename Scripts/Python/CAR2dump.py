#!/usr/bin/python
import sys

def header(fn,t=0,n=1,xmn=0,xmx=1,ymn=0,ymx=1,zmn=0,zmx=1): 
    with open(fn, 'w') as g:
        g.write('ITEM: TIMESTEP\n%s\n'%t)
        g.write('ITEM: NUMBER OF ATOMS\n%s\n'%n)
        g.write('ITEM: BOX BOUNDS pp pp pp\n')
        g.write('%s %s\n%s %s\n%s %s\n'%(xmn,xmx,ymn,ymx,zmn,zmx))
        g.write('ITEM: ATOMS id type x y z\n')


fn = str(sys.argv).split()[1][1:-2]

f = open(fn,'r')
L = f.readlines()
f.close()

x=[]; y=[]; z=[]; ID=[]

for l in L[5:]:
    if l.strip() == 'end':
        break
    ls = l.split()
    x.append(float(ls[1]))
    y.append(float(ls[2]))
    z.append(float(ls[3]))
    if ls[7] == 'B':
        ID.append(1)
    else:
        ID.append(2)

fn = '%sdump'%fn[:-3]
header(fn,0,len(x),min(x),max(x),min(y),max(y),min(z),max(z))
f = open('%s'%fn,'a')

for i in range(len(x)):
    f.write('%d %d %g %g %g\n'%(i+1,ID[i],x[i],y[i],z[i]))
f.close()
