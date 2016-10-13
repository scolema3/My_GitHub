import glob, os                           # Python standard library
import sys                                # Python standard library
import time                               # Python standard library
import unicodedata                        # Python standard library
import numpy as np                        # NEED to be installed
import json                               # Python standard library
import requests                           # NEED to be installed
import xml.dom.minidom                    # Python standard library
from xml.dom.minidom import parseString   # Python standard library
import collections                        # Python standard library
import xmltodict                          # attached
import dicttoxml                          # attached
import xlrd                               # NEED to be installed
import matplotlib.pyplot as plt
import matplotlib.colorbar as clb
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec
from matplotlib import ticker
from matplotlib.ticker import LogLocator
import mdcs

'''
This is just an very quick example of using the mdcs tools to extract information 
from the grain boundary database (GBD).



'''
#################################
USER='scolema3'
PSWD='MDCS1pswd'
MDCS_URL='https://mdcs1.nist.gov'
template_name='record-calculation-grain-boundary'

#################################

       
template_id=mdcs.templates.current_id(MDCS_URL,USER,PSWD,title=template_name)        
if template_id==None: 
    print "Error: template_name not found at MDCS_URL\n"
    print("       template_name: "+template_name)
    print("       MDCS_URL: "+MDCS_URL+"\n")
    quit()

print "Looking for records following this Schema"    
template_results=mdcs.explore.select(MDCS_URL,USER,PSWD,template=template_id)
n_records=len(template_results)


print("Found "+str(n_records)+" record entries for "+template_name)

record_conent=[]
x=[]
y=[]

for i in range(0,n_records): 
    record=(xmltodict.parse(xml.dom.minidom
                            .parseString(template_results[i]['content']).toprettyxml()))
    # Potential Info
    pot_name=record['calculation-grain-boundary']['potential']['id']
    pot_key=record['calculation-grain-boundary']['potential']['key']
    pot_elements=record['calculation-grain-boundary']['potential']['element']
    
    # Code Info
    comments_method=record['calculation-grain-boundary']['atomistic-method']['comment']
    code=record['calculation-grain-boundary']['atomistic-method']['code']['name']
    
    # GB Info
    gb_name=record['calculation-grain-boundary']['bicrystal']['id']
    lx=record['calculation-grain-boundary']['bicrystal']['dimension']['x']['value']
    ly=record['calculation-grain-boundary']['bicrystal']['dimension']['y']['value']
    lz=record['calculation-grain-boundary']['bicrystal']['dimension']['z']['value']
    lx_u=record['calculation-grain-boundary']['bicrystal']['dimension']['x']['unit']
    ly_u=record['calculation-grain-boundary']['bicrystal']['dimension']['y']['unit']
    lz_u=record['calculation-grain-boundary']['bicrystal']['dimension']['z']['unit']
    natoms=record['calculation-grain-boundary']['bicrystal']['number-of-atoms']
    
    try:
        # Two GB information listed (PBC)
        multiple=record['calculation-grain-boundary']['bicrystal'][0]['grain-boundary']['area']['value']
        gb_normal=record['calculation-grain-boundary']['bicrystal'][0]['grain-boundary']['normal']['value']
        gb_normal_u=record['calculation-grain-boundary']['bicrystal'][0]['grain-boundary']['normal']['unit']
        gb_area=record['calculation-grain-boundary']['bicrystal'][0]['grain-boundary']['area']['value']
        gb_area_u=record['calculation-grain-boundary']['bicrystal'][0]['grain-boundary']['area']['unit']
    except:
        # Single GB information listed
        gb_normal=record['calculation-grain-boundary']['bicrystal']['grain-boundary']['normal']['value']
        gb_normal_u=record['calculation-grain-boundary']['bicrystal']['grain-boundary']['normal']['unit']
        gb_area=record['calculation-grain-boundary']['bicrystal']['grain-boundary']['area']['value']
        gb_area_u=record['calculation-grain-boundary']['bicrystal']['grain-boundary']['area']['unit']

    # Grain 1
    ID=0
    grain1_prototype=record['calculation-grain-boundary']['bicrystal']['grain'][ID]['crystal-prototype']['name']
    # Row major rotation matrix: [x1 x2 x3; y1 y2 y3; z1 z2 z3]
    grain1_rot=np.matrix([[float(record['calculation-grain-boundary']['bicrystal']['grain'][ID]['orientation']['x']['h']), float(record['calculation-grain-boundary']['bicrystal']['grain'][ID]['orientation']['x']['k']), float(record['calculation-grain-boundary']['bicrystal']['grain'][ID]['orientation']['x']['l'])], [float(record['calculation-grain-boundary']['bicrystal']['grain'][ID]['orientation']['y']['h']), float(record['calculation-grain-boundary']['bicrystal']['grain'][ID]['orientation']['y']['k']), float(record['calculation-grain-boundary']['bicrystal']['grain'][ID]['orientation']['y']['l'])], [float(record['calculation-grain-boundary']['bicrystal']['grain'][ID]['orientation']['z']['h']), float(record['calculation-grain-boundary']['bicrystal']['grain'][ID]['orientation']['z']['k']), float(record['calculation-grain-boundary']['bicrystal']['grain'][ID]['orientation']['z']['l'])]])
    
    # Grain 2
    ID=1
    grain2_prototype=record['calculation-grain-boundary']['bicrystal']['grain'][ID]['crystal-prototype']['name']
    # Row major rotation matrix: [x1 x2 x3; y1 y2 y3; z1 z2 z3]
    sigma=record['calculation-grain-boundary']['bicrystal']['sigma-value']
    misorientation=record['calculation-grain-boundary']['bicrystal']['misorientation']['value']
    misorientation_u=record['calculation-grain-boundary']['bicrystal']['misorientation']['unit']
    
    # GBE
    gbe=record['calculation-grain-boundary']['grain-boundary-energy']['value']
    gbe_u=record['calculation-grain-boundary']['grain-boundary-energy']['unit']
    
    # Temperature
    T=record['calculation-grain-boundary']['temperature']['value']
    T_u=record['calculation-grain-boundary']['temperature']['unit']

    # Save data to plot into a list
    x.append(misorientation)
    y.append(gbe)
    
fig, ax = plt.subplots(figsize=(12, 9))
ax.plot(x,y)
plt.show()
  


