import DiffractionSim as DS
from numpy import *
import time
import winsound
import os
Freq = 12000  # Set Frequency To 2500 Hertz
Dur = 600  # Set Duration To 1000 ms == 1 second



A = DS.DiffractionSim()
A.read_data()

A.map2saed_screen(zone=[3,0,2.14], pixel_dim=512, shell_thickness=0.04)

winsound.Beep(Freq, Dur)

A.image_screen(screen=None, exposure_c=1, sigma=2, blur_steps=2,plot=True)
