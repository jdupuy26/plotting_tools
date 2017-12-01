#!/usr/bin/python
import numpy as np
import sys 
from matplotlib import pyplot as plt

#=====================================================
#
#  Code: read_otf.py
#
#  Purpose: Open/Read "otf" binary files produced in
#           from hvc simulations and plot results
#
#  Keywords: file - athena otf binary file to be opened
#                   (e.g. bar.0000.otf.bin)
#   (these are always written to root directory - id0)
#
#  Usage: import read_otf 
#
#          
#  Author: John Dupuy 
#          UNC Chapel Hill
#  Date:   11/7/17
#=====================================================


# Read binary file 
def readbin(fl,precision, **kwargs):
    
    # Set defaults for nx1_dom, Nang
    nx1_dom = 256
    Nang    = 8
    
    # User specified values for nx1_dom, Nang
    for key in kwargs:
        if key == 'nx1_dom':
            nx1_dom = kwargs[key]
        if key == 'Nang':
            Nang = kwargs[key]
    try: 
        file = open(fl,'rb')
    except:
        print('[read_otf]: failed to supply input file')
        quit()
    
    if precision==32:
        prec=np.float32
    elif precision==64:
        prec=np.float64
    else:
        print('[read_otf]: Invalid precision input')
        quit()
    
    file.seek(0,2)
    eof = file.tell()
    file.seek(0,0)

    # Read dat
    dat = np.fromfile(file,dtype=prec,count=7)
    # Parse dat
    t        = dat[0]   # [Myr]
    mhvc     = dat[1]   # [M_sun]
    rhvc     = dat[2]   # [pc]
    rpos     = dat[3]   # [pc]
    acc_rate = dat[4]   # [M_sun/Myr]
    facvhvc  = dat[5]   # [dimensionless]
    ahvc     = dat[6]   # [radians]
    
    # Read scal [M_sun] 
    scal = np.fromfile(file,dtype=prec,count=2)
    # mcR = mwin_circleR, mcL = mwin_circleL
    mcR, mcL = scal[0], scal[1]
    
    # Read radii [pc]
    r   = np.fromfile(file,dtype=prec,count=nx1_dom)
    
    # Read angles [rad]
    ang = np.fromfile(file,dtype=prec,count=Nang)

    # Read vrot [pc/Myr] 
    vrot = np.zeros((Nang, nx1_dom))
    for i in range(Nang):
        vrot[i] = np.fromfile(file,dtype=prec,count=nx1_dom)

    file.close()
    return t, mhvc, rhvc, rpos,\
           acc_rate, facvhvc, ahvc,\
           mcR, mcL,\
           r, ang, vrot
