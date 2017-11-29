#!/usr/bin/python
import numpy as np
import sys 
from matplotlib import pyplot as plt

#=====================================================
#
#  Code: read_lv.py
#
#  Purpose: Open/Read "lv" binary files produced in
#           from simulations and plot results
#
#  Keywords: file - athena lv binary file to be opened
#                   (e.g. bar.0000.lv.bin)
#   (these are always written to root directory - id0)
#
#  Usage: import read_lv 
#
#          
#  Author: John Dupuy 
#          UNC Chapel Hill
#  Date:   11/28/17
#=====================================================


# Read binary file 
def readbin(fl,precision, **kwargs):
    
    # User specified values for nx1_dom, Nang
    for key in kwargs:
        if key == 'nlong':
            nlong  = kwargs[key]
        if key == 'nvbins':
            nvbins = kwargs[key]
    try: 
        file = open(fl,'rb')
    except:
        print('[read_lv]: failed to supply input file')
        quit()
    
    if precision==32:
        prec=np.float32
    elif precision==64:
        prec=np.float64
    else:
        print('[read_lv]: Invalid precision input')
        quit()
    
    file.seek(0,2)
    eof = file.tell()
    file.seek(0,0)

    # Read integ
    integ = np.fromfile(file,dtype=np.uint32,count=2)
    nlong, nvbins = integ[0], integ[1]

    # Read dat
    dat = np.fromfile(file,dtype=prec,count=7)
    # Parse dat
    t        = dat[0]   # [Myr]
    minvlos  = dat[1]   
    maxvlos  = dat[2]   
    minl     = dat[3]*180./np.pi   
    maxl     = dat[4]*180./np.pi

   
    # Read lv diagram
    lvdiag = np.zeros((nlong,nvbins))
    for i in range(nlong):
        lvdiag[i] = np.fromfile(file,dtype=prec,count=nvbins)

    # construct velocity, longitude array
    lvals = np.linspace(minl   ,maxl   ,nlong )
    vvals = np.linspace(minvlos,maxvlos,nvbins)
    file.close()

    print(lvdiag) 
    
    plt.imshow(lvdiag.T, origin='lower', extent=[minl,maxl,minvlos,maxvlos],aspect='auto')
    #plt.savefig('nompi.eps')
    plt.show()
    return t, lvals, vvals, lvdiag

readbin('id0/bar.0001.lv.bin',32)
