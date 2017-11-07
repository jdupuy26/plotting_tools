#!/usr/bin/python
import numpy as np
import sys
from matplotlib import pyplot as plt

#=====================================================
#
#  Code: read_athena_bin.py
#
#  Purpose: Open/Read athena binary file and export 
#           return content to main code
#
#  Keywords: file - athena binary file to be opened
#
#  Usage: import read_athena_bin 
#         nx,ny,nz,x,y,z,\
#         d,Mx,My,Mz,e,s,\
#         bx,by,bz,phi,\
#         gamm1,cs,t,dt,nscalars   =  read_athena_bin.readbin('test.bin')
#
#  WARNING: This is compatible with the dump_binary.c file in the
#           git directory. To use with the default dump_binary file
#           the naddvar variable needs to be removed...
#
#          
#  Author: Chris Frazer
#          UNC Chapel Hill
#  Modified: (9/15/17) returns internal energy
#
#=====================================================

 
# Read binary file
def readbin(fl,precision):

    try: 
        file = open(fl,'rb')
    except:
        print 'read_athena_bin: failed to supply file input'
        quit()


    if precision == 32:
        prec= np.float32
    elif precision == 64:
        prec=np.float64
    else:
        print 'read_athena_bin: failed to assign appropriate precision'
        quit()


    #check location of EOF
    file.seek(0,2)
    eof = file.tell()
    file.seek(0,0)

    #Read variables
    coordsys  = np.fromfile(file,dtype=np.int32,count=1)[0]
    ndata     = np.fromfile(file,dtype=np.int32,count=8)         #ndata[7] is number of additional vars
    gamm1,cs  = np.fromfile(file,dtype=prec,count=2)
    t,dt      = np.fromfile(file,dtype=prec,count=2)

    #break up ndata
    nx        = ndata[0]
    ny        = ndata[1]
    nz        = ndata[2]
    nvar      = ndata[3]
    nscalars  = ndata[4]
    ngrav     = ndata[5]
    npart     = ndata[6]
    naddvar   = ndata[7]

    #x-y-z arrays
    x = np.fromfile(file,dtype=prec,count=nx)
    y = np.fromfile(file,dtype=prec,count=ny)
    z = np.fromfile(file,dtype=prec,count=nz)

    #determine shape of arrays
    shape = (nz,ny,nx)
    scalshape = (nz,ny,nx,nscalars)
    count = np.prod(shape)    

    #allocate all arrays
    d   = np.zeros(shape)
    Mx  = np.zeros(shape)
    My  = np.zeros(shape)
    Mz  = np.zeros(shape)
    e   = np.zeros(shape)
    s   = np.zeros(scalshape)
    bx  = np.zeros(shape)
    by  = np.zeros(shape)
    bz  = np.zeros(shape)
    phi = np.zeros(shape)
    ie  = np.zeros(shape)
    P   = np.zeros(shape)
    

    #read in arrays
    d  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
    Mx = np.fromfile(file,dtype=prec,count=count).reshape(shape)
    My = np.fromfile(file,dtype=prec,count=count).reshape(shape)
    Mz = np.fromfile(file,dtype=prec,count=count).reshape(shape)

    a = nvar-nscalars-naddvar
    if a == 5 or a == 8:   #not Barotropic
        e  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
    if a == 7 or a == 8:   #includes Magnetic fields
        bx  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
        by  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
        bz  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
    if nscalars > 0:
        for i in range(0,nscalars):
            s[:,:,:,i]  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
    if naddvar > 0:
        ie = np.fromfile(file,dtype=prec,count=count).reshape(shape)
    
    if ngrav > 0:
        phi = np.fromfile(file,dtype=prec,count=count).reshape(shape)
                         
    if gamm1 != 0:
        if naddvar > 0:
            P = gamm1*ie
        else:
            P = gamm1*(e-0.5*(Mx**2 + My**2 + Mz**2))/d 
            if a == 7: 
                P -= gamm1*0.5*(bx**2 + by**2 + bz**2)
    elif gamm1 == 0:
        P = cs*cs*d

    return nx,ny,nz,x,y,z,\
            d,Mx,My,Mz,e,ie,s,\
            bx,by,bz,phi,\
            gamm1,cs,t,dt,nscalars


                                                            
