#!/usr/bin/python
import numpy as np
import sys
import matplotlib.pyplot as plt

# This file contains all the functions to read binary files 
# dumped by Athena, including user-enrolled binary dumps
#=====================================================
#
#  Function: read_bin.py
#
#  Purpose: Open/Read athena binary file and export 
#           return content to main code
#
#  Keywords: file - athena binary file to be opened
#
#  Usage: import read_bin 
#         nx,ny,nz,x,y,z,\
#         d,Mx,My,Mz,e,s,\
#         bx,by,bz,phi,\
#         gamm1,cs,t,dt,nscalars   =  read_bin.read_bin('test.bin')
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
def read_bin(fl,precision):

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

#=====================================================
#
#  Function: read_otf
#
#  Purpose: Open/Read "otf" binary files produced in
#           from hvc simulations and plot results
#
#  Keywords: file - athena otf binary file to be opened
#                   (e.g. bar.0000.otf.bin)
#   (these are always written to root directory - id0)
#
#  Usage: import read_bin 
#         read_bin.read_otf(fl, precision)  
#
#          
#  Author: John Dupuy 
#          UNC Chapel Hill
#  Date:   11/7/17
#=====================================================


# Read binary file 
def read_otf(fl,precision, **kwargs):
    
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

#=====================================================
#
#  Function: read_lv
#
#  Purpose: Open/Read "lv" binary files produced in
#           from simulations and plot results
#
#  Keywords: file - athena lv binary file to be opened
#                   (e.g. bar.0000.lv.bin)
#   (these are always written to root directory - id0)
#
#  Usage: import read_bin
#         read_bin.read_lv(fl,precision) 
#
#          
#  Author: John Dupuy 
#          UNC Chapel Hill
#  Date:   11/28/17
#=====================================================


# Read binary file 
def read_lv(fl,precision, **kwargs):
    
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
    integ = np.fromfile(file,dtype=np.uint32,count=3)
    nlong, nvbins, ntracer = integ[0], integ[1], integ[2]

    # Read dat
    dat = np.fromfile(file,dtype=prec,count=5)
    # Parse dat
    t        = dat[0]   # [Myr]
    minvlos  = dat[1]   
    maxvlos  = dat[2]   
    minl     = dat[3]*180./np.pi   
    maxl     = dat[4]*180./np.pi

       
    # Read lv diagram
    lv = np.zeros((ntracer,nlong,nvbins))
    for it in range(ntracer):
        for i in range(nlong):
            lv[it][i]   = np.fromfile(file,dtype=prec,count=nvbins)
    
    # construct velocity, longitude array
    lvals = np.linspace(minl   ,maxl   ,nlong )
    # Handle case where no gas has T < 1e4
    if abs(minvlos) > 1e28:
        minvlos = -300
        maxvlos =  300
        nvbins  = 2

    vvals = np.linspace(minvlos,maxvlos,nvbins)
    
    file.close()

    return t, lvals, vvals, lv[0].T, lv[1].T

#=====================================================
#
#  Function: read_sii
#
#  Purpose: Open/Read "sii" binary files produced in
#           from simulations and plot results
#
#  Keywords: file - athena sii binary file to be opened
#                   (e.g. bar.0000.sii.bin)
#   (these are always written to root directory - id0)
#
#  Usage: import read_bin
#         read_bin.read_sii(fl,precision) 
#
#          
#  Author: John Dupuy 
#          UNC Chapel Hill
#  Date:   11/28/17
#=====================================================
# Read binary file 
def read_sii(fl,precision, **kwargs):
    
    try: 
        file = open(fl,'rb')
    except:
        print('[read_sii]: failed to supply input file')
        quit()
    
    if precision==32:
        prec=np.float32
    elif precision==64:
        prec=np.float64
    else:
        print('[read_sii]: Invalid precision input')
        quit()
    
    file.seek(0,2)
    eof = file.tell()
    file.seek(0,0)

    # Read Nx, Ny
    integ = np.fromfile(file,dtype=np.uint32,count=2)
    Nx = integ[0]
    Ny = integ[1] 

    # Read dat
    dat = np.fromfile(file,dtype=prec,count=5)
    # Parse dat
    t        = dat[0]   # [Myr]
    minx     = dat[1]   # [pc]   
    maxx     = dat[2]   # [pc]
    miny     = dat[3]   # [pc]
    maxy     = dat[4]   # [pc]

    # Read sii grid
    griddata = np.zeros((Ny,Nx))
    for j in range(Ny):
         griddata[j] = np.fromfile(file,dtype=prec,count=Nx)

    # construct x,y arrays
    x = np.linspace(minx, maxx, Nx)
    y = np.linspace(miny, maxy, Ny)

    
    file.close()

    return t, x, y, griddata 

