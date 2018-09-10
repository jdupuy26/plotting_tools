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
#         gamm1,cs,t,dt,nscalars,\
#         clesshd              =  read_bin.read_bin('test.bin')
#
#  WARNING: This is compatible with the dump_binary.c file in the
#           git directory. To use with the default dump_binary file
#           the naddvar variable needs to be removed...
#
#          
#  Author: Chris Frazer
#          UNC Chapel Hill
#  Modified: (9/15/17) returns internal energy
#            (5/14/18) returns collisionless variables, if present 
#
#=====================================================

 
# Read binary file
def read_bin(fl,precision):

    try: 
        file = open(fl,'rb')
    except:
        print('[read_bin]: failed to supply file input')
        quit()


    if precision == 32:
        prec= np.float32
    elif precision == 64:
        prec=np.float64
    else:
        print('[read_bin]: failed to assign appropriate precision')
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
    nvar      = ndata[3] # number of variables 
    nscalars  = ndata[4] # number of scalars
    ngrav     = ndata[5] # self-gravity (1: on, 0: off)
    npart     = ndata[6] # number of particles
    naddvar   = ndata[7] # number of additional vars, including CLESSHD 

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
    ie  = np.zeros(shape)
    s   = np.zeros(scalshape)
    bx  = np.zeros(shape)
    by  = np.zeros(shape)
    bz  = np.zeros(shape)
    phi = np.zeros(shape)
    P   = np.zeros(shape)
    
    # Possible additional variables 
        # allocate arrays for collisionless fluid 
    dcl  = np.zeros(shape)
    M1cl = np.zeros(shape)
    M2cl = np.zeros(shape)
    M3cl = np.zeros(shape)
    E11  = np.zeros(shape)
    E22  = np.zeros(shape)
    E33  = np.zeros(shape) 
    E12  = np.zeros(shape)
    E13  = np.zeros(shape)
    E23  = np.zeros(shape) 
    p11  = np.zeros(shape)
    p22  = np.zeros(shape)
    p33  = np.zeros(shape) 

    #read in arrays
    d  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
    Mx = np.fromfile(file,dtype=prec,count=count).reshape(shape)
    My = np.fromfile(file,dtype=prec,count=count).reshape(shape)
    Mz = np.fromfile(file,dtype=prec,count=count).reshape(shape)

    a = nvar - nscalars
    
    if (a == 5): # Adiabatic & no dual energy
        e  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
    if (a == 6): # Adiabatic & dual energy
        e   = np.fromfile(file,dtype=prec,count=count).reshape(shape)
        ie  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
    if (a == 8): # Adiabatic + MHD & no dual energy 
        e   = np.fromfile(file,dtype=prec,count=count).reshape(shape)
        bx  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
        by  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
        bz  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
    
    if nscalars > 0:
        for i in range(0,nscalars):
            s[:,:,:,i]  = np.fromfile(file,dtype=prec,count=count).reshape(shape)

    if ngrav > 0:
        phi = np.fromfile(file,dtype=prec,count=count).reshape(shape)
   
    if naddvar == 10: # CLESS and no dual energy 
        # Read in collisionless variables 
        dcl  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
        M1cl = np.fromfile(file,dtype=prec,count=count).reshape(shape) 
        M2cl = np.fromfile(file,dtype=prec,count=count).reshape(shape)
        M3cl = np.fromfile(file,dtype=prec,count=count).reshape(shape)
        E11  = np.fromfile(file,dtype=prec,count=count).reshape(shape) 
        E22  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
        E33  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
        E12  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
        E13  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
        E23  = np.fromfile(file,dtype=prec,count=count).reshape(shape)

    elif (naddvar == 13): # CLESS and dual energy 
        # Read in collisionless variables 
        dcl  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
        M1cl = np.fromfile(file,dtype=prec,count=count).reshape(shape) 
        M2cl = np.fromfile(file,dtype=prec,count=count).reshape(shape)
        M3cl = np.fromfile(file,dtype=prec,count=count).reshape(shape)
        E11  = np.fromfile(file,dtype=prec,count=count).reshape(shape) 
        E22  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
        E33  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
        E12  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
        E13  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
        E23  = np.fromfile(file,dtype=prec,count=count).reshape(shape)
        p11  = np.fromfile(file,dtype=prec,count=count).reshape(shape) 
        p22  = np.fromfile(file,dtype=prec,count=count).reshape(shape) 
        p33  = np.fromfile(file,dtype=prec,count=count).reshape(shape) 


    # Package the collisionless data in a tuple 
    clesshd = (dcl, M1cl, M2cl, M3cl, E11, E22, E33, E12, E13, E23, p11, p22, p33) 
                         
    return nx,ny,nz,x,y,z,\
            d,Mx,My,Mz,e,ie,s,\
            bx,by,bz,phi,\
            gamm1,cs,t,dt,nscalars,\
            clesshd, coordsys 

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
    # Read A1, A2 [unitless] 
    A1  = np.fromfile(file,dtype=prec,count=nx1_dom)
    A2  = np.fromfile(file,dtype=prec,count=nx1_dom)
    
    file.close()
    return t, mhvc, rhvc, rpos,\
           acc_rate, facvhvc, ahvc,\
           mcR, mcL,\
           r, ang, vrot,\
           A1, A2

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
    
    ipos = 0
    for key in kwargs:
        if key == 'ipos':
            ipos = kwargs[key] 
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
   
    # if ipos = 0, this will terminate after reading 
    #       the first (l,v) diagram 
    # if ipos = 1, then this loop runs twice and 
    #       the function will return the second (l,v) diagram  
    # And so on for greater ipos numbers 
    for i in range(ipos+1): 
        # Read integ
        integ = np.fromfile(file,dtype=np.uint32,count=3)
        nlong, nvbins, ntracer = integ[0], integ[1], integ[2]

        # Read dat
        dat = np.fromfile(file,dtype=prec,count=8)
        # Parse dat
        t        = dat[0]   # [Myr]
        minvlos  = dat[1]   
        maxvlos  = dat[2]   
        minl     = dat[3]*180./np.pi   
        maxl     = dat[4]*180./np.pi
        Robs     = dat[5]
        pobs     = dat[6]
        vobs     = dat[7]

        # Read lv diagram
        #lv    = np.zeros((ntracer,nlong,nvbins))
        lv = np.fromfile(file,dtype=prec,count=ntracer*nlong*nvbins)
        lv = np.reshape(lv, (ntracer,nlong,nvbins) ).astype(np.float64) # convert to float64  
        #for it in range(ntracer):
        #    for i in range(nlong):
        #        lv[it][i]   = np.fromfile(file,dtype=prec,count=nvbins)
        
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

    old   = False
    # way to still be able to plot old sii files 
    for key in kwargs:
        if key == 'old':
            old    = kwargs[key]
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
    if old:
        integ  = np.fromfile(file,dtype=np.uint32,count=2)
        nlines = 1
        Nx     = integ[0]
        Ny     = integ[1]
    else: 
        integ  = np.fromfile(file,dtype=np.uint32,count=3)
        nlines = integ[0]
        Nx     = integ[1]
        Ny     = integ[2] 

    # Read dat
    dat = np.fromfile(file,dtype=prec,count=5)
    # Parse dat
    t        = dat[0]   # [Myr]
    minx     = dat[1]   # [pc]   
    maxx     = dat[2]   # [pc]
    miny     = dat[3]   # [pc]
    maxy     = dat[4]   # [pc]

    # Read sii grid
    griddata = np.zeros((nlines,Ny,Nx))
    for iline in range(nlines):
        for j in range(Ny):
            griddata[iline,j] = np.fromfile(file,dtype=prec,count=Nx)

    # construct x,y arrays
    x = np.linspace(minx, maxx, Nx)
    y = np.linspace(miny, maxy, Ny)

    
    file.close()

    return t, x, y, griddata 

