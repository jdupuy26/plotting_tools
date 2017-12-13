import os
import sys
import time
import gobject
import argparse          
import numpy as np
# matplotlib imports
import matplotlib as mpl
import matplotlib.animation as animation
from matplotlib import pyplot as plt
from matplotlib import pylab as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.projections import PolarAxes
import mpl_toolkits.axisartist.floating_axes as floating_axes
from mpl_toolkits.axisartist.grid_finder import DictFormatter
#scipy imports 
from scipy.interpolate import RectBivariateSpline, interp1d, griddata, interp2d
from scipy.ndimage import map_coordinates
from scipy.stats import binned_statistic
from scipy.ndimage.filters import gaussian_filter, uniform_filter

# Import from correct directory
import socket as s
comp = s.gethostname()
if comp == 'thetis': 
    sys.path.insert(0,'/afs/cas.unc.edu/users/j/d/'
                      'jdupuy26/Johns_work/'
                      'misc_scripts/plotting_tools')
elif comp == 'debpad':
    sys.path.insert(0,'/home/jdupuy26/Johns_work/'
                      'Grad_Work/Research/codes/'
                      'athena/misc_scripts/'
                      'plotting_tools')
else: 
    print('[init]: Computer %s not recognized!' % comp) 

from read_bin import read_bin as read_athena_bin
import read_athinput
# import units class 
from units_class import * 
#=========================================================
"""
 Code: plot_bgsbu.py (previously smrplot3D.py)

 Purpose: Automated smr/mpi compatible plotter for Athena binary files
          intended for 3D simulations


 Keywords: Required:
           quant   [string]   - quantity to be visualized
           start   [int]      - starting frame
           fin     [int]      - final frame
           
           Optional:
           step    [int]      - step between dumps (default 1)
           lev     [int]      - smr level to be plotted
           ival    [int]      - cell value in x direction for vertical line
           jval    [int]      - cell value in y direction for vertical line
           qminmax [2x float] - 2 values to specify range of color bar
           proj    [flag]     - (only 2D, for cyl sims) create (l,v) projection plots for bgsbu simulations 
           proj_cart [flag]   - (only 2D, for cart sims) create (l,v) projection plots for bgsbu simulations
           cyl     [flag]     - (only 2D) plotting for cylindrical simulations
           p2c     [flag]     - (only 2D) projects cylindrical data onto cartesian grid for cartesian plotting
           iner    [flag]     - (only 2D) rotate inertial frame sims so that they appear stationary 
           log     [flag]     - log plotted quantity
           slice   [flag]     - take slice through middle or specified ival/jval
           dumpbin [flag]     - dump combined files to single binary
           units   [int]      - convert to normal units (options 1: CGS, 2: SI)
           slide   [flag]     - slide colorbar 
           gif     [flag]     - save animation as gif
		   mp4	   [flag]     - save animation as mp4 movie (requires ffmpeg)
           dval    [int]      - (3D only) cell value of slice
           view    [int]      - (3D only) 0:xy   1: xz  2: yz

 Requires:  read_athena_bin.py
            read_athinput.py

 Calling Sequence:  python plot_bgsbu.py 0 10 --quant 'd'

 Author: Chris Frazer
         UNC Chapel Hill
 Edited for bgsbu: 
				  -- added an option to make mp4 movie (5/30/17)
				  -- made code units be physical 
                  -- added options for cylindrical simulations (--cyl)
                  -- added option to make inertial frame simulations apepar stationary so that 
                     we can compare directly to rotating frame simulations. (--iner) 
                  -- added option to project cylindrical simulations onto a cartesian grid (--p2c) 
                  -- added option to create projection plots in the (l,v) plane (--proj and QUANT must be 'proj') 
 Current work: 
                  -- apply correct units conversions when --units flag is passed (9/20/17)
                  -- Implement [SII] emission maps 
"""
#=========================================================

#User Specified Priors
athdir = os.getcwd()+'/'
# this bit of python syntax will get athinput automatically
athin  = [fname for fname in os.listdir(athdir) if fname.startswith("athinput")][0]
xwinsz = 512
ywinsz = 512
ctable = 'magma'
precision = 32

# Location of files for [SII] emission
file_loc = '/afs/cas.unc.edu/users/j/d/jdupuy26/Johns_work/cooling/'
#Read in athinput 
base,params = read_athinput.readath(athdir+athin)

# mean molecular weight 
mu = params[0].mu
# Scale height 
scaleh = params[0].scaleh # pc
# set omg
omg = params[0].omg
# set gammma
gam = params[0].gam



# Line of interest [SII] emission
# 1: 6717  Ang
# 2: 4070  Ang
# 3: 10300 Ang
myline = 1

# Smoothing factor for Intensity 
# HERE
desired_smooth = 1000.  # [pc]

# for the projections
p0 = -20.*np.pi/180.   # [deg] viewing angle
r0 = 8000.  # [pc] Sun's orbital radius
V0 = 220.   # [km/s] Sun's azimuthal speed 



#============FUNCTIONS===========
def get_num(tempnum):
    #get currentnum
    if tempnum < 10  : 
        currentnum = '000'+str(tempnum)
    elif tempnum < 100 : 
        currentnum = '00'+str(tempnum)
    elif tempnum < 1000: 
        currentnum = '0'+str(tempnum)
    elif tempnum >= 1000:             
        currentnum = str(tempnum) 

    return currentnum
                                    


def get_file(dirs,base,np,lev,tempnum,nprocs):
    
    currentnum=get_num(tempnum)
    
    #generate file name
    if (nprocs == 0 ):
        idtag=''
        iddir=''
        levdir=''
        levtag=''
    elif (np==0 and lev==0) :
        idtag=''
        iddir='id'+str(np)+'/'
        levdir=''
        levtag=''
    elif (lev ==0):
        idtag='-id'+str(np)
        iddir='id'+str(np)+'/'
        levdir=''
        levtag=''
    else:
        idtag='-id'+str(np)
        iddir='id'+str(np)+'/'
        levdir='lev'+str(lev)+'/'
        levtag='-lev'+str(lev)
        

    fl = dirs+iddir+levdir+base+idtag+levtag+'.'+currentnum+'.bin'
    return fl


def get_xyzt(file,units):
    #read in binary file
    nx,ny,nz,x,y,z,\
    d,Mx,My,Mz,e,ie,s,\
    bx,by,bz,phi, \
    gamm1,cs,t,dt,nscalars = read_athena_bin(file,precision)
 
    #if   units == 1:
    #    u = units_CGS()
    #elif units == 2: 
    #    u = units_SI()

    #if units != 0:
    #    x *= u.pc  # cm or m
    #    y *= u.pc  # cm or m
    #    z *= u.pc  # cm or m 
    #    t *= u.myr # s  
    return x,y,z,t

def get_quant(file,quant,units):
    #read in binary file
    nx,ny,nz,x,y,z,\
    d,Mx,My,Mz,e,ie,s,\
    bx,by,bz,phi, \
    gamm1,cs,t,dt,nscalars = read_athena_bin(file,precision)
    
    # Set units
    
    if   units == 1:
        u = units_CGS()
    elif units == 2:
        u = units_SI()
    elif units == 0:
        u = units_COMP()

    # Do the unit conversion
    e  *= u.esdens
    ie *= u.esdens
    d  *= u.rhos
    Mx *= u.momsdens
    My *= u.momsdens
    Mz *= u.momsdens
    cs *= u.v

    if quant == 'E':
        arr = e           # surface nrg dens
    elif quant == 'ie':
        arr = ie
    elif quant == 'd':
        arr = d
    elif quant == 'n':
        arr = d/(mu*u.m_h)  
    elif quant == 'p':
        arr = gamm1*ie
    elif quant == 'T':
        d  /= mu*u.m_h     
        arr = gamm1*ie/(u.k_b * d)  # in [K]
    elif quant == 'M':
        arr = np.sqrt(Mx**2 + My**2 + Mz**2)
    elif quant == 'vx':
        arr = Mx/d    
    elif quant == 'vy':
        arr = My/d  
    elif quant == 'vz':
        arr = Mz/d  
    elif quant == 'Mx':
        arr = Mx  
    elif quant == 'My':
        arr = My  
    elif quant == 'Mz':
        arr = Mz  
    elif quant == 'V':
        arr = np.sqrt(Mx**2 + My**2 + Mz**2)/d  
    elif quant =='proj' or quant == 'SII':
        n = d/(mu*u.m_h)   # surface number density
        n /= scaleh*u.pc   # volume number density 
        T = gamm1*ie/(u.k_b * (n*scaleh*u.pc))
        arr = n, Mx/d, My/d + omg*x*u.v, T  # note velocities are in lab frame 
    elif quant[0] == 's':
        num = int(quant[1])
        arr = s[:,:,:,num]
    elif quant == 'x':
        arr = s[:,:,:,0]/d
    elif quant == 'iso_mach':   
        arr = np.sqrt(Mx**2 + My**2 + Mz**2)/(d*cs)
    elif quant == 'adia_mach':  
        vel = np.sqrt(Mx**2 + My**2 + Mz**2)/d
        p   = gamm1*ie
        csound = np.sqrt((gamm1+1)*p/d)
        arr = vel/csound
    elif quant == 'div_v':  
        arr = Mx/d, My/d + omg*x*u.v  # put back into lab frame
        
    return arr

def get_title(quant,units,log,column):

    if units ==0:
        if quant == 'E':
            title = 'Energy Density [code units]'
        elif quant =='ie':
            title = 'Internal Energy Density [code units]'
        elif quant == 'd':
            title = 'Density [code units]'
        elif quant == 'n':
            title = 'Number Density [code units]'
        elif quant == 'p':
            title = 'Pressure [code units]'
        elif quant == 'T':
            title = 'Temperature [code units]'
        elif quant == 'M':
            title = 'Momenutum Density [code units]'
        elif quant == 'vx':
            title = 'X Velocity [code units]'
        elif quant == 'vy':
            title = 'Y Velocity [code units]'
        elif quant == 'vz':
            title = 'Z Velocity [code units]'
        elif quant == 'Mx':
            title = 'X Momentum Density [code units]'
        elif quant == 'My':
            title = 'Y Momentum Density [code units]'
        elif quant == 'Mz':
            title = 'Z Momentum Desnity [code units]'
        elif quant == 'V':
            title = 'Velocity [code units]'
        elif quant[0] == 's':
            num = quant[1]
            title = 'Scalar '+num
        elif quant[0] == 'x':
            title = 'x'
        elif quant == 'proj':
            title = 'Intensity'  
        elif quant == 'iso_mach':
            title = 'Mach #'
        elif quant == 'adia_mach':
            title = 'Mach #'
        elif quant == 'div_v': 
            title = 'div $V$ [km/s]'
        elif quant == 'SII':
            title = '[SII] Emission'

    elif units == 1:
        if quant == 'E':
            title = 'Energy Density [erg cm$^{-2}$]'
        elif quant == 'ie':
            title = 'Internal Energy Density [erg cm$^{-2}$]' 
        elif quant == 'd':
            title = 'Density [g  cm$^{-2}$]'
        elif quant == 'n':
            title = 'Number Density [cm$^{-2}$]'
        elif quant == 'p':
            title = 'Pressure [erg cm$^{-2}$]'
        elif quant == 'T':
            title = 'Temperature [K]'
        elif quant == 'M':
            title = 'Momenutum Density [g cm$^{-1}$ s$^{-1}$]'
        elif quant == 'vx':
            title = 'X Velocity [cm s$^{-1}$]'
        elif quant == 'vy':
            title = 'Y Velocity [cm s$^{-1}$]'
        elif quant == 'vz':
            title = 'Z Velocity [cm s$^{-1}$]'
        elif quant == 'Mx':
            title = 'X Momentum Density [g cm$^{-1}$ s^${-1}$]'
        elif quant == 'My':
            title = 'Y Momentum Density [g cm$^{-1}$ s^${-1}$]'
        elif quant == 'Mz':
            title = 'Z Momentum Desnity [g cm$^{-1}$ s^${-1}$]'
        elif quant == 'V':
            title = 'Velocity [cm s$^{-1}$]'
        elif quant == 'proj':
            title = 'Intensity' 
        elif quant[0] == 's':
            num = quant[1]
            title = 'Scalar '+num
        elif quant[0] == 'x':
            title = 'x'
        elif quant == 'div_v': 
            title = 'Intensity'
        elif quant == 'SII':
            # HERE
            #title = '[SII] Emission' 
            title = 'log(Intensity [arb. units])'
            #title = 'Velocity [km s$^{-1}$]'
            #title = 'div V [km s$^{-1}$]'
            #title = ' log(T$_{shk}$ / T$_{sim}$)'
            

    elif units == 2:
        if quant == 'E':
            title = 'Energy Density [J m$^{-2}$]'
        elif quant == 'ie':
            title = 'Internal Energy Density [J m$^{-2}$]' 
        elif quant == 'd':
            title = 'Density [kg  m$^{-2}$]'
        elif quant == 'n':
            title = 'Number Density [m$^{-2}$]'
        elif quant == 'p':
            title = 'Pressure [J m$^{-2}$]'
        elif quant == 'T':
            title = 'Temperature [K]'
        elif quant == 'M':
            title = 'Momenutum Density [kg m$^{-1}$ s$^{-1}$]'
        elif quant == 'vx':
            title = 'X Velocity [m s$^{-1}$]'
        elif quant == 'vy':
            title = 'Y Velocity [m s$^{-1}$]'
        elif quant == 'vz':
            title = 'Z Velocity [m s$^{-1}$]'
        elif quant == 'Mx':
            title = 'X Momentum Density [kg m$^{-1}$ s^${-1}$]'
        elif quant == 'My':
            title = 'Y Momentum Density [kg m$^{-1}$ s^${-1}$]'
        elif quant == 'Mz':
            title = 'Z Momentum Desnity [kg m$^{-1}$ s^${-1}$]'
        elif quant == 'V':
            title = 'Velocity [m s$^{-1}$]'
        elif quant[0] == 's':
            num = quant[1]
            title = 'Scalar '+num
        elif quant[0] == 'x':
            title = 'x'
        elif quant == 'div_v': 
            title = 'div $V$ [km/s]'
        elif quant == 'SII':
            title = '[SII] Emission'


    if column:
        title = 'column '+title

    if log:
        title = 'log(' + title+')'  
                          
    return title  

def get_SIIemis(line):
    ''' Purpose: Get the emissivity j_v from computed values 
        in Mathematica for [SII] 
        Inputs: line (1,2,3): integer value that 
                corresponds to \lambda = 6717, 4070, 10300 Ang 
                respectively
            Look @ [SII] in Appendix E of Draine. 
            Everything is logged, but units are CGS 
        Output: Interpolated grid in (n, T) space of emissivity
        Note: This gives j_v, to get I_v, must multiply by scaleh
    '''
    
    # these are logged, but in [K], [1/cm^3]
    logTarr, lognarr = np.loadtxt(file_loc+"n_Tarr.dat")

    # logged j_v data
    if line == 1:
        data = np.loadtxt(file_loc+"SII_l6717.dat")
    elif line==2:
        data = np.loadtxt(file_loc+"SII_l4070.dat")
    else: 
        data = np.loadtxt(file_loc+"SII_l10300.dat")
    
    # Interpolate on the data
    data = RectBivariateSpline(lognarr, logTarr, data)
    
    return data


def cyl_axes(fig, rect, mn1, mx1, mn2, mx2):
    ''' Purpose: Create wedge plots 
        fig: figure object 
        rect: 111
        mn1/mx1: min/max r value
        mn2/mx2: min/max theta value 
        
        Adapted from https://stackoverflow.com/questions/
        44189043/plotting-sector-of-polar-plot-so-a-wedge-in-matplotlib
    '''

    tr = PolarAxes.PolarTransform()
    # Set tick labels 
    r_ticks = np.linspace(mn1,mx1,1)
    r_ticks_str = [r'%1.1f' % x for x in r_ticks]
    th_ticks = np.linspace(mn2,mx2,1)
    th_ticks_str = [r'%1.1f' % x for x in th_ticks]
    # Create a dictionary
    tick_formatter1 = DictFormatter(dict(zip(th_ticks,th_ticks_str)))
    tick_formatter2 = DictFormatter(dict(zip(r_ticks,r_ticks_str)))
             
    grid_helper = floating_axes.GridHelperCurveLinear(
        tr, extremes=(mn2, mx2, mn1, mx1),
        tick_formatter2 = tick_formatter2,
        tick_formatter1 = tick_formatter1)

    ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
    fig.add_subplot(ax1)

    # adjust axis
    ax1.axis["left"].set_axis_direction("bottom")
    ax1.axis["right"].set_axis_direction("top")
    ax1.axis["bottom"].toggle(ticklabels=False, label=False)
    ax1.axis["top"].set_axis_direction("bottom")
    ax1.axis["top"].toggle(ticklabels=True, label=True)
    ax1.axis["top"].major_ticklabels.set_axis_direction("top")
    ax1.axis["top"].major_ticklabels.set_pad(10)
    ax1.axis["top"].label.set_axis_direction("top")
    ax1.axis["left"].label.set_text(r"r")
    ax1.axis["top"].label.set_text(r"$\theta$ [rad]")

    # create a parasite axes whose transData in theta, r
    aux_ax = ax1.get_aux_axes(tr)

    aux_ax.patch = ax1.patch  
    ax1.patch.zorder = 0.9  

    return ax1, aux_ax


''' 
Old version
def polar2cartesian(r, t, grid, x, y, order=3):
     Maps polar data (including log grid) to a cartesian grid 
        -------- Inputs ----------  
        r:     radial grid   
        t:     theta grid 
        grid:  r, t data         
        x:     x grid
        y:     y grid 
        order: order of interpolation

        Adapted from: https://stackoverflow.com/questions/2164570/reprojecting-polar-to-cartesian-grid 
    X, Y = np.meshgrid(x,y,indexing='xy') 

    new_r = np.sqrt(X*X+Y*Y)
    new_t = np.arctan2(Y,X)
    # correct negative angles 
    new_t[new_t < 0] = 2*np.pi + new_t[new_t < 0]
    
    ir    = interp1d(r, np.arange(len(r)), bounds_error=False, fill_value=0)
    it    = interp1d(t, np.arange(len(t)), bounds_error=False, fill_value=0)

    new_ir = ir(new_r.ravel())
    new_it = it(new_t.ravel())

    new_ir[new_r.ravel() > r.max()] = len(r)
    new_ir[new_r.ravel() < r.min()] = 0

    #new_it[new_t.ravel() > t.max()] = len(t)
    #new_it[new_t.ravel() < t.min()] = 0

    return map_coordinates(grid, np.array([new_it,new_ir]),order=order).reshape(new_r.shape)
'''



# New version -- needs to be compatible with intensity (l,v) diagrams
def polar2cartesian(r, t, grid, x, y, method='linear'):
    """ Maps polar data (including log grid) to a cartesian grid 
        -------- Inputs ----------  
        r:     radial grid   
        t:     theta grid 
        grid:  r, t data         
        x:     x grid, desired for interpolation
        y:     y grid, desired for interpolation 
        method: method of interpolation
         """
    
    # Make meshgrid of points where we want to know what grid looks like
    X, Y = np.meshgrid(x,y,indexing='xy') 
    
    # Create meshgrid of radial and theta points  
    R, T = np.meshgrid(r,t,indexing='xy')
    # Convert to cartesian meshgrid
    Xold, Yold = R*np.cos(T), R*np.sin(T)

    # Ravel so that the data can be input into griddata 
    Xold = Xold.ravel()
    Yold = Yold.ravel()
    grid = grid.ravel()

    return griddata((Yold,Xold), grid, (Y,X), method=method)


def vel_p2c(vr,vt,x,y):
    """ Converts polar velocity data ON A CARTESIAN GRID 
        to cartesian velocity data on that same CARTESIAN GRID
        ------ Inputs -----------
        vr: velocity in r-dir
        vt: velocity in t-dir
        x,y:  x,y-coordinates of cartesian grid
        ------ Output -----------
        vx, vy: velocities in (x,y) dir
        """ 
    X, Y = np.meshgrid(x,y,indexing='xy')
    R = np.sqrt(X**2 + Y**2)

    vx = vr*X/R - vt*Y/R
    vy = vt*X/R + vr*Y/R

    return vx, vy


def iner_rot(grid, tvec, T):
    """ Rotate inertial data so that it is stationary when plotting
    -------- Inputs ----------                    
    grid: cylindrical grid data
    tvec: vector of time points
    T:    rotation period
    """
    dtvec   = tvec[1]-tvec[0]
    rotated = np.zeros(grid.shape)
    factor1 = simtime/T - (simtime/T)%1.0
    # deterimine interpolation points 
    z = int((simtime-factor1*T)/dtvec)
    # Normalize z 
    if z == len(tvec)-1:
        z = len(tvec)-2

    # Get interpolation points 
    t1 = factor1*T + tvec[z  ]
    t2 = factor1*T + tvec[z+1]
    for j in range(0,int(nx2)):
        jlow = j-z
        jup  = jlow-1
        # normalize
        if (jlow < 0):
            jlow += int(nx2)
        if (jup < 0):
            jup  += int(nx2)
        # interpolate img for times in between t1 and t2:/
        rotated[j,:] = ((t2-simtime)*grid[jlow,:] + 
                        (simtime-t1)*grid[jup ,:])/(t2-t1)
    return rotated


def rad_trans_iso(d,v,vscan,dx):
    ''' Does 1d isothermal radiative transfer for HI 21 cm emission
        along a line of sight (LoS)
        ------- Inputs ----------
        d: density along LoS
        v: v along LoS (i.e. towards or away from the sun)
        vscan: velocity "bins" (scanning frequency of telescope)
        dx: array distance b/w sampling points along LoS
        T: Temperature array along LoS 
        ------- Output -----------
        I: intensity along LoS '''
    
    shape   = np.ones(len(v))
    
    sigma = 10.  # [km/s] 

    invnorm = 1.0/np.sqrt(2.0*np.pi)
    facHI   = 2.2008102e-19             # front factor for optically thin HI line in K (see Draine Ch. 8) 
    # convert d to kg/km^2, dx to km 
    d  *= 1.98855e30/(3.086e13**2)
    dx *= 3.086e13
    # This is a compact way to do what is done in fheitsch's hvc.c radtrans func 
    #I = np.array([ sum(d * dx * np.exp(-(0.5*(vs*shape - v)/sigma)**2)) for vs in vscan ]) 
    I = []
    for vs in vscan:
        doppler = 1.0 - vs/2.9979248e5
        phi = invnorm*np.exp(-(0.5*(vs*shape - v)/sigma)**2)/sigma
        val = facHI * sum(phi * d * dx)/doppler**2
        I.append(val) 
        
    return np.array(I) 


def rad_trans(d,v,T,vscan,dx):
    ''' Does 1d radiative transfer for HI 21 cm emission
        along a line of sight (LoS)
        ------- Inputs ----------
        d: density along LoS  [ g/cm^2 ]
        v: v along LoS (i.e. towards or away from the sun) [km/s]
        T: temperature along LoS  [K]
        vscan: velocity "bins" (scanning frequency of telescope) [km/s]
        dx: array distance b/w sampling points along LoS [pc]
        ------- Output -----------
        I: intensity along LoS 
        
        THIS ASSUMES CGS UNITS AS INPUT
        '''
    u = units_CGS()
    clight = u.c*1e-5                 # speed of light in [km/s]
    invnorm = 1.0/np.sqrt(2.0*np.pi)
    facHI   = 2.2008102e-19             # front factor for optically thin HI line in K (see Draine Ch. 8) 
    facHI   *= 1.e10  # 1/cm^2 to 1/km^2
    
    # This is a compact way to do what is done in fheitsch's hvc.c radtrans func 
    #I = np.array([ sum(d * dx * np.exp(-(0.5*(vs*shape - v)/sigma)**2)) for vs in vscan ]) 
    I = np.zeros(len(vscan))

    # Only include T <= 1e4 K
    ind1 = np.where(T <= 1e4) 
    d  =  d[ind1]
    v  =  v[ind1]
    T  =  T[ind1]
    dx = dx[ind1]
    # Only include T > 1 K
    ind2 = np.where(T > 1)
    d  =  d[ind2]
    v  =  v[ind2]
    T  =  T[ind2]
    dx = dx[ind2] 
    
    shape   = np.ones(len(v))
    #d /= u.m_h  # make a number density
    #d *= 1e10 # 1/cm^2 to 1/km^2 
    d *= 10
    dx *= u.pc*1e-5 # pc to km 
    
    #print dx.shape, T.shape, v.shape, d.shape
    # Now loop over velocity bins
    if len(T) != 0:
        for iv in range(0,len(vscan)):
            sigma = np.sqrt(u.k_b*T/u.m_h)*1e-5  # in km/s 
            doppler = 1.0 - vscan[iv]/clight
            phi = invnorm*np.exp(-(0.5*(vscan[iv]*shape - v)/sigma)**2)/sigma
            I[iv] = facHI*np.sum(phi * d * dx)/doppler**2
    return I 

def yf(px1,py1,px0,py0,px):
    """ Given two points (px1,py1) and (px0,py0)
        this function returns the line that goes
        through those two points. 
    """
    return py0 + ((py1-py0)/(px1-px0))*(px-px0) 

def lv_proj_los(d, vr, vp, T, x, y, iner_frame=None, cart_grid=None):
    """ Create l,v projection plots via HI 21 cm emission
    ------- Inputs ---------
    d: density on Cartesian grid
    vr: radial velocity on Cartesian grid 
        (if cart_grid = True, vr = vx)
    vp: azimuthal velocity on Cartesian grid
        (if cart_grid = True, vp = vy)
    x,y: Cartisian coordinate arrays 
    
    ------- Outputs ----------
    l: longitude array
    v: line of sight velocity array
    intensity: intensity on l,v grid
    
    THIS ASSUMES CGS UNITS AS INPUT
    """
# Part 1: create data arrays    
    # set necessary constants
    x1min, x1max = x.min(), x.max()
    x2min, x2max = y.min(), y.max()
    nx1, nx2 = len(x), len(y)
    
    # Create meshgrids
    X, Y = np.meshgrid(x,y, indexing='xy')
    R    = np.sqrt(X*X + Y*Y)
    Phi  = np.arctan(Y/X)
    # Sun's position
    R0       = r0*np.ones((nx2,nx1))
    X0, Y0   = R0*np.cos(p0), R0*np.sin(p0)
    x0, y0   = r0*np.cos(p0), r0*np.sin(p0)

    # Distance from Sun to each grid point
    X1, Y1 = X-X0, Y-Y0
    R1     = np.sqrt(X1*X1 + Y1*Y1)
    # Unit vector pointing from Sun to each cell center 
    uX1, uY1 = X1/R1, Y1/R1
     
# Part 2: Get longitudes 
    # Determine longitude @ each point
    L = np.arccos((R0*R0 + R1*R1 - R*R)/(2.*R0*R1))
    # Correct longitdues that should be negative
    L[Phi < p0] = -L[Phi < p0]
    L[X   <  0] = -L[X   <  0]

# Part 3: Get line of sight velocities 
    # Correct vp velocities in rotating frame
    
    u = units_CGS()
    omgcgs = omg/u.myr # convert to 1/s
    Rcgs = R*u.pc # convert to cm 
    
    if iner_frame == False:
        if cart_grid == True:
            vr -= omg*Y
            vp += omg*X
        else:
            if np.mean(vp) > 0:
                vp += omgcgs*Rcgs
            else: 
                vp += omgcgs*Rcgs
    # Convert vr, vp from cm/s to km/s
    vr *= 1e-5 
    vp *= 1e-5 
    # Get x,y velocities from vr and vp
    if cart_grid == True:
        VX = vr
        VY = vp
    else: 
        VX, VY = vr*(X/R) - vp*(Y/R), vr*(Y/R) + vp*(X/R)
    # Get line of sight velocity (V "dot" u)
    Vlos  = VX*uX1 + VY*uY1
    # Correct for Sun's motion
    if cart_grid == True:
        Vlos -= V0*np.sin(L)
    else:
        if np.mean(vp) > 0:
            # check this
            Vlos -= V0*np.sin(L)
        else: 
            Vlos += V0*np.sin(L)
# Part 4: Interpolate arrays
    #d_interp    = interp2d(y,x,d   , kind='cubic', bounds_error=False,fill_value=0)
    #V_losinterp = interp2d(y,x,Vlos, kind='cubic', bounds_error=False,fill_value=0)
    d_interp    = RectBivariateSpline(y,x,d   )
    V_losinterp = RectBivariateSpline(y,x,Vlos)
    T_interp    = RectBivariateSpline(y,x,T   )
# Part 5: Create rays and velocity bins
   
    # note: nL and nV set the resolution of our intensity plot
    nL = int((L.max()    - L.min()   )/0.0024) #0.0044) 
    nV = int((Vlos.max() - Vlos.min())/1.5)     # 2.5   )
   
    L_rays = np.linspace(L.min()   ,L.max()   ,nL)
    V_bins = np.linspace(Vlos.min(),Vlos.max(),nV)
    
    # If we want a greater resolution than the simulation resolution
    # must redefine x (this becomes our sampling resolution)  
    # note this x goes from x.min() to r0  
    x = np.linspace(x.min(),r0, nL)
    
    # this is some geometry stuff required to create the rays
    pp = (p0 + np.pi/2) * np.ones(nL) - L_rays
    r2 = r0*np.sin(L_rays)
    # x,y point that makes a 90 deg angle with the center
    x2, y2 = r2*np.cos(pp), r2*np.sin(pp)
# Part 6: Compute intensity by looping over rays 
    I = np.zeros((nV, nL))
    t0 = time.time()
    for l in range(0, nL):
        # get the coordinates for the LoS
        ylos = yf(x2[l], y2[l], x0, y0, x)
        # get d and v along LoS
        dval =    d_interp.ev(ylos,x) #,grid=False)
        vval = V_losinterp.ev(ylos,x) #,grid=False)
        Tval =    T_interp.ev(ylos,x)
        # Put zeros in regions outside of the simulation box
        dval[abs(ylos) > x1max] = 0
        vval[abs(ylos) > x1max] = 0
        dval[x > x1max] = 0
        vval[x > x1max] = 0
        
        # LoS position array
        rlos       = np.sqrt( x**2 + ylos**2 )
        dxlos      = np.zeros(len(rlos)) 
        dxlos[-1 ] = abs(rlos[-1] - rlos[-2])
        dxlos[:-1] = abs(rlos[1:] - rlos[0:-1])
        # Do the radiative transfer
        I[:,l] = rad_trans(dval,vval,Tval,V_bins,dxlos)
# Part 7: Return results
    t1 = time.time()
    print "Part 6 takes:", t1 - t0, "s" 
    return L, Vlos, I



def lv_proj_bin(dens, vr, vp, x, y, iner_frame=None, cart_grid=None):
    """ Create l,v projection plots via a binning 
     routine given 
    ------- Inputs ---------
    dens: density on Cartesian grid
    vr: radial velocity on Cartesian grid 
        (if cart_grid = True, vr = vx)
    vp: azimuthal velocity on Cartesian grid
        (if cart_grid = True, vp = vy)
    x,y: Cartisian coordinate arrays 
    
    ------- Outputs ----------
    l: longitude array
    v: line of sight velocity array
    intensity: intensity on l,v grid
    """
# Part 1: create data arrays    
    # set necessary constants
    x1min, x1max = x.min(), x.max()
    x2min, x2max = y.min(), y.max()
    nx1, nx2 = len(x), len(y)
    
    # Create meshgrids
    X, Y = np.meshgrid(x,y, indexing='xy')
    R    = np.sqrt(X*X + Y*Y)
    Phi  = np.arctan(Y/X)
    # Sun's position
    R0       = r0*np.ones((nx2,nx1))
    X0, Y0   = R0*np.cos(p0), R0*np.sin(p0)

    # Distance from Sun to each grid point
    X1, Y1 = X-X0, Y-Y0
    R1     = np.sqrt(X1*X1 + Y1*Y1)
    # Unit vector pointing from Sun to each cell center 
    uX1, uY1 = X1/R1, Y1/R1
     
# Part 2: Get longitudes 
    # Determine longitude @ each point
    L = np.arccos((R0*R0 + R1*R1 - R*R)/(2.*R0*R1))
    # Correct longitdues that should be negative
    L[Phi < p0] = -L[Phi < p0]
    L[X   <  0] = -L[X   <  0]

# Part 3: Get line of sight velocities 
    # Correct vp velocities in rotating frame
    if iner_frame == False:
        if cart_grid == True:
            vr -= omg*Y
            vp += omg*X
        else:
            if np.mean(vp) > 0:
                vp += omg*R
            else: 
                vp += omg*R
    # Convert vr, vp from pc/Myr to km/s
    vr /= (3.154e13*3.24078e-14)
    vp /= (3.154e13*3.24078e-14)
    # Get x,y velocities from vr and vp
    if cart_grid == True:
        VX = vr
        VY = vp
    else: 
        VX, VY = vr*(X/R) - vp*(Y/R), vr*(Y/R) + vp*(X/R)
    # Get line of sight velocity (V "dot" u)
    Vlos  = VX*uX1 + VY*uY1
    # Correct for Sun's motion
    if cart_grid == True:
        Vlos -= V0*np.sin(L)
    else:
        if np.mean(vp) > 0:
            Vlos -= V0*np.sin(L)
        else: 
            Vlos += V0*np.sin(L)
# Part 4: Interpolate arrays (3rd order)
    d_interp    = RectBivariateSpline(y,x,dens)
    L_interp    = RectBivariateSpline(y,x,L   )
    V_losinterp = RectBivariateSpline(y,x,Vlos)
    R1_interp   = RectBivariateSpline(y,x,R1  )
    # Create sample points 
    n_samp = (x1max-x1min)/5  # this gives a sample point every 10 pc
    x_samp = np.linspace(x1min,x1max,n_samp)
    y_samp = np.linspace(x2min,x2max,n_samp)
    # Arrays at sample points 
    dinterp     = d_interp(y_samp,x_samp)
    Linterp     = L_interp(y_samp,x_samp)
    Vlosinterp  = V_losinterp(y_samp,x_samp)
    R1interp    = R1_interp(y_samp,x_samp)

    # Get area of each cell (on the sampled grid)
    area = (x_samp[1]-x_samp[0])*(y_samp[1]-y_samp[0])
    # Convert densities to mass
    dinterp *= area
# Part 5: Line of sight binning portion 
    # Create bins
    nbins_l = int((L.max() - L.min())/0.0044)  # gives 0.25 deg resolution
        
    L_bins = np.linspace(L.min()   , L.max()   , nbins_l)
    # Place longitude points into bins
    L_digit = np.digitize(Linterp, L_bins)
    
    # Place mass, velocity, position into line of sight bins
    t1 = time.time()
    # Create boolean array
    b_a = [ L_digit==q for q in xrange(0,nbins_l+1) ]
   
    l_arrays = [ [q,       dinterp[b_a[q]], 
                        Vlosinterp[b_a[q]],
                          R1interp[b_a[q]] ] for q in xrange(0,nbins_l+1) ] 
# Part 6: Velocity binning along lines of sight portion
    t2 = time.time()
    nbins_v = int((Vlos.max()-Vlos.min())/2.5) # gives 2.5 km/s resolution
    V_bins = np.linspace(Vlos.min(), Vlos.max(), nbins_v)
    # Create intensity array
    intensity = np.zeros((nbins_v-1,nbins_l-1))
    for j in xrange(0,nbins_l-1):
        # l_arrays[j][3]: distance to sun values along each line of sight
        # l_arrays[j][2]: velocity values along each line of sight
        # l_arrays[j][1]: density values along each line of sight
        # This portion of code bins the velocity values and takes the sum of mass of each bin
        #       along EACH line of sight 
       
        try:  
            myd, bin_edges, binnumber = binned_statistic(l_arrays[j][2], l_arrays[j][1], statistic='sum', bins = V_bins)
        # not sure if this position stuff is right
        #myp, bin_edges, binnumber = binned_statistic(l_arrays[j][2], l_arrays[j][3], statistic='mean',bins = V_bins)
        # write to intensity array
        except ValueError:
            myd = np.zeros(nbins_v-1)
        intensity[:,j] = myd 
    
    # this is just a "band-aid" fix for the NaN's in myp 
    intensity[intensity == np.nan] = 1e-20
    t3 = time.time()
    print "Part 5 takes:", t2-t1, "s, Part 6 takes:", t3-t2, "s" 
# Part 7: Return stuff
    return L, Vlos, intensity

#================================

#Read in system arguments
parser = argparse.ArgumentParser(description='Check for plotting flags')
parser.add_argument('start',nargs=1,type=int);
parser.add_argument('fin',nargs=1,type=int);
parser.add_argument('--quant'  ,dest='quant'  ,required=True)
parser.add_argument('--step'   ,dest='step'   ,default=1 ,type=int)
parser.add_argument('--lev'    ,dest='lev'    ,default=0 ,type=int)
parser.add_argument('--ival'   ,dest='ival'   ,default=-1,type=int)
parser.add_argument('--jval'   ,dest='jval'   ,default=-1,type=int)
parser.add_argument('--dval'   ,dest='dval'   ,default=-1,type=int)       #3D only
parser.add_argument('--view'   ,dest='view'   ,default=0, type=int)       #3D only
parser.add_argument('--units'  ,dest='units'  ,default=0,type=int)
parser.add_argument('--qminmax',dest='qminmax',default=-1,nargs='+',type=float)
parser.add_argument('--proj'   ,dest='proj'   ,action='store_true',default=False)  # 2D only
parser.add_argument('--proj_cart',dest='proj_cart',action='store_true',default=False) # 2D only
parser.add_argument('--cyl'    ,dest='cyl'    ,action='store_true',default=False)  # 2D only, cyl coords
parser.add_argument('--p2c'    ,dest='p2c'    ,action='store_true',default=False)  # 2D only 
parser.add_argument('--iner'   ,dest='iner'   ,action='store_true',default=False)  # 2D only
parser.add_argument('--log'    ,dest='log'    ,action='store_true',default=False)
parser.add_argument('--idlcol' ,dest='idlcol' ,action='store_true',default=False)
parser.add_argument('--column' ,dest='column' ,action='store_true',default=False)
parser.add_argument('--slice'  ,dest='slice'  ,action='store_true',default=False)
parser.add_argument('--dumpbin',dest='dumpbin',action='store_true',default=False)
parser.add_argument('--slide'  ,dest='slide'  ,action='store_true',default=False)
parser.add_argument('--gif'    ,dest='gif'    ,action='store_true',default=False)
parser.add_argument('--mp4'    ,dest='mp4'    ,action='store_true',default=False)
parser.add_argument('--eps'    ,dest='eps'    ,action='store_true',default=False)
parser.add_argument('--smrtrace',dest='smrtrace'  ,action='store_true',default=False)
parser.add_argument('--gridtrace',dest='gridtrace',action='store_true',default=False)
parser.add_argument('--nosmr',dest='nosmr',action='store_true',default=False)
args = parser.parse_args()

qmmflag = 1 if np.size(args.qminmax) >1 else 0


if args.idlcol:
    mycmap=IDL_color.readIDLctab('IDLctable/IDL_ctable_4')  
    mpl.cm.register_cmap(name='IDLcolor', cmap=mycmap)
    plt.rcParams['image.cmap'] = 'IDLcolor'
else:
    plt.rcParams['image.cmap'] = ctable

numlevs= len(params)
if (args.nosmr == True):
    numlevs=1

# Cylindrical stuff
if args.cyl or args.p2c or args.proj or args.quant == 'div_v' or args.quant == 'SII':
# Create log grid for plotting 
# Get face max/min values, and sizes 
    mn1 = params[0].x1min
    mx1 = params[0].x1max
    nx1 = int(params[0].nx1)
    mn2 = params[0].x2min
    mx2 = params[0].x2max
    nx2 = int(params[0].nx2)
    
# Create cell-faces array for x1
    ri  = np.zeros(nx1+1)
    dx1 = np.zeros(nx1+1)
# Create cell centered array for x1
    r   = np.zeros(nx1)
    ri[0] = mn1
    if int(params[0].ilog) == 0:
        x1rat = params[0].x1rat
    else:
        x1rat = np.power(mx1/mn1, 1.0/float(nx1))
    print 'x1rat = ', x1rat
    dx1[0] = (mx1-mn1)*(x1rat-1.0)/(x1rat**float(nx1)-1.0)
    for i in range(1,len(ri)):
        ri[i] = ri[i-1] + dx1[i-1]
        dx1[i] = dx1[i-1]*x1rat
    for i in range(0,len(r)):
        r[i] = ri[i] + 0.5*dx1[i]
# Create cell-faces array for x2 (necessary for pcolormesh)
    th = np.linspace(mn2,mx2,nx2)
# Create cell-centered array for x2
    dx2 = (mx2-mn2)/(float(nx2))
    th_cent = np.linspace(mn2+0.5*dx2,mx2-0.5*dx2,nx2)
  # define variables necessary for inertial frame
    if args.iner:
        T     = 2*np.pi/(omg)
        tvec  = np.linspace(0,T,nx2)

# p2c, proj stuff
if args.p2c or args.proj or args.quant == 'div_v' or args.quant == 'SII':
# Create x,y grid, we want cell centers for interpolation to Cartesian grid 
    mn1 = -max(r)
    mx1 =  max(r)
    mn2 = -max(r)
    mx2 =  max(r)
    if args.quant == 'div_v' or args.quant == 'SII': # give more resolution
        nx1 *= 4
        nx2 *= 4  
    x_grid = np.linspace(mn1,mx1,nx1)
    y_grid = np.linspace(mn2,mx2,nx2)
# proj_cart stuff
if args.proj_cart:
    mn1 = params[0].x1min
    mx1 = params[0].x1max
    nx1 = int(params[0].nx1)
    mn2 = params[0].x2min
    mx2 = params[0].x2max
    nx2 = int(params[0].nx2)
    # Create cell center array
    dx = (mx1-mn1)/(float(nx1))
    dy = (mx2-mn2)/(float(nx2))
    x_grid = np.linspace(mn1+0.5*dx,mx1-0.5*dx,nx1)
    y_grid = np.linspace(mn2+0.5*dx,mx1-0.5*dx,nx2)
    if args.iner:
        print "--proj_cart and --iner are not compatible currently! {'_'}"

# SII emission stuff
if args.quant == 'SII':
    SIIemis_interp = get_SIIemis(myline)

# Error checking
if args.p2c and args.cyl:
    print("Options '--p2c' and '--cyl' are not compatible with one another! {'_'}")
    quit()
if args.proj or args.proj_cart:
    if args.quant != 'proj':
        print("For projection plots MUST use '--quant=proj' {'_'}")
        quit()

if params[0].nx3 == 1:
    if args.view != 0:
        print "must have view = 0 for 2D simulations"
        quit()
    if  args.dval != -1:
        print "dval must not be set for 2D simulations"
        quit()

#Determine number of procs
nprocs=0
for file in os.listdir(athdir):
        if file.startswith("id"):
                nprocs += 1

#START LOOPING THROUGH ALL IMAGES
iters=0
top = int((args.fin[0]-args.start[0])/args.step)


if args.gridtrace:
    gridcols = ["hotpink","blue","gray","red","purple"]
if args.smrtrace:
    smrcols = ["hotpink","blue","gray","red","purple"]


#OPEN FIGURE
fig = plt.figure(figsize=(7.0,5.5),facecolor='white')

for imnum in range(top+1):
    starttime = time.time()
    tempnum    = args.start[0]+args.step*iters

    print 'image {} of {}'.format(imnum, top)
        
    for levs in range (args.lev,numlevs):
    
        minx = 1e20
        miny = 1e20
        minz = 1e20
        maxx = -1e20
        maxy = -1e20
        maxz = -1e20
        imgactive = 1

        resx     = int(params[levs].nx1)
        resy     = int(params[levs].nx2)
        resz     = int(params[levs].nx3)
        ntilesx1 = int(np.max([1,params[levs].ngridx1]))
        ntilesx2 = int(np.max([1,params[levs].ngridx2]))
        ntilesx3 = int(np.max([1,params[levs].ngridx3]))
        
        if args.proj or args.proj_cart or args.quant == 'SII':
            # we need 4 to deal with d,vr,vphi,T
            lev_ims1 = np.zeros((int(resz/ntilesx3),int(resy/ntilesx2),int(resx/ntilesx1),int(ntilesx1*ntilesx2*ntilesx3)))
            lev_ims2 = np.zeros((int(resz/ntilesx3),int(resy/ntilesx2),int(resx/ntilesx1),int(ntilesx1*ntilesx2*ntilesx3)))
            lev_ims3 = np.zeros((int(resz/ntilesx3),int(resy/ntilesx2),int(resx/ntilesx1),int(ntilesx1*ntilesx2*ntilesx3)))
            lev_ims4 = np.zeros((int(resz/ntilesx3),int(resy/ntilesx2),int(resx/ntilesx1),int(ntilesx1*ntilesx2*ntilesx3))) 
        elif args.quant == 'div_v':
            lev_ims1 = np.zeros((int(resz/ntilesx3),int(resy/ntilesx2),int(resx/ntilesx1),int(ntilesx1*ntilesx2*ntilesx3)))
            lev_ims2 = np.zeros((int(resz/ntilesx3),int(resy/ntilesx2),int(resx/ntilesx1),int(ntilesx1*ntilesx2*ntilesx3)))
        else: 
            lev_ims  = np.zeros((int(resz/ntilesx3),int(resy/ntilesx2),int(resx/ntilesx1),int(ntilesx1*ntilesx2*ntilesx3)))
        

        j=0
        if args.gridtrace:
            gxb = np.zeros((nprocs,2))
            gyb = np.zeros((nprocs,2))
            gzb = np.zeros((nprocs,2))
        for nps in range(np.max([1,nprocs])):
            file = get_file(athdir,base,nps,levs,tempnum,nprocs)

            if os.path.isfile(file):
                if args.proj or args.proj_cart or args.quant == 'SII':
                    temp_d, temp_vr, temp_vp, temp_T = get_quant(file,args.quant,args.units)
                    lev_ims1[:,:,:,j] = temp_d[:,:,:]
                    lev_ims2[:,:,:,j] = temp_vr[:,:,:]
                    lev_ims3[:,:,:,j] = temp_vp[:,:,:]
                    lev_ims4[:,:,:,j] = temp_T[:,:,:]
                elif args.quant == 'div_v':
                    temp_vr, temp_vp = get_quant(file,args.quant,args.units)
                    lev_ims1[:,:,:,j] = temp_vr[:,:,:]
                    lev_ims2[:,:,:,j] = temp_vp[:,:,:]
                else:
                    temp =get_quant(file,args.quant,args.units)
                    lev_ims[:,:,:,j] = temp[:,:,:]
                j +=1 
                
                # Cell centers
                x,y,z,simtime = get_xyzt(file,args.units)
                dx = x[1]-x[0]
                dy = x[1]-x[0]
                dz = x[1]-x[0]
                if args.cyl or args.proj or args.p2c:
                    minx =np.min([x[0]-0.5*dx1[0],minx]) 
                    maxx =np.max([x[len(x)-1]+0.5*dx1[-2],maxx]) 
                else:
                    minx =np.min([x[0]-0.5*dx,minx])
                    maxx =np.max([x[len(x)-1]+0.5*dx,maxy])   
                miny =np.min([y[0]-0.5*dy,miny])        
                maxy =np.max([y[len(y)-1]+0.5*dy,maxy]) 
                minz =np.min([z[0]-0.5*dz,minz])        
                maxz =np.max([z[len(z)-1]+0.5*dz,maxz])        
                if args.gridtrace:
                    gxb[nps,0] = x[0]        - 0.5*dx
                    gxb[nps,1] = x[len(x)-1] + 0.5*dx
                    gyb[nps,0] = y[0]        - 0.5*dy
                    gyb[nps,1] = y[len(y)-1] + 0.5*dy
                    gzb[nps,0] = z[0]        - 0.5*dz
                    gzb[nps,1] = z[len(z)-1] + 0.5*dz
        #RECONSTRUCT 3D CUBE
        p=0
        if args.proj or args.proj_cart or args.quant == 'SII':
            lev_im1 = np.zeros((resz,resy,resx))
            lev_im2 = np.zeros((resz,resy,resx))
            lev_im3 = np.zeros((resz,resy,resx))
            lev_im4 = np.zeros((resz,resy,resx))
            for k in range(ntilesx3):
                for j in range(ntilesx2):
                    for i in range(ntilesx1):
                        x0 = int(resx*(i*1./ntilesx1))
                        y0 = int(resy*(j*1./ntilesx2))
                        z0 = int(resz*(k*1./ntilesx3))
                        x1 = int(resx*(i+1.)/ntilesx1)
                        y1 = int(resy*(j+1.)/ntilesx2)
                        z1 = int(resz*(k+1.)/ntilesx3)
                        lev_im1[z0:z1,y0:y1,x0:x1] = lev_ims1[:,:,:,p]
                        lev_im2[z0:z1,y0:y1,x0:x1] = lev_ims2[:,:,:,p]
                        lev_im3[z0:z1,y0:y1,x0:x1] = lev_ims3[:,:,:,p]
                        lev_im4[z0:z1,y0:y1,x0:x1] = lev_ims4[:,:,:,p]
                        p +=1
        elif args.quant == 'div_v':
            lev_im1 = np.zeros((resz,resy,resx))
            lev_im2 = np.zeros((resz,resy,resx))
            for k in range(ntilesx3):
                for j in range(ntilesx2):
                    for i in range(ntilesx1):
                        x0 = int(resx*(i*1./ntilesx1))
                        y0 = int(resy*(j*1./ntilesx2))
                        z0 = int(resz*(k*1./ntilesx3))
                        x1 = int(resx*(i+1.)/ntilesx1)
                        y1 = int(resy*(j+1.)/ntilesx2)
                        z1 = int(resz*(k+1.)/ntilesx3)
                        lev_im1[z0:z1,y0:y1,x0:x1] = lev_ims1[:,:,:,p]
                        lev_im2[z0:z1,y0:y1,x0:x1] = lev_ims2[:,:,:,p]
                        p +=1
        else:
            lev_im = np.zeros((resz,resy,resx))
            for k in range(ntilesx3):
                for j in range(ntilesx2):
                    for i in range(ntilesx1):
                        x0 = int(resx*(i*1./ntilesx1))
                        y0 = int(resy*(j*1./ntilesx2))
                        z0 = int(resz*(k*1./ntilesx3))
                        x1 = int(resx*(i+1.)/ntilesx1)
                        y1 = int(resy*(j+1.)/ntilesx2)
                        z1 = int(resz*(k+1.)/ntilesx3)
                        lev_im[z0:z1,y0:y1,x0:x1] = lev_ims[:,:,:,p]
                        p +=1


        # dump combined binary file
        if args.dumpbin:
            binfile = athdir+base+'-'+args.quant+'-'+str(tempnum)+'-'+str(ntilesx1)+'x'+  \
            str(ntilesx2)+'x'+str(ntilesx3)+'.bin'
            binobj = open(binfile, mode='wb')
            
            resarr=np.array((resx,resy,resz))
            resarr.tofile(binobj)
            lev_im.tofile(binobj)
            binobj.close()

        #2D only 
        
        if args.proj or args.proj_cart or args.quant == 'SII':
            img_d  = lev_im1[0,:,:]
            img_vr = lev_im2[0,:,:] # recall img_vr = vx for cart sims
            img_vp = lev_im3[0,:,:] # recall img_vp = vy for cart sims
            img_T  = lev_im4[0,:,:] # temperature
        elif args.quant == 'div_v':
            img_vr = lev_im1[0,:,:]
            img_vp = lev_im2[0,:,:]
        else: 
            img = lev_im[0,:,:]
        #print minx,maxx,miny,maxy
        
        
        #TAKE THE LOG 
        if args.log and args.proj == False and args.quant != 'div_v' and args.quant != 'SII':
            inds = np.where(img == 0) 
            img[inds] = 1e-20
            img = np.abs(img)
            img = np.log10(img)

        # Uncomment (mn1,mx1,mn2,mx2) for original implementation by cfrazer
        if args.p2c == False and args.proj == False and args.cyl == False:
            mn1,mx1,mn2,mx2 = minx,maxx,miny,maxy   
        
        lab1,lab2 = "pc","pc"  # "[code length units]", "[code length units]"
       
    
        #PLOT IMAGE ON INITIAL AXIS
        if levs==args.lev:
            #determine data range
            if args.slide:
                if args.log:
                    inds = np.where(img > -19)
                    if(np.size(inds) == 0):
                        print "warning no cell values with quant > 0"
                        maxval = 1
                        minval = 0
                    else:
                        maxval = np.max(img[inds])
                        minval = np.min(img[inds])
                else:
                    maxval = np.max(img)
                    minval = np.min(img)

            elif imnum == 0:
                
                if qmmflag:
                    minval = args.qminmax[0]
                    maxval = args.qminmax[1]
                else:
                    if args.log and args.proj == False:
                        inds = np.where(img > -19)
                        maxval = np.max(img[inds])
                        minval = np.min(img[inds])
                    else:
                        if args.proj or args.proj_cart or args.quant == 'SII':
                            maxval = 0.0
                            minval = 0.0
                        elif args.quant == 'div_v':
                            maxval = np.max(img_vp)
                            minval = np.min(img_vp)
                        else:
                            maxval = np.max(img)
                            minval = np.min(img)

            
            valrange  = maxval-minval
            if valrange == 0:
                valrange = 0.5*maxval
                minval = minval - valrange
                maxval = maxval + valrange
            if valrange == 0:
                valrange = 0.5
                minval = minval 
                maxval = maxval + valrange


            fig.clear()
            
            if args.cyl:
                # Create cylindrical axes
                ax1, aux_ax1 = cyl_axes(fig,111,mn1,mx1,mn2,mx2)

                #------- FOR THE POLAR PLOT ---------------------------------------
                rad, theta = np.meshgrid(ri, th, indexing ='xy')
                
                # Check whether we want to do inertial rotation for comparison
                if args.iner:
                    img = iner_rot(img,tvec,T)
                
                im = aux_ax1.pcolormesh(theta,rad, img, vmin=minval, vmax=maxval) # shading='gourand')

                #timeunit = ' s' if args.units==1 or args.units==2 else ' [Myr]'#' [code units]' 
                
                timeunit = ' [Myr]' 
                ax1.set_title('t = '+'{0:.2f}'.format(simtime)+timeunit)
               
               # ADD color bar
                cbar= plt.colorbar(im,pad=0.1,label=get_title(args.quant,args.units,args.log,args.column))

            else:
                ax1 = fig.add_subplot(111)
                # Transform cyl to cartesian data
                if args.p2c:
                    # rotation for inertial frame 
                    if args.iner:
                        img = iner_rot(img,tvec,T)
                    # note that we act on the transpose of img and interpolate on cell centers!
                    # recall th_cent = \phi for cylindrical simulations
                    img = polar2cartesian(r, th_cent, img, x_grid, y_grid)
                # If doing proj plots 
                if args.proj:
                    if args.iner:
                        img_d  = iner_rot(img_d ,tvec,T)
                        img_vr = iner_rot(img_vr,tvec,T)
                        img_vp = iner_rot(img_vp,tvec,T)
                    # Now interpolate to cartesian grid
                    img_d  = polar2cartesian(r, th_cent, img_d , x_grid, y_grid)
                    img_vr = polar2cartesian(r, th_cent, img_vr, x_grid, y_grid)
                    img_vp = polar2cartesian(r, th_cent, img_vp, x_grid, y_grid)
                    img_T  = polar2cartesian(r, th_cent, img_T , x_grid, y_grid)
                    # Do the projection
                    L, Vlos, img = lv_proj_los(img_d, img_vr, img_vp, img_T, x_grid, y_grid, args.iner, args.proj_cart)
                    
                    # Get shape and extent of image
                    print np.shape(img)
                    myext = [L.min(),L.max(),Vlos.min(),Vlos.max()]
                    
                    # Take log of image 
                    if args.log:
                        img[img <= 0] = 1e-20
                        img = np.log10(img)
                    # Display the intensity
                    im =ax1.imshow(img,interpolation='None',
                                   extent=myext, origin='lower',
                                   aspect='auto',
                                   vmin = -6, vmax = -1)                     

                    ax1.set_xlim(-0.192,0.192)
                    ax1.set_ylim(Vlos.min(),Vlos.max())
                    ax1.set_xlabel("l [Rad]")
                    ax1.set_ylabel("v [km/s]")
                # If plotting div V
                elif args.quant == 'div_v':
                    
                    #img_vr = polar2cartesian(r, th_cent, img_vr, x_grid, y_grid)
                    #img_vp = polar2cartesian(r, th_cent, img_vp, x_grid, y_grid)

                    #img_vx, img_vy = vel_p2c(img_vr, img_vp, x_grid, y_grid)
                    
                    # Take the divergence 
                    #gradx = img_vx[2:,2:] - img_vx[:-2,:-2]
                    #grady = img_vy[2:,2:] - img_vy[:-2,:-2]
                   
                    grad_r = img_vr[1:,1:] - img_vr[1:,:-1] 
                    grad_p = img_vp[1:,1:] - img_vp[:-1,1:]

                    div_v = grad_r + grad_p
                    
                    div_v = polar2cartesian(r[1:], th_cent[1:], div_v, x_grid, y_grid)
                    #div_v = gradx + grady
                    #if args.units == 0:
                     #   u = units_CGS()
                     #   div_v *= u.v/1e5 # convert to km/s
                    #elif args.units == 1:
                     #   div_v /= 1e5  # convert to km/s
                    #else: 
                      #  div_v /= 1e3  # convert to km/s

                    #div_v[np.sqrt(x_grid[1:-1]**2 + y_grid[1:-1]**2) > 4000.] = 0.0

                    mn1 = x_grid.min()/1000
                    mx1 = x_grid.max()/1000
                    mn2 = y_grid.min()/1000
                    mx2 = y_grid.max()/1000

                    maxval = div_v.max()
                    minval = div_v.min()
                        
                    # only negative values contribute to emission
                    #div_v[div_v > 0] = 0.
                    #div_v = abs(div_v)

                    #div_v[div_v < 10] = 0.
                    #div_v[(div_v > 10) & (div_v < 20)] = 1.
                    #div_v[(div_v > 20) & (div_v < 30)] = 2.
                    #div_v[(div_v > 30) & (div_v < 40)] = 3.
                    #div_v[(div_v > 10) & (div_v < 100)] =  9.*(div_v[(div_v > 10) & (div_v < 100)]/div_v[(div_v > 10) & (div_v < 100)].max())
                    #div_v[div_v > 40] = 4. 
                        
                    # apply filter
                    pixel_size     = (x_grid[-1]-x_grid[0])/len(x_grid)
                    desired_smooth = 500. # pc
                    sigma = desired_smooth/pixel_size # gives sigma in pixels 
                    print sigma
                    
                    #div_v = gaussian_filter(div_v,sigma)  # gives smoothing 
                    
                    maxval = 2. # div_v.max()
                    minval = -2. #div_v.min()

                    im = ax1.imshow(div_v, interpolation='None',
                                            origin='lower',
                                            extent = [mn1,mx1,mn2,mx2],
                                            vmin=minval,vmax=maxval)
                    ax1.set_xlim(-1.5,1.5)
                    ax1.set_ylim(-.5,.5)
                    ax1.set_xlabel("x [kpc]")
                    ax1.set_ylabel("y [kpc]")
                
                # Plotting SII emission
                elif args.quant == 'SII': 
                    # Set units
                    if args.units == 1:
                        u = units_CGS()
                        vconv = 1.e-5
                    elif args.units == 2:
                        u = units_SI()
                        vconv = 1.e-3
                    else:
                        u = units_COMP()
                
                    # Take the log of these quantities
                    n     = img_d
                    n     = np.log10(n)
                    img_T = np.log10(img_T) 
                                                            
                    # Convert to cartesian grid   
                    n      = polar2cartesian(r, th_cent, n     , x_grid, y_grid)
                    '''
                    img_T  = polar2cartesian(r, th_cent, img_T , x_grid, y_grid)
                    
                    img_vr = polar2cartesian(r, th_cent, img_vr, x_grid, y_grid)
                    img_vp = polar2cartesian(r, th_cent, img_vp, x_grid, y_grid)

                    # Convert velocity to cartesian quantities 
                    img_vx, img_vy = vel_p2c(img_vr, img_vp, x_grid, y_grid)
                    
                    # Take the divergence 
                    gradx = img_vx[2:,2:] - img_vx[:-2,:-2]
                    grady = img_vy[2:,2:] - img_vy[:-2,:-2]

                    div_v = gradx + grady
                    '''
                    
                    grad_r = img_vr[1:,1:] - img_vr[1:,:-1] 
                    grad_p = img_vp[1:,1:] - img_vp[:-1,1:]

                    div_v = grad_r + grad_p
                    
                    div_v = polar2cartesian(r[1:], th_cent[1:], div_v, x_grid, y_grid)
                    
                    div_v *= vconv       # convert to [km/s], this will be v_shock 
                    div_v[div_v > 0] = 1.e-20 # only negative div_v counts towards emission 
                    div_v = abs(div_v) 

                    #----- DIFFERENT INTENSITY MEASURES --------# 

                    # Get T_shock (cf. Pg. 401 of Tielens (2005))
                    #Tshk = np.log10( 1.e5 * (div_v/1.e2)**2. ) 
                    Tshk = np.log10( (2.*(gam - 1.)/(gam + 1)**2.) * 
                                        (mu*u.m_h/u.k_b) * 
                                        (div_v*1.e5)**2. )   # eqn. 11.20 of Tielens (2005)

                    # 1) From interpolation table and T_shock
                    # Use n and T to compute intensity
                    int1 = 10.**(SIIemis_interp.ev(n, Tshk))*scaleh*u.pc # [erg cm^{-2} s^{-1} sr^{-1}]
                    int1 = np.log10(int1)
                    
                    # 2) From a simple prescription using div_v 
                    '''
                    div_v = abs(div_v)
                    
                    int2  = np.zeros(div_v.shape)
                    int2[div_v < 10] = 0.
                    int2[(div_v > 10) & (div_v < 20)] = 1.
                    int2[(div_v > 20) & (div_v < 30)] = 2.
                    int2[(div_v > 30) & (div_v < 40)] = 3.
                    int2[div_v > 40] = 4. 
                    # 3) Intensity \propto v_s ^3 
                    int3 = np.power(div_v,3.)
                    int3[int3 == np.nan] = 1e-20 
                    #int3 = np.log10(int3)
                    '''
                    #------------------------------------------#
                    img = int1

                    # apply filter
                    pixel_size = (x_grid[-1]-x_grid[0])/len(x_grid)
                    sigma = desired_smooth/pixel_size # gives sigma in pixels 
                    print sigma
                    
                    #img = gaussian_filter(img,sigma,truncate=2,cval=1.0)  # gives smoothing 
                    #img = np.log10(img)

                    # set min, max for plotting               
                    #mask = np.isnan(img)
                    #img[mask] = -20
                    #maxval = img[mask==False].max() 
                    # HERE
                    #minval = 3 #img[mask==False].min()

                    minval = -15
                    maxval = 0
                    
                    # set extent for plotting 
                    mn1 = x_grid.min()/1000
                    mx1 = x_grid.max()/1000
                    mn2 = y_grid.min()/1000
                    mx2 = y_grid.max()/1000

                    im = ax1.imshow(img, interpolation='None',
                                            origin='lower',
                                            extent = [mn1,mx1,mn2,mx2],
                                            vmin=-15,vmax=0,
                                            aspect= 0.77 )

                    # HERE
                    ax1.set_xlim(-1.5,1.5)
                    ax1.set_ylim(-.5,.5)
                    ax1.set_xlabel("x [kpc]")
                    ax1.set_ylabel("y [kpc]")

                elif args.proj_cart:
                    # Do the projections
                    L, Vlos, img = lv_proj1(img_d, img_vr, img_vp, x_grid, y_grid, args.iner, args.proj_cart)
                    print np.shape(img)
                    myext = [L.min(), L.max(), Vlos.min(), Vlos.max()]
                    # Take log of image
                    img[img <= 0] = 1e-20
                    img = np.log10(img)

                    # Display
                    im = ax1.imshow(img,interpolation='None',
                                    origin='lower',
                                    aspect='auto',
                                    extent=myext, 
                                    vmin = -1.5, vmax = 1)
                    ax1.set_xlim(-0.192,0.192)
                    ax1.set_ylim(Vlos.min(),Vlos.max())
                    ax1.set_xlabel("l [Rad]")
                    ax1.set_ylabel("v [km/s]")
                
                # Not proj plots
                else: 
                   im = ax1.imshow(img,interpolation='None',extent=[mn1,mx1,mn2,mx2],origin='lower',vmin=minval,vmax=maxval)
                   ax1.set_xlim(mn1,mx1)
                   ax1.set_ylim(mn2,mx2)
                   ax1.set_xlabel(lab1) 
                   ax1.set_ylabel(lab2) 
                
                
               # timeunit = ' s' if args.units==1 or args.units==2 else ' [Myr]'#' [code units]'
                timeunit = ' [Myr]' 
                ax1.set_title('t = '+'{0:.2f}'.format(simtime)+timeunit)
                           
                 #ax1.set_title(get_title(quant,units)) 
               #ADD COLORBAR (using make axes locateable)
                divider = make_axes_locatable(ax1)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                cbar= plt.colorbar(im, cax=cax,label = get_title(args.quant,args.units,args.log,args.column))
                #plt.gca().invert_yaxis()   

                minx0 = minx
                maxx0 = maxx
                miny0 = miny
                maxy0 = maxy
            
        elif(imgactive ==1):
            
            im = ax1.imshow(img,interpolation='None',extent=[mn1,mx1,mn2,mx2],origin='lower',vmin=minval,vmax=maxval)
            plt.gca().invert_yaxis()
            ax1.set_xlim(minx0,maxx0)
            ax1.set_ylim(miny0,maxy0)
      
            if args.smrtrace:
                if args.view ==0:
                    xdom_bounds = np.array([minx,minx,maxx,maxx,minx])
                    ydom_bounds = np.array([miny,maxy,maxy,miny,miny])
                    ax1.plot(xdom_bounds,ydom_bounds,color=smrcols[levs],marker='None',linewidth=2)
                if args.view ==1:
                    xdom_bounds = np.array([minx,minx,maxx,maxx,minx])
                    zdom_bounds = np.array([minz,maxz,maxz,minz,minz])
                    ax1.plot(xdom_bounds,zdom_bounds,color=smrcols[levs],marker='None',linewidth=2)
                if args.view ==2:
                    ydom_bounds = np.array([miny,miny,maxy,maxy,miny])
                    zdom_bounds = np.array([minz,maxz,maxz,minz,minz])
                    ax1.plot(ydom_bounds,zdom_bounds,color=smrcols[levs],marker='None',linewidth=2)
        if args.gridtrace:
            for nps in range(np.max([1,nprocs])):
                if args.view == 0:
                    xgrid_bounds = np.array([gxb[nps,0],gxb[nps,0],gxb[nps,1],gxb[nps,1],gxb[nps,0]])
                    ygrid_bounds = np.array([gyb[nps,0],gyb[nps,1],gyb[nps,1],gyb[nps,0],gyb[nps,0]])
                    ax1.plot(xgrid_bounds,ygrid_bounds,color=gridcols[levs],marker='None')
                if args.view == 1:
                    xgrid_bounds = np.array([gxb[nps,0],gxb[nps,0],gxb[nps,1],gxb[nps,1],gxb[nps,0]])
                    zgrid_bounds = np.array([gzb[nps,0],gzb[nps,1],gzb[nps,1],gzb[nps,0],gzb[nps,0]])
                    ax1.plot(xgrid_bounds,zgrid_bounds,color=gridcols[levs],marker='None')
                if args.view == 2:
                    ygrid_bounds = np.array([gyb[nps,0],gyb[nps,0],gyb[nps,1],gyb[nps,1],gyb[nps,0]])
                    zgrid_bounds = np.array([gzb[nps,0],gzb[nps,1],gzb[nps,1],gzb[nps,0],gzb[nps,0]])
                    ax1.plot(ygrid_bounds,zgrid_bounds,color=gridcols[levs],marker='None')

        if args.eps ==1 and levs==numlevs-1:
            print "HERE"
            plt.savefig("int_smooth"+str(desired_smooth)+"_bar6kpc.eps",format='eps',dpi=1000,bbox_inches='tight')  
            #plt.savefig("div_v_bar6kpc.eps",format='eps',dpi=1000,bbox_inches='tight')  

        if (args.gif == 1 or args.mp4 == 1) and imnum == 0 : 
            test= os.path.exists(athdir+'snapshots/')
            if test == False:
                os.makedirs(athdir+'snapshots/')  
            if test == True:       
                filelist = [f for f in os.listdir(athdir+'snapshots/') if f.endswith(".png") ]
                for f in filelist:
                    os.remove(athdir+'snapshots/'+f)

        if imnum < top:
            if levs == numlevs-1:
                if args.gif == 1 or args.mp4 == 1:
                    plt.savefig(athdir+"snapshots/"+base+get_num(tempnum)+'.png')
                else:
                    plt.show(block=False)  
                    plt.draw()
                    plt.pause(0.001)
        else:
            if levs == numlevs-1:
                if args.gif == 1 or args.mp4 == 1:
                    plt.savefig(athdir+"snapshots/"+base+get_num(tempnum)+'.png')
                else:
                    plt.show(block=True)

        


    iters +=1

    endtime = time.time()
    totaltime = starttime-endtime
        
if args.gif:                               
    print "Converting png's to gif"
    os.system("convert -loop 0 -delay 0 " +athdir+"snapshots/*.png "+athdir+base+"_"+args.quant+".gif")
    filelist = [f for f in os.listdir(athdir+'snapshots/') if f.endswith(".png") ]
    for f in filelist:
        os.remove(athdir+'snapshots/'+f)
    os.rmdir(athdir+"snapshots/")
    print "Gif completed"
if args.mp4:
	print "Converting png's to mp4 movie"
	os.system("cat "+athdir+"snapshots/*.png | ffmpeg -y -framerate 10 -f image2pipe -i - "+athdir+base+"_"+args.quant+".mp4")
	filelist = [f for f in os.listdir(athdir+'snapshots/') if f.endswith(".png") ]
	for f in filelist: 
		os.remove(athdir+'snapshots/'+f)
	os.rmdir(athdir+"snapshots/")
	print "mp4 movie created"


