#!/usr/bin/python
import numpy as np
import sys
import os
import subprocess as sbp
from functools import reduce
# matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.projections import PolarAxes
import mpl_toolkits.axisartist.floating_axes as floating_axes
from mpl_toolkits.axisartist.grid_finder import DictFormatter
from matplotlib import _cntr as cntr 
# scipy
from scipy.interpolate import griddata 
from scipy.interpolate import interp1d 
import argparse
from argparse import RawTextHelpFormatter
import time as t

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
from units_class import * 
from read_bin import read_bin 
import read_athinput


#=====================================================
#
#  Code: plot_sims.py
#
#  Purpose: Reads in athena binary dumps and plots
#           quantities. Designed to replace 
#           plot_cylbgsbu.py eventually  
#
#  Keywords: python plot_sims.py -h   
#
#  Usage: python plot_sims.py quant   
#
#  WARNING: THIS MUST BE RUN FROM SIMULATION DIRECTORY 
#          
#  Author: John Dupuy 
#          UNC Chapel Hill
#  Date:    02/15/18 
#  Updated: 03/01/18 
#=====================================================


#============ FUNCTIONS ==============================

#------------ FUNCTIONS FOR READING DATA -------------#
#\func get_files()
# returns dictionary of the form {'id0': bin_files0, 'id1': bin_files1,...}
def get_files():
    # Get processor directories 
    proc  = os.walk('.').next()[1]
    # Get rid of unwanted directories 
    proc  = [x for x in proc if 'id' in x]
    nproc = len(proc)
    # Define processor range
    p_no  = range(nproc)
    
    # Now create the dictionary 
    pdict = {} 

    #   id0 directory is special since it contains other stuff
    pdict['id0'] = sorted([fname for fname in
                           os.listdir('id0'+'/')
                           if  '.bin' in fname
                           and 'lv'  not in fname
                           and 'sii' not in fname
                           and 'otf' not in fname])
    for i in p_no[1:]:
        pno = 'id'+str(i) 
        try:
            pdict[pno] = sorted([fname for fname in
                                 os.listdir(pno+'/')
                                 if '.bin' in fname])
        except OSError:
            print('[get_files]: ERROR, sim files not present for %s ' % (pno))  
            quit() 

    return pdict 

#\func get_quant 
# gets quantity 
def get_quant(file,quant,units,precision=32,omg=0.064424,ctrace=False):
    #read in binary file
    nx1,nx2,nx3,x1,x2,x3,\
    d,M1,M2,M3,e,ie,s,\
    bx,by,bz,phi, \
    gamm1,cs,t,dt,nscalars = read_bin(file,precision)
    
    # Set units
    
    if   units == 1:
        u = units_CGS()
    elif units == 2:
        u = units_SI()
    elif units == 0:
        u = units_COMP()

    ucgs = units_CGS() 

    # Do the unit conversion
    '''
    e  *= u.esdens
    ie *= u.esdens
    d  *= u.rhos
    M1 *= u.momsdens
    M2 *= u.momsdens
    M3 *= u.momsdens
    cs *= u.v
    '''
    
    # Deal with potentially empty nscalar
    if s.size:
        s0 = s[:,:,:,0]
        s1 = s[:,:,:,1]
    else: 
        s0 = np.zeros( d.shape )
        s1 = np.zeros( d.shape ) 
    
    # Define dictionary for quants
    all_quants = {'E':e, 'ie':ie, 'd':d,
                  'n':d*ucgs.rhos/(ucgs.m_h), 'pie':gamm1*ie,
                  'p':gamm1*(e - 0.5*(M1**2.+M2**2.+M3**2.)/d),
                  'T':gamm1*ie*ucgs.esdens/(ucgs.k_b * ( d*ucgs.rhos/ucgs.m_h)),
                  'M':np.sqrt(M1**2.+M2**2.+M3**2.),
                  'v1':M1/d,'v2':M2/d,'v3':M3/d,
                  'M1':M1,'M2':M2,'M3':M3,
                  'V':np.sqrt(M1**2.+M2**2.+M3**2.)/d,
                  'cs': gamm1*ie/d, 
                  's1': s1}
    if ctrace:
        return t, x1, x2, all_quants[quant], all_quants['s1']
    else:
        return t, x1, x2, all_quants[quant] 

#\func get_athinput():
# reads athinput and returns base, params
def get_athinput():
    athin = os.getcwd()+'/'+[fnm for fnm in os.listdir(os.getcwd()+'/')
                                        if fnm.startswith('athinput')][0]
    base, params = read_athinput.readath(athin) 
    return base, params    

#\func get_stitch
# given pdict, this will stitch together an image
# Here pdict is the dictionary from get_files()
def get_stitch(pdict,quant,myfrms,units=0,ctrace=False):
    
    # read athinput
    base, params = get_athinput() 
   
    # get no. of points 
    nx1 = int(params[0].nx1)
    nx2 = int(params[0].nx2)
    # get extent of grid
    mn1 = params[0].x1min
    mx1 = params[0].x1max
    mn2 = params[0].x2min
    mx2 = params[0].x2max
    # get no. of processors
    nproc = len(pdict) 
    # get no. of processors in x1 dir
    ngridx1 = int(params[0].ngridx1)
    ngridx2 = int(params[0].ngridx2) 
    # get size of grid on each processor
    npx1    = nx1/ngridx1
    npx2    = nx2/ngridx2 
    # get no. of files in each processor
    #nf    = len(pdict['id0']) # assumes all are same 
    nf     = len(myfrms) 

    # Make array to hold quant for all time
    quant_arr = np.zeros((nf, nx2, nx1))
    if ctrace:
        c_arr = np.zeros((nf, nx2, nx1)) 
    # Make time, x1, x2 arrays 
    tarr = np.zeros(nf)
    x1   = np.zeros(nx1)
    x2   = np.zeros(nx2) 

    # Define processor index
    ip = 0
    for ip2 in range(ngridx2):
        fac2     = ip2%ngridx2 
        for ip1 in range(ngridx1): 
            fac1     = ip1%ngridx1
            # Define processor key
            pkey     = 'id'+str(ip)
            i        = 0            # index for filling arrays
            for iff in myfrms: 
                # Get filename
                fnm      = pkey + '/' + pdict[pkey][iff]
                # Fill data array
                data = get_quant(fnm,quant,units,omg=params[0].omg,ctrace=ctrace) 
                # parse data 
                tarr[i]                                = data[0]
                x1[fac1*npx1:npx1*(1+fac1)]            = data[1]
                x2[fac2*npx2:npx2*(1+fac2)]            = data[2]
                quant_arr[i,fac2*npx2:npx2*(1+fac2),
                              fac1*npx1:npx1*(1+fac1)] = data[3]
                if ctrace:
                    c_arr[i,fac2*npx2:npx2*(1+fac2),
                              fac1*npx1:npx1*(1+fac1)] = data[4]
                i += 1
            ip += 1

    # Note x1 is converted to kpc
    if ctrace:
        aprc, appc = get_fc(False, True) 
        c_arr[:] *= aprc               # get cloud mass  
        return tarr, x1/1.e3, x2, quant_arr, c_arr
    else:
        return tarr, x1/1.e3, x2, quant_arr
#------------------------------------------------------#

#---------- MISCELLANIOUS FUNCTIONS -------------------# 
# \func get_fc()
# given cell centered x1, x2, returns face centered grid 
# THIS ASSUMES ilog=1 in x1-dir (i.e. uniform log spacing) 
def get_fc(cart, apc=False):
    base, params = get_athinput() 
   
    if cart:
        x1f = np.linspace(params[0].x1min/1.e3,
                          params[0].x1max/1.e3,params[0].nx1+1)
        x2f = np.linspace(params[0].x2min/1.e3,
                          params[0].x2max/1.e3,params[0].nx2+1) 
    else:
        # x1 direction (convert to kpc) 
        lx1f = np.linspace(np.log(params[0].x1min/1.e3),
                           np.log(params[0].x1max/1.e3),params[0].nx1+1)
        # x2 direction
        x2f  = np.linspace(       params[0].x2min,
                                  params[0].x2max, params[0].nx2+1) 
        x1f  = np.exp(lx1f) 
    if apc:      
        # get area per cell 
        x1c = 0.5*(x1f[1:] + x1f[0:-1])*1.e3
        ar  =     (x1f[1:] - x1f[0:-1])*1.e3  # convert to pc 
        ap  = x1c*(x2f[1:] - x2f[0:-1])

        # Area per cell
        return np.meshgrid(ar*ap, ar*ap, indexing='xy') 
    else:
        return np.meshgrid(x1f, x2f, indexing='xy')   

# \func get_fc2()
# gets fc2 values with nx1 cells, not correct, 
#   only used for cloud tracing,
#   so that full range of 2pi is plotted
# THIS ASSUMES ilog=1 in x1-dir (i.e. uniform log spacing) 
def get_fc2(cart):
    base, params = get_athinput() 
    
    if cart:
        x1f = np.linspace(params[0].x1min/1.e3,
                          params[0].x1max/1.e3,params[0].nx1)
        x2f = np.linspace(params[0].x2min/1.e3,
                          params[0].x2max/1.e3,params[0].nx2) 
    else:
        # x1 direction (convert to kpc) 
        lx1f = np.linspace(np.log(params[0].x1min/1.e3),
                           np.log(params[0].x1max/1.e3),params[0].nx1)
        # x2 direction
        x2f  = np.linspace(       params[0].x2min,
                                  params[0].x2max, params[0].nx2) 
        x1f  = np.exp(lx1f) 
    
    return np.meshgrid(x1f, x2f, indexing='xy') 

# \func get_levels() 
# gets the levels for the contour tracing of the cloud
# This returns levels that draw contours encompassing given
# percentages of the cloud mass
def get_levels(img,pcts):
    # img : array of cloud mass 
    # pcts: array of desired percentages (in descending order) 

    # define 1D array of values inside img 
    n    = 1000
    t    = np.linspace(0, img.max(), n)  
    # Compute integral within each t 
    integral = (( img >= t[:, None, None]) * img).sum(axis=(1,2)) 
    # Interpolate on the integral 
    f        = interp1d(integral, t)
    try: 
        levels = list(f(pcts))
    except ValueError: 
        levels = [1.0]  
        
    return levels   

#------------------------------------------------------#

#---------- PLOTTING FUNCTIONS ------------------------#
#\func polar2cartesian
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

def get_labels(quant,iunit=0,log=False):

    # Define dictionary for quants
    lab_quants = {'E':'Energy density', 'ie':'Internal energy density', 
                  'd':'$\Sigma$',
                  'n':'Column density', 'p':'Surface pressure',
                  'pie':'Surface pressure (from U.IE)', 
                  'T':'T',
                  'M':'M$_{tot}$',
                  'v1':'v$_R$','v2':'v$_{\phi}$','v3':'v$_z$',
                  'M1':'$\Sigma$ v$_R$','M2':'$\Sigma$ v$_{\phi}$','M3':'$\Sigma$ v$_{z}$',
                  'cs':'c$_s$', 
                  'V':'V$_{tot}$',
                  's1':'$\Sigma_c$'}

    xlabel = 'x [kpc]'
    ylabel = 'y [kpc]' 
    cbar_l = lab_quants[quant]

    # Define dictionary for units 
    lab_unitsCOMP = {'E':' [M$_{\odot}$ Myr$^{-2}$]','ie':' [M$_{\odot}$ Myr$^{-2}$]',
                     'd':' [M$_{\odot}$ pc$^{-2}$]', 'p':' [M$_{\odot}$ Myr$^{-2}$]',
                     'pie':' [M$_{\odot}$ Myr$^{-2}$]',
                     'T':' [K]', 'M': ' [M$_{\odot}$ pc$^{-1}$ Myr$^{-1}]$',
                     'v1':' [pc Myr$^{-1}$]', 'v2':' [pc Myr$^{-1}$]','v3':' [pc Myr$^{-1}$]', 
                     'M1': ' [M$_{\odot}$ pc$^{-1}$ Myr$^{-1}]$', 'M2': ' [M$_{\odot}$ pc$^{-1}$ Myr$^{-1}]$', 
                     'M3': ' [M$_{\odot}$ pc$^{-1}$ Myr$^{-1}]$','V':' [pc Myr$^{-1}$]',
                     'cs': ' [pc Myr$^{-1}$]',
                     's1': ' [M$_{\odot}$ pc$^{-2}$]'}

    lab_unitsCGS  = {'E':' [g s$^{-2}$]','ie':' [g s$^{-2}$]',
                     'd':' [g cm$^{-2}$]', 'p':' [g s$^{-2}$]',
                     'pie':' [g s$^{-2}$]',
                     'T':' [K]', 'M': ' [g cm$^{-1}$ s$^{-1}]$',
                     'v1':' [cm s$^{-1}$]', 'v2':' [cm s$^{-1}$]','v3':' [cm s$^{-1}$]', 
                     'M1': ' [g cm$^{-1}$ s$^{-1}]$', 'M2': ' [g cm$^{-1}$ s$^{-1}]$', 
                     'M3': ' [g cm$^{-1}$ s$^{-1}]$', 'V':' [cm s$^{-1}$]',
                     'cs': ' [cm s$^{-1}$]',
                     's1': ' [g cm$^{-2}$]'}
   
    lab_unitsSI   = {'E':' [kg s$^{-2}$]','ie':' [kg s$^{-2}$]',
                     'd':' [kg m$^{-2}$]', 'p':' [kg s$^{-2}$]',
                     'pie':' [kg s$^{-2}$]',
                     'T':' [K]', 'M': ' [kg m$^{-1}$ s$^{-1}]$',
                     'v1':' [m s$^{-1}$]', 'v2':' [m s$^{-1}$]','v3':' [m s$^{-1}$]', 
                     'M1': ' [kg m$^{-1}$ s$^{-1}]$', 'M2': ' [kg m$^{-1}$ s$^{-1}]$', 
                     'M3': ' [kg m$^{-1}$ s$^{-1}]$', 'V':' [m s$^{-1}$]',
                     'cs': ' [m s$^{-1}$]',
                     's1': ' [kg m$^{-2}$]'}

    # Handle units 
    if iunit == 0:
        cbar_l += lab_unitsCOMP[quant]
    elif iunit == 1:
        cbar_l += lab_unitsCGS[quant]
    elif iunit == 2:
        cbar_l += lab_unitsSI[quant] 

    if log:
        cbar_l = 'log$_{10}$('+cbar_l+')'

    return xlabel, ylabel, cbar_l 

#\func factors()
# given n, this returns the factors 
def factors(n): 
    return sorted(reduce(list.__add__,
                 ([i,n//i] for i in 
                 range(1,int(pow(n,0.5)+1)) if n%i == 0)))

#------------- MAIN FUNCTION ------------------------#
def main():

    ctable = 'magma'
    plt.rcParams['image.cmap'] = ctable
    base, params = get_athinput() 
    mc           = params[0].mhvc


    # Read in system arguments
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter) 
    parser.add_argument("quant",type=str,
                        help="Plotting options:\n"
                              " Any quant available in get_quant\n")
    parser.add_argument("--anim", dest="anim",action='store_true',
                        default=False,
                        help="Switch to do animation\n")
    parser.add_argument("--iani", dest="iani",nargs=2,required=False,
                        default=[0,0],type=int,
                        help="Animate from frame iani[0] to iani[1]\n")
    parser.add_argument("--qminmax", dest="qminmax",nargs='+',required=False,
                        default=-1,type=float,
                        help="Min/max value for imshow")
    parser.add_argument("--ifrm", dest="ifrm",type=int,default=[0],
                        nargs='+',
                        help="Frame of simulation to plot:\n"
                             "  0: tsim = 0 \n"
                             "  1: tsim = dt_dump\n"
                             "  2: tsim = 2*dt_dump\n"
                             "  .               \n"
                             "  .               \n"
                             "  .               \n"
                             "  n: tsim = n*dt_dump\n"
                             " Note: by entering multiple integers, will make\n"
                             " panel plots for each ifrm")
    parser.add_argument("--mnmx",dest="mnmx", type=float,nargs=2,
                        required=False,default=[-3.5,3.5],
                        help="Plotting range in x and y\n"
                             "Note: assumes x, y range are equivalent")
    parser.add_argument("--save", dest="save",action='store_true',
                        default=False,
                        help="Switch to save anim or figure")
    parser.add_argument("--log",dest="log",action='store_true',
                        default=False,
                        help="Switch to take log images, default is False")
    parser.add_argument("--ctrace",dest="ctrace",action='store_true',
                        default=False, help="Switch to overplot a contour for the cloud location\n")
    parser.add_argument("--units",dest="units",type=int,required=False,
                        default=0, help="units: 0-comp,1-cgs,2-SI")  
    parser.add_argument("--levels",dest="levels",type=int,required=False,
                        default=1,help="no. oflevels for cloud tracing")
    parser.add_argument("--cart", dest="cart",action='store_true',
                        default=False, help="Switch for cartesian simulations") 
    parser.add_argument("--grid", dest="grid",action='store_true',
                        default=False, help="Switch to make plot to show grid")
    parser.add_argument("--fmt", dest="fmt", default='eps',
                        type=str, help='format for saving graphics, default: eps') 

    
    # parsing arguments            
    args  = parser.parse_args()
    quant = args.quant
    anim  = args.anim
    iani  = args.iani
    ifrm  = args.ifrm
    save  = args.save
    log   = args.log
    mnx, mxx = args.mnmx
    iunit    = args.units 
    qminmax  = args.qminmax
    ctrace   = args.ctrace
    levels   = args.levels
    cart     = args.cart
    grid     = args.grid
    fmt      = args.fmt
    # Get qminmax flag 
    qflag = True if np.size(qminmax) > 1 else False
    # Get panel flag
    pflag = True if np.size(ifrm) > 1 else False 
    
    if np.size(ifrm) == 1: ifrm = ifrm[0]

    
    # get files 
    pdict = get_files() 
    # Change default iani values
    if iani[1] == 0:
        iani[1] = len(pdict['id0']) 

    # determine myframes
    if anim:
        myfrms = range(iani[0],iani[1])
    elif pflag:
        myfrms = ifrm
    else: 
        myfrms = [ifrm] 
    # get data
    if ctrace:
        tarr, x1, x2, imgs, imgc = get_stitch(pdict,quant,myfrms,iunit,ctrace=ctrace)  
        # Normalize imgc by hvc mass
        imgc /= mc 
        # Set contour parameters for ctrace
        colors=['#40E0D0','#00C957']#,'#FF3030']
        if pflag:
            lw = 1.
        else:
            lw = 2.
        pcts = np.array([0.95,0.5])  # define percentages  
        alpha=1.0
    else:
        tarr, x1, x2, imgs       = get_stitch(pdict,quant,myfrms,iunit)

    # get face centered meshgrid
    x1f, x2f           = get_fc(cart) 
    # get cartesian version of face centered meshgrid
    if cart:
        x1cf, x2cf = x1f, x2f
    else:
        x1cf, x2cf = x1f*np.cos(x2f), x1f*np.sin(x2f) 
    # get cartesian version of 'cell centered' meshgrid
    x1c, x2c           = get_fc2(cart)
    if cart:
        x1cc, x2cc = x1c, x2c
    else:
        x1cc, x2cc = x1c*np.cos(x2c),  x1c*np.sin(x2c)
    # get total number of images
    n = len(imgs) 
    # get labels for plots
    xlab, ylab, clab = get_labels(quant,iunit,log) 

    # Conflicting argument checks 
    if anim and ifrm != 0:
        print("[main]: specifying a single frame while animating doesn't make sense")
        quit()

    if anim:
        print("\n[main]: Animating from t = %1.1f [Myr]\n" 
              "                    to t = %1.1f [Myr]\n"
                    %( tarr[0], tarr[-1] ) )  
    # Take log of data
    if log: 
        imgs[imgs < 0] = 1e-20
        imgs           = np.log10(imgs)  
   
    if not pflag: 
        # Open figure 
        fig = plt.figure(figsize=(7.5,7.0),facecolor='white')
        ax1 = fig.add_subplot(111) 
            # Set labels for x,y
        ax1.set_xlabel(xlab)
        ax1.set_ylabel(ylab)
            # Set xmin, xmax
        ax1.set_xlim(mnx,mxx)
        ax1.set_ylim(mnx,mxx)
            # Set aspect ratio
        ax1.set_aspect('equal') 
            # Colorbar stuff
        if not grid:
            div = make_axes_locatable(ax1)
            cax = div.append_axes('right', '5%', '5%')


   
    # Handle animation
    if anim: 
        if qflag:
            qmin, qmax = qminmax[0], qminmax[1]
        else:
            qmin, qmax = np.min(imgs), np.max(imgs) 

        ax1.set_title('t = %1.1f [Myr]' %(tarr[0]) )
        im   = ax1.pcolorfast(x1cf,x2cf, imgs[0], vmin=qmin,vmax=qmax)
        im.set_rasterized(True) 
        if ctrace:
            if tarr[iani[0]] > params[0].thvcs:
                imc  = ax1.contour(x1cc,x2cc,imgc[0],levels=get_levels(imgc[0],pcts)
                                                    ,colors=colors,linewidths=lw,alpha=alpha)
        cbar = fig.colorbar(im,label=clab,cax=cax) 
        

        def animate(ifrm):
            # Clear figure 
            ax1.cla()
            cax.cla()
            
            # Set title 
            ax1.set_title('t = %1.1f [Myr]' %(tarr[ifrm]) )
            # Plot 
            im   = ax1.pcolorfast(x1cf,x2cf, imgs[ifrm], vmin=qmin,vmax=qmax)
            im.set_rasterized(True) 
            if ctrace:
                if tarr[ifrm] > params[0].thvcs:
                    imc  = ax1.contour(x1cc,x2cc,imgc[ifrm],levels=get_levels(imgc[ifrm],pcts),
                                                            linewidths=lw,colors=colors,alpha=alpha)
            # Set labels for x,y
            ax1.set_xlabel(xlab)
            ax1.set_ylabel(ylab)
            # Set xmin, xmax
            ax1.set_xlim(mnx,mxx)
            ax1.set_ylim(mnx,mxx)
            # Set aspect ratio
            ax1.set_aspect('equal')
            # Set colorbar
            cbar = fig.colorbar(im,label=clab, cax=cax)
            return   

        ani = animation.FuncAnimation(fig, animate, range(len(myfrms)),
                                      repeat=False)
    
    # Handle plotting a panel of subplots
    elif pflag:
        # Get the factors 
        fact    = factors(len(ifrm))
        # Set the number of panels in x and y direction 
        if len(fact) == 2:
            nxp = np.max(fact) 
            nyp = np.min(fact)
        else:
            nxp = np.max(fact[1:-1])
            nyp = np.min(fact[1:-1])
        
        # Set qminmax 
        if qflag:
            qmin, qmax = qminmax[0], qminmax[1]
        else:
            qmin, qmax = np.min(imgs), np.max(imgs) 
        
        # define ratio
        ratio = float(nxp)/float(nyp) 
        fsize = (ratio*7., 7.)
        # Create figure object
        fig, axes = plt.subplots(nyp,nxp, sharex='col', sharey='row',facecolor='white',
                                 figsize=(ratio*7.,7.)) 
        fig.subplots_adjust(hspace=0, wspace=0)
        # define the colorbar
        fig.subplots_adjust(top=0.8) 
        cax = fig.add_axes([0.16,0.85,0.7,0.02])


        # Now plot
        for (ax, iff) in zip(axes.flat, range(len(ifrm))):
            im = ax.pcolorfast(x1cf,x2cf,imgs[iff],vmin=qmin,vmax=qmax) 
            if ctrace:
                if tarr[iff] > params[0].thvcs:  
                    imc  = ax.contour(x1cc,x2cc,imgc[iff],levels=get_levels(imgc[iff],pcts)
                                                         ,colors=colors,alpha=alpha,linewidths=lw)
            ax.set_xlim(mnx,mxx)
            ax.set_ylim(mnx,mxx)
            #ax.set_aspect('equal')  
            ax.text(0.9*mnx, 0.8*mxx, 't = %1.1f [Myr]' % (tarr[iff]),
                         bbox={'facecolor':'white', 'alpha':0.9, 'pad':5})
            # This makes .eps files manageable 
            im.set_rasterized(True) 

        # define global labels
        fig.text(0.5, 0.04, xlab, ha='center')
        fig.text(0.06, 0.5, ylab, va='center', rotation='vertical') 

        # Set colorbar
        cb = fig.colorbar(im,cax=cax, orientation='horizontal') 
        cax.xaxis.set_ticks_position('bottom') 
        cb.ax.set_title(clab) 

    elif grid:
        imgs = np.zeros( imgs[0].shape ) 
        im = ax1.pcolor(x1cf,x2cf, imgs, facecolor='none',edgecolor='k')


        
   
    # Handle plotting a single frame 
    else: 
        # Set qminmax 
        if qflag:
            qmin, qmax = qminmax[0], qminmax[1]
        else:
            qmin, qmax = np.min(imgs), np.max(imgs) 

        ax1.set_title('t = %1.1f [Myr]' %(tarr[0]) )
        im   = ax1.pcolorfast(x1cf,x2cf, imgs[0], vmin=qmin,vmax=qmax)
        if ctrace:
            if tarr[0] > params[0].thvcs:
                imc  = ax1.contour(x1cc,x2cc,imgc[0],levels=get_levels(imgc[0],pcts)
                                                    ,colors=colors,alpha=alpha,linewidths=lw)
        # Same deal here 
        im.set_rasterized(True) 
        cbar = fig.colorbar(im,label=clab,cax=cax) 
    
        
    if save:
        mydir  = '/srv/analysis/jdupuy26/figures/'
        #mydir = os.getcwd()+'/'
        # Create file name (have to make sure it is unique for each sim to avoid overwrites)  
        #myname = os.path.basename(os.path.dirname(os.path.realpath('bgsbu.log')))
        myname = os.getcwd().split('longevity_study/',1)[1].replace('/','_') 
        if anim:
            print("[main]: Saving animation...")
            ani.save(mydir+myname+'_'+base+'_'+quant+'.mp4',fps=7.,
                     writer='imagemagick')

        elif pflag:
            print('[main]: Saving panel plot...')

            plt.savefig(mydir+myname+'_'+base+'_'+quant+'.'+fmt, dpi=80, format=fmt,bbox_inches='tight',writer='imagemagick')
        else:
            print("[main]: Saving frame...")
            plt.savefig(mydir+myname+'_'+base+'_'+quant+str(ifrm)+'.'+fmt, format=fmt,bbox_inches='tight')
    else:
        plt.show() 
    
main() 
print('Program complete.') 



