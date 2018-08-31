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
from scipy.interpolate import interp2d 

import argparse
from argparse import RawTextHelpFormatter
import time as t
import re

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
elif comp == 'xueboyudeair.wireless-1x.unc.edu':
    sys.path.insert(0,'/Users/xueboyu/desktop/plt')

else: 
    print('[init]: Computer %s not recognized!' % comp)
from units_class import * 
from read_bin import read_bin 
import read_athinput



#===================================================================
#
#  Code: gen_plot.py
#
#  Purpose: Reads in athena binary dumps and plots
#           quantities. Designed to be a general plotter
#           that plots all quantities, and any coordinate system,
#           including collisionless variables. Very similar to 
#           plot_sims.py 
#
#  Keywords: python gen_plot.py -h   
#
#  Usage: python plot_sims.py quant   
#
#  WARNING: THIS MUST BE RUN FROM SIMULATION DIRECTORY 
#          
#  Author: John Dupuy 
#          UNC Chapel Hill
#  Date:    05/14/18 
#  Updated: 05/14/18 
#====================================================================

#============ FUNCTIONS ==============================

#\func natural_sort()
# does natural (i.e. human) sorting of a list 
# of strings
# cf. https://stackoverflow.com/questions/4836710/
#            does-python-have-a-built-in-function-for-string-natural-sort
def natural_sort(l):
    convert      = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)  

#------------ FUNCTIONS FOR READING DATA -------------#
#\func get_files()
# returns dictionary of the form {'id0': bin_files0, 'id1': bin_files1,...}
# If not using mpi, treats CWD as 'id0' 
def get_files(mydir='./'):
    # Get processor directories 
    proc  = os.walk(mydir).__next__()[1]
    # Get rid of unwanted directories 
    proc  = sorted([x for x in proc if 'id' in x])
   
    # If not using mpi, proc will be empty, so append dummy
    #   directory to it  
    if not proc:
        proc.append('.') 
        mydir = '.'
    else:
        mydir = 'id0' 

    nproc = len(proc)
    # Define processor range
    p_no  = range(nproc)


    # Now create the dictionary 
    pdict = {} 

    #   id0 directory is special since it contains other stuff
    pdict[mydir] = sorted([fname for fname in
                           os.listdir(mydir+'/')
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

# \func get_log():
# takes the log of data 
def get_log(data):
    data[np.where(data <= 0.0)] = 1e-80
    data = np.log10(data)
    return data 

#\func get_athinput():
# reads athinput and returns base, params
def get_athinput(cwd=-1):
    if cwd == -1:
        cwd = os.getcwd()
    athin = cwd + '/'+[fnm for fnm in os.listdir(os.getcwd()+'/')
                                        if fnm.startswith('athinput')][0]
    base, params = read_athinput.readath(athin) 
    return base, params 

#\func get_quant 
# gets quantity 
def get_quant(file,quant,units,precision=32):
    #read in binary file
    nx1,nx2,nx3,x1,x2,x3,\
    d,M1,M2,M3,e,ie,s,\
    bx,by,bz,phi, \
    gamm1,cs,t,dt,nscalars,\
    clesshd, coordsys = read_bin(file,precision)
    
    # Set units
    if   units == 1:
        u = units_CGS()
    elif units == 2:
        u = units_SI()
    elif units == 0:
        u = units_COMP()

    ucgs  = units_CGS() 
    ucomp = units_COMP()

    T = 0
    if quant == 'T':
        T      = gamm1*ie*ucgs.esdens/(ucgs.k_b * ( d*ucgs.rhos/ucgs.m_h)) 

    # Parse collisionless variables  
    dcl, M1cl, M2cl, M3cl, E11, E22, E33, E12, E13, E23, p11, p22, p33 = clesshd 
    # Test whether CLESSHD is used 
    cless = dcl.any() 

    
    # Define dictionary for quants
    hydro_quants = {'E':e, 'ie':ie, 'd':d, 'mcent':d,
                  'n':d*ucgs.rhos/(ucgs.m_h), 'pie':gamm1*ie,
                  'p':gamm1*(e - 0.5*(M1**2.+M2**2.+M3**2.)/d),
                  'T':T,
                  'M':np.sqrt(M1**2.+M2**2.+M3**2.),
                  'v1':M1/d,'v2':M2/d,'v3':M3/d,
                  'M1':M1,'M2':M2,'M3':M3,
                  'v':np.sqrt(M1**2.+M2**2.+M3**2.)/d,
                  'cs': np.sqrt(gamm1*ie/d),
                  'phi':phi}
    if cless:
        cless_quants = {'dcl':dcl, 'M1cl':M1cl, 'M2cl':M2cl, 'M3cl':M3cl,
                      'Mcl':np.sqrt(M1cl**2.+M2cl**2.+M3cl**2.), 
                      'v1cl':M1cl/dcl,'v2cl':M2cl/dcl,'v3cl':M3cl/dcl,
                      'vcl':np.sqrt(M1cl**2.+M2cl**2.+M3cl**2.)/dcl,
                      'E11':E11, 'E22':E22, 'E33':E33, 'E12':E12, 'E13':E13, 'E23':E23, 
                      'P11':E11-(M1cl)*(M1cl/dcl),'P22':E22-(M2cl)*(M2cl/dcl),
                      'P33':E33-(M3cl)*(M3cl/dcl),'P12':E12-(M1cl)*(M2cl/dcl),
                      'P23':E23-(M2cl)*(M3cl/dcl),'P13':E13-(M1cl)*(M3cl/dcl),
                      'p11':p11, 'p22': p22, 'p33': p33}
    
    if cless:
        all_quants = dict(hydro_quants, **cless_quants)
    else:
        all_quants = hydro_quants 

    if quant == 'detP' or quant == 'detE':
        all_quants['detP'] = ( all_quants['P11']*all_quants['P22'] -
                               all_quants['P12']*all_quants['P12'] )
        all_quants['detE'] = ( all_quants['E11']*all_quants['E22'] -
                               all_quants['E12']*all_quants['E12'] )
    if quant == 'normP' or quant == 'normE':
        all_quants['normP'] = np.sqrt( all_quants['P11']**2. + 2.*all_quants['P12']**2. +
                                       all_quants['P22']**2. + 2.*all_quants['P13']**2. +
                                       all_quants['P33']**2. + 2.*all_quants['P23']**2. ) 
        all_quants['normE'] = np.sqrt( all_quants['E11']**2. + 2.*all_quants['E12']**2. +
                                       all_quants['E22']**2. + 2.*all_quants['E13']**2. +
                                       all_quants['E33']**2. + 2.*all_quants['E23']**2. )
    if quant == 'coordsys':
        return coordsys 

    return t, x1, x2, x3, all_quants[quant] 

#\func get_stitch
# given pdict, this will stitch together an image
# Here pdict is the dictionary from get_files()
def get_stitch(pdict,quant,myfrms,**kwargs):
    mydir    = './'
    vvec     = False
    for key in kwargs:
        if key == 'mydir':
            mydir = kwargs[key]
        if key == 'vvecs':
            vvec = kwargs[key]
    
    # read athinput
    base, params = get_athinput() 
   
    # get no. of points 
    nx1 = int(params[0].nx1)
    nx2 = int(params[0].nx2)
    nx3 = int(params[0].nx3) 

    if nx3 == 0:
        nx3 = 1
    if nx2 == 0:
        nx2 = 1 

    # get extent of grid
    mn1 = params[0].x1min
    mx1 = params[0].x1max
    mn2 = params[0].x2min
    mx2 = params[0].x2max
    mn3 = params[0].x3min
    mx3 = params[0].x3max 

    # get no. of processors
    nproc = len(pdict) 
    # get no. of processors in x1 dir
    if nproc > 1:
        ngridx1 = int(params[0].ngridx1)
        ngridx2 = int(params[0].ngridx2) 
        ngridx3 = int(params[0].ngridx3) 
        if ngridx2 == 0:
            ngridx2 = 1
        if ngridx3 == 0:
            ngridx3 = 1
    else: 
        ngridx1 = 1
        ngridx2 = 1
        ngridx3 = 1

    # get size of grid on each processor
    npx1    = nx1/ngridx1
    npx2    = nx2/ngridx2 
    npx3    = nx3/ngridx3 

    # get no. of files in each processor
    nf     = len(myfrms) 

    # Make array to hold quant for all time
    quant_arr = np.zeros((nf, nx3, nx2, nx1))
    if vvec:
        v1_arr = np.zeros((nf, nx3, nx2, nx1))
        v2_arr = np.zeros((nf, nx3, nx2, nx1))
    # Make time, x1, x2 arrays 
    tarr = np.zeros(nf)
    x1   = np.zeros(nx1)
    x2   = np.zeros(nx2) 
    x3   = np.zeros(nx3) 

    keys = natural_sort(list(pdict.keys())) 

    # Define processor index
    ip = 0
    for ip3 in range(ngridx3): 
        fac3    = ip3%ngridx3 
        for ip2 in range(ngridx2):
            fac2     = ip2%ngridx2 
            for ip1 in range(ngridx1): 
                fac1     = ip1%ngridx1
                # Define processor key
                pkey     = keys[ip]  
                i        = 0            # index for filling arrays
                
                for iff in myfrms: 
                    # Get filename
                    fnm      = pkey + '/' + pdict[pkey][iff]
                    # Fill data array
                    data = get_quant(fnm,quant,0) 
                    # parse data 
                    tarr[i]                                = data[0]
                    x1[int(fac1*npx1):int(npx1*(1+fac1))]            = data[1]
                    x2[int(fac2*npx2):int(npx2*(1+fac2))]            = data[2]
                    x3[int(fac3*npx3):int(npx3*(1+fac3))]            = data[3]
                    quant_arr[i,  int(fac3*npx3):int(npx3*(1+fac3)),
                                  int(fac2*npx2):int(npx2*(1+fac2)),
                                  int(fac1*npx1):int(npx1*(1+fac1))] = data[4]
                    
                    # Get velocity data 
                    if vvec:
                        data = get_quant(fnm,'v1',units,params)
                        # Get vr
                        v1_arr[i, fac2*npx2:npx2*(1+fac2),
                                  fac1*npx1:npx1*(1+fac1)] = data[3]
                        data = get_quant(fnm,'v2',units,params)
                        # Get v\phi
                        v2_arr[i, fac2*npx2:npx2*(1+fac2),
                                  fac1*npx1:npx1*(1+fac1)] = data[3] 
                    i += 1
                ip += 1

    if vvec:
        # Convert velocities to vx, vy
        v1_arr, v2_arr = vel_p2c(v1_arr, v2_arr, x1, x2)
        return tarr, x1, x2, quant_arr, v1_arr, v2_arr 
    else:
        return tarr, x1, x2, x3, quant_arr

# \func get_fc()
# given cell centered x1, x2, returns face centered grid 
# currently only necessary for 2d 
def get_fc(params, lg=False, apc=False):
    
    # if using log grid (ilog=1)  
    if lg:
        # x1 direction  
        lx1f = np.linspace(np.log(params[0].x1min),
                           np.log(params[0].x1max),params[0].nx1+1)
        # x2 direction
        x2f  = np.linspace(       params[0].x2min,
                                  params[0].x2max, params[0].nx2+1) 
        x1f  = np.exp(lx1f) 
    else: 
        # x1 direction  
        x1f =  np.linspace(params[0].x1min,
                           params[0].x1max,params[0].nx1+1)
        # x2 direction
        x2f  = np.linspace(params[0].x2min,
                           params[0].x2max, params[0].nx2+1) 

    if apc:      
        # get area per cell 
        x1c = 0.5*(x1f[1:] + x1f[0:-1])*1.e3
        ar  =     (x1f[1:] - x1f[0:-1])*1.e3  # convert to pc 
        ap  = x1c*(x2f[1:] - x2f[0:-1])

        # Area per cell
        return np.meshgrid(ar*ap, ar*ap, indexing='xy') 
    else:
        return np.meshgrid(x1f, x2f, indexing='xy')   



def get_labels(quant,dim,log=False):

    # Define dictionary for quants
    lab_quants = {  # Fluid variables  
                  'E':'Energy density', 'ie':'Internal energy density', 
                  'd':'$\\rho$',
                  'n':'Column density', 'p':'Surface pressure',
                  'pie':'Surface pressure (from U.IE)', 
                  'T':'T',
                  'M':'M$_{tot}$',
                  'v1':'v$_1$','v2':'v$_2$','v3':'v$_3$',
                  'M1':'M$_1$','M2':'M$_2$','M3':'M$_3$',
                  'phi':'$\Phi_G$',
                  'cs':'c$_s$', 
                  'v':'v$_{tot}$',
                  's1':'$\Sigma_c$','s1c':'$M_c (R < 0.5 \; {\\rm [kpc])}/M_c$',
                  'jl':'$\lambda_J$','vlos':'$v_{\\rm los}$',
                    # Collisionless variables 
                  'dcl':'$\\rho_{\\rm cl}$', 
                  'v1cl':'v$_{1,\\rm cl}$','v2cl':'v$_{2,\\rm cl}$',
                  'v3cl':'v$_{3,\\rm cl}$','M1cl':'M$_{1,\\rm cl}$',
                  'vcl':'v$_{\\rm cl}$','Mcl':'M$_{\\rm cl}$', 
                  'M2cl':'M$_{2,\\rm cl}$','M3cl':'M$_{3,\\rm cl}$',
                  'P11':'P$_{11}$','P22':'P$_{22}$','P33':'P$_{33}$',
                  'p11':'P$_{ie,11}$','p22':'P$_{ie,22}$','p33':'P$_{ie,33}$',
                  'P12':'P$_{12}$','P23':'P$_{23}$','P13':'P$_{13}$',
                  'E11':'E$_{11}$','E22':'E$_{22}$','E33':'P$_{33}$',
                  'E12':'E$_{12}$','E23':'E$_{23}$','E13':'P$_{13}$',
                  'detP':'det(P$_{ij}$', 'detE':'det(E$_{ij}$',
                  'normP':'$|| P_{ij} ||_F$', 'normE':'$|| E_{ij} ||_F$'}

    # Define dictionary for units 
    units = ' [comp units]' 

    if dim == 1:
        xlabel = 'x'+units
        ylabel = ''
    if dim == 2 or dim == 3:
        xlabel = 'x'+units
        ylabel = 'y'+units

    cbar_l = lab_quants[quant]
    cbar_l += units 

    

    if log:
        cbar_l = 'log$_{10}$('+cbar_l+')'

    return xlabel, ylabel, cbar_l 


#\func get_args()
# this function parses CMD line args
def get_args():
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
                        required=False,default=[-np.pi,np.pi],
                        help="Plotting range in x and y\n"
                             "Note: assumes x, y range are equivalent")
    parser.add_argument("--save", dest="save",action='store_true',
                        default=False,
                        help="Switch to save anim or figure")
    parser.add_argument("--log",dest="log",action='store_true',
                        default=False,
                        help="Switch to take log of images, default is False")
    parser.add_argument("--units",dest="units",type=int,required=False,
                        default=0, help="units: 0-comp,1-cgs,2-SI")  
    parser.add_argument("--grid", dest="grid",action='store_true',
                        default=False, help="Switch to make plot to show grid")
    parser.add_argument("--fmt", dest="fmt", default='eps',
                        type=str, help='format for saving graphics, default: eps') 
    parser.add_argument("--vvec",dest="vvec",action='store_true',
                        default=False, required=False,
                        help="Overplot velocity vectors\n")
    parser.add_argument("--noplot",dest="noplot",action='store_true',
                        default=False, required=False,
                        help="Switch to return only stitched together array\n"
                             "To be used if this file is imported from another file\n")
    parser.add_argument("--sliced1d",dest="sliced1d",action='store_true',
                        default=False, required=False,
                        help="Switch to take 1D slice of 2D array along DIAGONAL\n")
    parser.add_argument("--slicel1d",dest="slicel1d",action='store_true',
                        default=False, required=False,
                        help="Switch to take 1D slice of 2D array along a line\n")
    parser.add_argument("--col",dest="col",action='store_true',
                        default=False, required=False,
                        help="Sum 3d simuations along the z-axis, so they can be viewed as" 
                             "2d plots\n") 
    parser.add_argument("--slice",dest="slc",type=int, required=False,
                        default='-1', help="Slice 3d array into 2d along INT axis\n")
    parser.add_argument("--lg",dest="lg",action='store_true',
                        default=False, required=False, 
                        help="Switch for turning on log-grid")
    parser.add_argument("--nocyl",dest="nocyl",action='store_true',
                        default=False, required=False,
                        help="Switch for plotting cylindrical simulations as R, \phi") 
    return parser.parse_args() 
    
 

#------------- MAIN FUNCTION ------------------------#
def main(args):
    ctable = 'magma'
    plt.rcParams['image.cmap'] = ctable
    base, params = get_athinput() 
    mc           = params[0].mhvc

    # parsing arguments            
    quant = args.quant
    anim  = args.anim
    iani  = args.iani
    ifrm  = args.ifrm
    save  = args.save
    log   = args.log
    mnx, mxx = args.mnmx
    iunit    = args.units 
    qminmax  = args.qminmax
    grid     = args.grid
    fmt      = args.fmt
    vvec     = args.vvec
    noplot   = args.noplot
    sliced1d = args.sliced1d
    slicel1d = args.slicel1d
    col      = args.col
    slc      = args.slc   
    lg       = args.lg 
    nocyl    = args.nocyl

    # get files 
    pdict = get_files()
    # get the coordinate system 
    mykey = list(pdict.keys())[0] 
    fnm   = mykey + '/' + pdict[mykey][0]
    coordsys = get_quant(fnm, 'coordsys', iunit)

    # Get qminmax flag 
    qflag = True if np.size(qminmax) > 1 else False
    # Get panel flag
    pflag = True if np.size(ifrm) > 1 else False 
    # Get mnmx flag
    mnmxflag = False if mxx == np.pi else True  
    # Get slice flag
    slice2d  = False if slc == -1 else True 
    # Set the coordinate system flag
    cyl      = True if coordsys == -2 and not nocyl else False 
    
    if np.size(ifrm) == 1: ifrm = ifrm[0]
    
    # Change default iani values
    if iani[1] == 0:
        iani[1] = len(pdict[list(pdict.keys())[0]]) 

    # determine myframes
    if anim:
        myfrms = range(iani[0],iani[1])
    elif pflag:
        myfrms = ifrm
    else: 
        myfrms = [ifrm] 
    # get data
    if vvec:
        tarr, x1, x2, imgs, vx_imgs, vy_imgs = get_stitch(pdict,quant,myfrms,vvecs=vvec)
    else:
        tarr, x1, x2, x3, imgs = get_stitch(pdict,quant,myfrms)
    if log:
        imgs = get_log(imgs) 


    # Get dimensions of data
    nf, nx3, nx2, nx1 = imgs.shape 
    # Determine dimensional plotting flag
    flag3d, flag2d, flag1d = False, False, False  
    if nx3 > 1:
        flag3d = True
        dim    = 3
    elif (nx3 == 1) and (nx2 > 1):
        flag2d = True
        dim    = 2
    else:
        flag1d = True 
        dim    = 1

    # Change flags for slicing 
    if sliced1d or slicel1d:
        if flag2d:
            flag2d = False
            flag1d = True 
        if flag3d:
            flag3d = False
            flag1d = True 
    if col:
        flag3d = False
        flag2d = True
    if slice2d:
        flag3d = False
        flag2d = True 
    
    # Get face-centered grid if cylindrical  
    if flag2d and cyl: 
        x1f , x2f  = get_fc(params, lg) 
        x1cf, x2cf = x1f*np.cos(x2f), x1f*np.sin(x2f) 

    # Determine labels 
    xlab, ylab, clab = get_labels(quant,dim,log)

    # Now plot the data 
    if flag3d:
        from mpl_toolkits.mplot3d import Axes3D 
        fig = plt.figure()
        ax1 = fig.add_subplot(111, projection='3d')  
    else:
        fig = plt.figure(figsize=(7.5,5.5),facecolor='white') 
        ax1 = fig.add_subplot(111) 

    if flag1d: 
        if sliced1d:
            # Note this assumes nx1 = nx2 = nx3 
            if dim == 2:
                imgs = imgs[:,0,:,:]
                # Take slice 
                vals = []
                for j in range(len(imgs)):
                    vals.append(np.array([imgs[j,i,i] for i in range(len(x1))]))
                imgs = vals 
                # Get distance along diagonal 
                rad   = np.sqrt(x1**2. + x2**2.) 
                rad[x1 < 0] = -rad[x1 < 0] 
                x1    = rad.copy() 
            if dim == 3:
                # Take slice
                vals = []
                for j in range(len(imgs)):
                    vals.append(np.array([imgs[j,i,i,i] for i in range(len(x1))]))
                imgs = vals 
                # Get distance along diagonal
                rad = np.sqrt(x1**2. + x2**2. + x3**2.)
                rad[x1 < 0] = -rad[x1 < 0]
                x1 = rad.copy() 

        elif slicel1d:
            if dim == 2:
                imgs = imgs[:,0,:,:]
                # Take slice 
                vals = []
                for j in range(len(imgs)):
                    vals.append(np.array([imgs[j,nx2/2,i] for i in range(len(x1))]))
                imgs = vals 
            if dim == 3:
                # Take slice
                vals = []
                for j in range(len(imgs)):
                    vals.append(np.array([imgs[j,nx3/2,nx2/2,i] for i in range(len(x1))]))
                imgs = vals 
            
        # Get rid of unnecessary dimensions 
        else:
            imgs = imgs[:,0,0,:] 

        # Handle animation
        if anim:
            if qflag:
                qmin, qmax = qminmax[0], qminmax[1]
            else:
                qmin, qmax = np.min(imgs), np.max(imgs) 
            
            # Set labels 
            ax1.set_title('t = %1.2f' % (tarr[0]) ) 
            ax1.set_xlabel(xlab)
            ax1.set_ylabel(clab) 
            # Plot first frame 
            ax1.plot(x1, imgs[0], '.')
            
            def animate(ifrm):
                # Clear figure
                ax1.cla()
                # Set title & labels 
                ax1.set_title('t = %1.2f' % (tarr[ifrm]) )
                ax1.set_xlabel(xlab)
                ax1.set_ylabel(clab) 
                # Set xmin, xmax 
                if mnmxflag:
                    ax1.set_xlim(mnx, mxx)
                #ax1.set_ylim(qmin,qmax)
                # Now plot
                ax1.plot(x1, imgs[ifrm],'.')
                return
            ani = animation.FuncAnimation(fig, animate, range(len(myfrms)),
                                          repeat=False) 

        # Handle plotting a single frame 
        else:
            # plot 
            ax1.plot(x1, imgs[0],'.')  
            ax1.set_xlabel(xlab)
            ax1.set_ylabel(clab)
            ax1.set_title('t = %1.2f' % (tarr[0]) ) 
    
    elif flag2d:
        # Get extent of grid 
        if cyl:
            mnx1, mxx1, mnx2, mxx2 = ( np.min(x1cf), np.max(x1cf),
                                       np.min(x2cf), np.max(x2cf) )
        else: 
            mnx1, mxx1, mnx2, mxx2 = ( np.min(x1), np.max(x1),
                                       np.min(x2), np.max(x2) ) 
            dx1, dx2 = x1[1]-x1[0], x2[1]-x2[0]

            mnx1 -= 0.5*dx1; mxx1 += 0.5*dx1
            mnx2 -= 0.5*dx2; mxx2 += 0.5*dx2
        # Determine colorbar 
        div = make_axes_locatable(ax1)
        cax = div.append_axes('right', '5%', '5%') 


        # Get rid of unnecessary dimensions 
        if col:
            imgs = np.sum(imgs, axis=1) 
        elif slice2d:
            imgs = imgs[:,slc,:,:] 
        else:
            imgs = imgs[:,0,:,:] 

        if qflag:
            qmin, qmax = qminmax[0], qminmax[1]
        else:
            qmin, qmax = np.min(imgs), np.max(imgs) 
    
        # Handle animation
        if anim: 
            ax1.set_title('t = %1.2f' %(tarr[0]) )
            if cyl:
                im   = ax1.pcolorfast(x1cf,x2cf, imgs[0], vmin=qmin,vmax=qmax)
                ax1.set_aspect('equal')
            else: 
                im   = ax1.imshow(imgs[0], extent=(mnx1, mxx1, mnx2, mxx2),
                                           vmin=qmin,vmax=qmax, origin='lower',
                                           interpolation='None')
            im.set_rasterized(True) 
            cbar = fig.colorbar(im,label=clab,cax=cax) 

            def animate(ifrm):
                # Clear figure 
                ax1.cla()
                cax.cla()
                
                # Set title 
                ax1.set_title('t = %1.2f' %(tarr[ifrm]) )
                # Plot 
                if cyl:
                    im   = ax1.pcolorfast(x1cf,x2cf, imgs[ifrm], vmin=qmin,vmax=qmax)
                else:
                    im   = ax1.imshow(imgs[ifrm], extent=(mnx1, mxx1, mnx2, mxx2),
                                                  vmin=qmin,vmax=qmax, origin='lower',
                                                  interpolation='None')
                im.set_rasterized(True) 
                # Set labels for x,y
                ax1.set_xlabel(xlab)
                ax1.set_ylabel(ylab)
                # Set xmin, xmax
                if mnmxflag:
                    ax1.set_xlim(mnx,mxx)
                    ax1.set_ylim(mnx,mxx)
                else:
                    ax1.set_xlim(mnx1,mxx1)
                    ax1.set_ylim(mnx2,mxx2) 

                # Set aspect ratio
                ax1.set_aspect('equal')
                # Set colorbar
                cbar = fig.colorbar(im,label=clab, cax=cax)
                return   

            ani = animation.FuncAnimation(fig, animate, range(len(myfrms)),
                                          repeat=False)
        # Handle a single frame 
        else:
            print(np.mean(imgs[0]))  
            ax1.set_xlabel(xlab)
            ax1.set_ylabel(ylab)

            ax1.set_title('t = %1.2f' %(tarr[0]) )
            if cyl:
                im   = ax1.pcolorfast(x1cf,x2cf, imgs[0], vmin=qmin,vmax=qmax)
                ax1.set_aspect('equal')
            else:
                im   = ax1.imshow(imgs[0], extent=(mnx1, mxx1, mnx2, mxx2),
                                           vmin=qmin,vmax=qmax, origin='lower',
                                           interpolation='None')
            im.set_rasterized(True) 
            cbar = fig.colorbar(im,label=clab,cax=cax) 

    # Handle 3d plotting (using mayavi) 
    elif flag3d: 
        # Get extent of grid 
        mnx1, mxx1, mnx2, mxx2, mnx3, mxx3 = ( np.min(x1), np.max(x1),
                                               np.min(x2), np.max(x2),
                                               np.min(x3), np.max(x3) ) 
        dx1, dx2, dx3 = x1[1]-x1[0], x2[1]-x2[0], x3[1]-x3[0]

        mnx1 -= 0.5*dx1; mxx1 += 0.5*dx1
        mnx2 -= 0.5*dx2; mxx2 += 0.5*dx2
        mnx3 -= 0.5*dx3; mxx3 += 0.5*dx3

        x1, x2, x3 = np.meshgrid(x1,x2,x3) 
        
        ax1.scatter(x1,x2,x3, c=imgs[0].ravel(),cmap=plt.hot()) 
    else: 
        print("[main]: Unsure what to plot  :( ")  
            

    if save:
        #mydir  = '/srv/analysis/jdupuy26/figures/'
        mydir = os.getcwd()+'/'
        # Create file name (have to make sure it is unique for each sim to avoid overwrites)  
        #myname = os.path.basename(os.path.dirname(os.path.realpath('bgsbu.log')))
        #myname = os.getcwd().split('longevity_study/',1)[1].replace('/','_') 
        myname = '' 
        if anim:
            print("[main]: Saving animation...")
            ani.save(mydir+myname+'_'+base+'_'+quant+'.gif',fps=15.
                     ,writer='imagemagick')

        else:
            print("[main]: Saving frame...")
            plt.savefig(mydir+myname+'_'+base+'_'+quant+str(ifrm)+'.'+fmt, format=fmt,bbox_inches='tight')
    else:
        plt.show() 

        
if __name__ == '__main__':
   # If this file is called from cmd line 
   args = get_args() 
   main(args) 

