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

#------------ FUNCTIONS FOR READING DATA -------------#
#\func get_files()
# returns dictionary of the form {'id0': bin_files0, 'id1': bin_files1,...}
# If not using mpi, treats CWD as 'id0' 
def get_files(mydir='./'):
    # Get processor directories 
    proc  = os.walk(mydir).next()[1]
    # Get rid of unwanted directories 
    proc  = [x for x in proc if 'id' in x]
   
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
def get_quant(file,quant,units,precision=64):
    #read in binary file
    nx1,nx2,nx3,x1,x2,x3,\
    d,M1,M2,M3,e,ie,s,\
    bx,by,bz,phi, \
    gamm1,cs,t,dt,nscalars,\
    clesshd                  = read_bin(file,precision)
    
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
    dcl, M1cl, M2cl, M3cl, E11, E22, E33, E12, E13, E23 = clesshd 

    
    # Define dictionary for quants
    all_quants = {'E':e, 'ie':ie, 'd':d, 'mcent':d,
                  'n':d*ucgs.rhos/(ucgs.m_h), 'pie':gamm1*ie,
                  'p':gamm1*(e - 0.5*(M1**2.+M2**2.+M3**2.)/d),
                  'T':T,
                  'M':np.sqrt(M1**2.+M2**2.+M3**2.),
                  'v1':M1/d,'v2':M2/d,'v3':M3/d,
                  'M1':M1,'M2':M2,'M3':M3,
                  'V':np.sqrt(M1**2.+M2**2.+M3**2.)/d,
                  'cs': np.sqrt(gamm1*ie/d),
                  'dcl':dcl, 'M1cl':M1cl, 'M2cl':M2cl, 'M3cl':M3cl,
                  'v1cl':M1cl/dcl,'v2cl':M2cl/dcl,'v3cl':M3cl/dcl,
                  'E11':E11, 'E22':E22, 'E33':E33, 'E12':E12, 'E13':E13, 'E23':E23, 
                  'P11':E11-(M1cl)*(M1cl/dcl),'P22':E22-(M2cl)*(M2cl/dcl),
                  'P33':E33-(M3cl)*(M3cl/dcl),'P12':E12-(M1cl)*(M2cl/dcl),
                  'P23':E23-(M2cl)*(M3cl/dcl),'P13':E13-(M1cl)*(M3cl/dcl)} 

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

    keys = list(pdict.keys()) 

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
                    x1[fac1*npx1:npx1*(1+fac1)]            = data[1]
                    x2[fac2*npx2:npx2*(1+fac2)]            = data[2]
                    x3[fac3*npx3:npx3*(1+fac3)]            = data[3]
                    quant_arr[i,  fac3*npx3:npx3*(1+fac3),
                                  fac2*npx2:npx2*(1+fac2),
                                  fac1*npx1:npx1*(1+fac1)] = data[4]
                    
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
                  'cs':'c$_s$', 
                  'V':'V$_{tot}$',
                  's1':'$\Sigma_c$','s1c':'$M_c (R < 0.5 \; {\\rm [kpc])}/M_c$',
                  'jl':'$\lambda_J$','vlos':'$v_{\\rm los}$',
                    # Collisionless variables 
                  'dcl':'$\\rho_{\\rm cl}$', 
                  'v1cl':'v$_{1,\\rm cl}$','v2cl':'v$_{2,\\rm cl}$',
                  'v3cl':'v$_{3,\\rm cl}$','M1cl':'M$_{1,\\rm cl}$',
                  'M2cl':'M$_{2,\\rm cl}$','M3cl':'M$_{3,\\rm cl}$',
                  'P11':'P$_{11}$','P22':'P$_{22}$','P33':'P$_{33}$',
                  'P12':'P$_{12}$','P23':'P$_{23}$','P13':'P$_{13}$',
                  'E11':'E$_{11}$','E22':'E$_{22}$','E33':'P$_{33}$',
                  'E12':'E$_{12}$','E23':'E$_{23}$','E13':'P$_{13}$',
                  }

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
                        required=False,default=[-15.0,15.0],
                        help="Plotting range in x and y\n"
                             "Note: assumes x, y range are equivalent")
    parser.add_argument("--save", dest="save",action='store_true',
                        default=False,
                        help="Switch to save anim or figure")
    parser.add_argument("--log",dest="log",action='store_true',
                        default=False,
                        help="Switch to take log images, default is False")
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

    # Get qminmax flag 
    qflag = True if np.size(qminmax) > 1 else False
    # Get panel flag
    pflag = True if np.size(ifrm) > 1 else False 
    
    if np.size(ifrm) == 1: ifrm = ifrm[0]
    
    # get files 
    pdict = get_files() 
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

    # Determine labels 
    xlab, ylab, clab = get_labels(quant,dim,log)
    

    # Now plot the data 
    fig = plt.figure(figsize=(7.5,5.5),facecolor='white') 
    ax1 = fig.add_subplot(111) 

    if flag1d: 
        # Get rid of unnecessary dimensions 
        imgs = imgs[:,0,0,:] 
        # plot 
        ax1.plot(x1, imgs[0],'.')  
        ax1.set_xlabel(xlab)
        ax1.set_ylabel(clab)
        ax1.set_title('t = %1.2f' % (tarr[0]) ) 


    plt.show() 
        
if __name__ == '__main__':
   # If this file is called from cmd line 
   args = get_args() 
   main(args) 

