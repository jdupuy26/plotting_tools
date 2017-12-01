#!/usr/bin/python
import numpy as np
import sys 
import os
import fnmatch as fnm
import argparse
from argparse import RawTextHelpFormatter

# matplotlib imports
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib.pyplot import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable 

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

from read_bin import read_lv as rlv
import read_athinput
from units_class import * 

#=====================================================
#
#  Code: plot_lv.py
#
#  Purpose: Plot files read in from lv dumps 
#
#  Keywords: python plot_lv.py -h  
#
#  Usage: python plot_lv.py quant (--ifrm, --anim, --iani) --iang  
#          
#  Author: John Dupuy 
#          UNC Chapel Hill
#  Date:    12/01/17
#  Updated: 12/01/17 
#=====================================================

#============FUNCTIONS===========
                        
#==============================
# get_files() function 
def get_files(dirs,proc,flag):
    # input: ALL TYPE str
    #       dirs: directory that contains athinput
    #       proc: processor no.
    #       flag: unique flag for your file 
    # get the files you want to do analysis on 
    return [dirs+proc+'/'+fname for 
            fname in os.listdir(dirs+proc+'/')
            if fnm.fnmatch(fname,flag)] 

#===========================
# get_data function
def get_data(file,**kwargs):
    # Set units 
    iunit = 0
    # Set precision
    prec  = 32
    for key in kwargs: 
        if key == 'iunit':
            iunit = kwargs[key] 
        if key == 'prec':
            prec  = kwargs[key]
    
    t, lvals, vvals,lvdiag = rlv(file,prec, **kwargs)

    if iunit == 1: 
        u = units_CGS()
    elif iunit == 2:
        u = units_SI()
    else: 
        u = units_COMP()

    # Do unit conversion
    t        *= u.myr

    return t, lvals, vvals, lvdiag

#================================
# get_title() function
def get_title(quant):
    lab = 'def'
    if quant == 'lv':
        lab = 'Intensity of HI emission'
    return lab

#================================
# get_xlabel() function
def get_xlabel(quant):
    lab = 'def'
    if quant == 'lv':
        lab = 'l [deg]'
    return lab

#================================
# get_ylabel() function
def get_ylabel(quant):
    lab = 'def'
    if quant == 'lv':
        lab = 'v [km/s]'
    return lab

#================================
# initialize() function
def init(quant, files, **kwargs):
    ctrl_files = 1
    # parse keyword args
    for key in kwargs:
        if key == 'ctrl_files':
            ctrl_files = kwargs[key]
            nctrl      = len(ctrl_files)

    # set length of simulation from file list
    n    = len(files)
    # setup time array
    tarr = np.zeros(n)

    # choose quant
    if quant == "lv":
        # Set up lists to store data
        lvdiags = []
        vvals   = [] 
        
        # Read in data
        for i in range(n):
            t, l, v, lv = get_data(files[i], **kwargs)
            
            tarr[i]     = t
            vvals.append(v)
            lvdiags.append(lv) 
        
        # lvals doesn't change w/ time 
        lvals   = l 

    
    return tarr, lvals, vvals, lvdiags 


#===============================
# main()
def main():
    # Get input file
    athdir = os.getcwd()+'/'
    athin = [fname for fname in os.listdir(athdir) if fname.startswith("athinput")][0]
    # Read input file
    base, params = read_athinput.readath(athdir+athin) 
    # Set precision
    prec  = 32

    # Get otf files 
    lv_files = get_files(athdir,'id0','*.lv.*')
    # Sort them
    lv_files.sort()
    
    n   = len(lv_files)

    # Read in system arguments
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter) 
    parser.add_argument("quant",type=str,
                        help="Plotting options:\n"
                             "  lv: HI emission f(l,v)\n")
    parser.add_argument("--anim", dest="anim",action='store_true',
                        default=False,
                        help="Switch to do animation\n")
    parser.add_argument("--iani", dest="iani",nargs=2,required=False,
                        default=(0,n-1),type=int,
                        help="Animate from frame iani[0] to iani[1]\n")
    parser.add_argument("--ifrm", dest="ifrm",type=int,default=None,
                        help="Frame of simulation to plot:\n"
                             "  0: tsim = 0 \n"
                             "  1: tsim = dt_dump\n"
                             "  2: tsim = 2*dt_dump\n"
                             "  .               \n"
                             "  .               \n"
                             "  .               \n"
                             "  n: tsim = n*dt_dump\n")
    parser.add_argument("--save", dest="save",action='store_true',
                        default=False,
                        help="Switch to save anim or figure")

    # parsing arguments            
    args  = parser.parse_args()
    quant = args.quant
    anim  = args.anim
    iani  = args.iani
    ifrm  = args.ifrm
    save  = args.save

    # Conflicting argument checks 
    if anim and ifrm != None:
        print("[main]: specifying a single frame while animating doesn't make sense")
        quit()
    if ifrm > (n-1):
        print("[main]: frame out of range of simulation!")
        quit()

    # Get the data
    tarr, lvals, vvals, lvdiags = init(quant,lv_files,prec=prec)

        # Get labels 
    title = get_title(quant)
    xlab  = get_xlabel(quant)
    ylab  = get_ylabel(quant)

    if anim:
        print("[main]: Animating from t = %1.1f [Myr]\n" 
              "                    to t = %1.1f [Myr]\n"
                    %( tarr[iani[0]], tarr[iani[1]] ) )  
    # Open figure
    fig = plt.figure(figsize=(7.0,5.5),facecolor='white')
    ax1 = fig.add_subplot(111)
    
    # Handle animation
    if anim:
        # Get spacing
        dl = lvals[1] - lvals[0]
        dv = vvals[iani[0]][1] - vvals[iani[1]][0]
        #  get extent
        mnl = lvals[ 0] - 0.5*dl
        mxl = lvals[-1] + 0.5*dl
        mnv = vvals[iani[0]][ 0] - 0.5*dv
        mxv = vvals[iani[1]][-1] + 0.5*dv

        im = ax1.imshow(lvdiags[iani[0]], extent=[mnl,mxl,mnv,mxv],
                        origin='lower',aspect='auto')
        ax1.set_xlim(mnl,mxl)
        ax1.set_ylim(mnv,mxv)
        ax1.set_xlabel(xlab)
        ax1.set_ylabel(ylab)
        ax1.set_title('t = %1.1f [Myr]' % tarr[0]) 

        div  = make_axes_locatable(ax1)
        cax  = div.append_axes('right', '5%', '5%')
        cbar = fig.colorbar(im,label=title, cax=cax)


        def animate(ifrm):
            dv  = vvals[ifrm][1 ] - vvals[ifrm][0]
            mnv = vvals[ifrm][ 0] - 0.5*dv
            mxv = vvals[ifrm][-1] + 0.5*dv
            # Update ylim
            ax1.set_ylim(mnv,mxv)
            # update title
            ax1.set_title('t = %1.1f [Myr]' % tarr[ifrm])
            # make the plot
            im = ax1.imshow(lvdiags[ifrm], extent=[mnl,mxl,mnv,mxv],
                        origin='lower',aspect='auto')
            cax.cla()
            fig.colorbar(im,label=title,cax=cax)

            return im 
        # do the animation
        ani = animation.FuncAnimation(fig, animate, range(iani[0],iani[1]+1), repeat=False)
        if save:
            print("[main]: Saving animation")
            ani.save(quant+".mp4")
    else:
        # Get spacing
        dl = lvals[1] - lvals[0]
        dv = vvals[ifrm][1] - vvals[ifrm][0]
        #  get extent
        mnl = lvals[ 0] - 0.5*dl
        mxl = lvals[-1] + 0.5*dl
        mnv = vvals[ifrm][ 0] - 0.5*dv
        mxv = vvals[ifrm][-1] + 0.5*dv

        im = ax1.imshow(lvdiags[ifrm], extent=[mnl,mxl,mnv,mxv],
                        origin='lower',aspect='auto')
        ax1.set_xlim(mnl,mxl)
        ax1.set_ylim(mnv,mxv)
        ax1.set_xlabel(xlab)
        ax1.set_ylabel(ylab)
        ax1.set_title('t = %1.1f [Myr]' % tarr[ifrm]) 

        plt.colorbar(im,label=title)

        if save:
            print("[main]: Saving figure")
            plt.savefig(quant+str(ifrm)+".eps")
    
    if save:
        print("[main]: Program complete")
    else:
        plt.show()
    #================================================================================#
            


        
main()

