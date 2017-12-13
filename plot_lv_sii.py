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
from scipy.ndimage.filters import gaussian_filter

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

from read_bin import * 
import read_athinput
from units_class import * 

# Set the colormap
plt.rcParams['image.cmap'] = 'magma'

#=====================================================
#
#  Code: plot_lv_sii.py
#
#  Purpose: Plot files read in from lv and sii dumps 
#
#  Keywords: python plot_lv_sii.py -h  
#
#  Usage: python plot_lv_sii.py quant (--ifrm, --anim, --iani)   
#          
#  Author: John Dupuy 
#          UNC Chapel Hill
#  Date:    12/01/17
#  Updated: 12/12/17 
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
    # Set quant
    quant = 'lv'
    iline = 0
    for key in kwargs: 
        if key == 'iunit':
            iunit = kwargs[key] 
        if key == 'prec':
            prec  = kwargs[key]
        if key == 'quant':
            quant = kwargs[key]
        if key == 'iline':
            iline = kwargs[key]
    if quant == 'lv':
        t, x, y, lv, lvc = read_lv(file,prec, **kwargs)
        data = lv
    elif quant == 'lvc':
        t, x, y, lv, lvc = read_lv(file,prec, **kwargs)
        data = lvc
    elif quant == 'sii':
        t, x, y, data = read_sii(file,prec,**kwargs)
        if iline == 2:
            data = data[0]/data[1] # line ratio
        else:
            data = data[iline]

    if iunit == 1: 
        u = units_CGS()
    elif iunit == 2:
        u = units_SI()
    else: 
        u = units_COMP()

    # Do unit conversion
    t *= u.myr

    return t, x, y, data 

#================================
# get_log function
def get_log(data):
    data[ data <= 0] = 1e-80
    data = np.log10(data)
    return data 

#================================
# get_title() function
def get_title(quant,nolog):
    lab = 'def'
    if quant == 'lv' or quant == 'lvc':
        lab = 'Intensity of HI emission'
    elif quant == 'sii':
        lab = 'Intensity of [SII] emission'
    if not nolog:
        lab = 'log('+lab+')'
    return lab

#================================
# get_xlabel() function
def get_xlabel(quant):
    lab = 'def'
    if quant == 'lv' or quant == 'lvc':
        lab = 'l [deg]'
    elif quant == 'sii':
        lab = 'x [pc]'
    return lab

#================================
# get_ylabel() function
def get_ylabel(quant):
    lab = 'def'
    if quant == 'lv' or quant == 'lvc':
        lab = 'v [km/s]'
    elif quant == 'sii':
        lab = 'y [pc]'
    return lab

#================================
# get_aspect() function
def get_aspect(quant):
    asp = None
    if quant =='lv' or quant == 'lvc':
        asp = 'auto'
    return asp

#================================
# get_smooth() function
def get_smooth(data,dx,smooth,**kwargs):
    # compute sigma
    sig    = smooth/dx # sig in pixels
    return gaussian_filter(data,sig) 

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

    # setup lists to store data
    grdt  = []
    yvals = [] # if quant = lv, y changes w/ time
    
    for i in range(n):
        t, x, y, data = get_data(files[i], quant=quant,**kwargs)
        # Set time dependent data 
        yvals.append(y)
        grdt.append(data)
        tarr[i] = t
           
    return tarr, x, yvals, grdt 


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
    
    # Read in system arguments
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter) 
    parser.add_argument("quant",type=str,
                        help="Plotting options:\n"
                             "  lv: HI emission f(l,v)\n")
    parser.add_argument("--anim", dest="anim",action='store_true',
                        default=False,
                        help="Switch to do animation\n")
    parser.add_argument("--iani", dest="iani",nargs=2,required=False,
                        default=[0,0],type=int,
                        help="Animate from frame iani[0] to iani[1]\n")
    parser.add_argument("--qminmax", dest="qminmax",nargs=2,required=False,
                        default=[-2,2],type=float,
                        help="Min/max value for imshow")
    parser.add_argument("--ifrm", dest="ifrm",type=int,default=0,
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
    parser.add_argument("--interp",dest="interp",type=str,default='None',
                        help="What type of interpolation for imshow?\n" 
                             "Default is 'None', cf. \n"
                             "https://matplotlib.org/devdocs/gallery/"
                             "images_contours_and_fields/interpolation_methods.html\n")
    parser.add_argument("--nolog",dest="nolog",action='store_true',
                        default=False,
                        help="Switch to take log of intensity, default is False")
    parser.add_argument("--smooth",dest="smooth",type=float,default=None,
                        help="Perform a gaussian smoothing of [SII] emission.\n"
                             "Enter smoothing scale in pc\n")
    parser.add_argument("--iline",dest="iline",type=int,default=0,
                        help="Line of interest for [SII] emission.\n"
                             "0: 6717 Angstrom \n"
                             "1: 6731 Angstrom \n"
                             "2: Line ratio \n")
    parser.add_argument("--old", dest="old",action='store_true',
                        default=False, help="Switch if plotting old sii files, containing only 6717 emission\n")

    # parsing arguments            
    args  = parser.parse_args()
    quant = args.quant
    anim  = args.anim
    iani  = args.iani
    qmin, qmax = args.qminmax[0], args.qminmax[1] 
    ifrm  = args.ifrm
    save  = args.save
    nolog = args.nolog
    interp= args.interp
    smooth= args.smooth
    iline = args.iline
    old   = args.old

    # Get otf files 
    if quant == 'lv' or quant == 'lvc':
        files = get_files(athdir,'id0','*.lv.*')
    else:
        files = get_files(athdir,'id0','*.'+quant+'.*')

    # Sort them
    files.sort()
    # Get total number of files 
    n   = len(files)

    # modify iani from default value 
    if iani[1] == 0:
        iani[1] = n-1
    # Conflicting argument checks 
    if anim and ifrm != 0:
        print("[main]: specifying a single frame while animating doesn't make sense")
        quit()
    if ifrm > (n-1):
        print("[main]: frame out of range of simulation!")
        quit()
    if smooth and quant =='lv':
        print("[main]: No need to smooth lv diagrams")
        quit()
    if smooth and quant =='lvc':
        print("[main]: No need to smooth lvc diagrams")
        quit()
    
    # Get the data
    tarr, x, y, data = init(quant,files,prec=prec,iline=iline,old=old)
    
        # Get labels 
    title = get_title(quant,nolog)
    xlab  = get_xlabel(quant)
    ylab  = get_ylabel(quant)
    asp   = get_aspect(quant)

    if anim:
        print("[main]: Animating from t = %1.1f [Myr]\n" 
              "                    to t = %1.1f [Myr]\n"
                    %( tarr[iani[0]], tarr[iani[1]] ) )  
    if interp != 'None':
        print("[main]: Interpolation is on! Using %s to"
              " interpolate\n" %( interp ) )
    # Open figure
    fig = plt.figure(figsize=(7.0,5.5),facecolor='white')
    ax1 = fig.add_subplot(111)
    
    # Handle animation
    if anim:
        # Get spacing
        dx = x[1] - x[0]
        dy = y[iani[0]][1] - y[iani[1]][0]
        #  get extent
        mnx = x[ 0] - 0.5*dx
        mxx = x[-1] + 0.5*dx
        mny = y[iani[0]][ 0] - 0.5*dy
        mxy = y[iani[1]][-1] + 0.5*dy
        
        # set mny, mxy for whole animation
        mny0 = 0.9*mny 
        mxy0 = 0.9*mxy
        
        if not nolog:
            data[iani[0]] = get_log(data[iani[0]])
        if smooth:
            data[iani[0]] = get_smooth(data[iani[0]],dx,smooth)

        im = ax1.imshow(data[iani[0]], extent=[mnx,mxx,mny,mxy],
                        origin='lower',aspect=asp,
                        vmin = qmin, vmax = qmax,
                        interpolation = interp)
        
        ax1.set_xlim(mnx,mxx)
        #if quant == 'lv':
        #    ax1.set_ylim(mny0,mxy0)
        ax1.set_xlabel(xlab)
        ax1.set_ylabel(ylab)
        ax1.set_title('t = %1.1f [Myr]' % tarr[iani[0]]) 

        div  = make_axes_locatable(ax1)
        cax  = div.append_axes('right', '5%', '5%')
        cbar = fig.colorbar(im,label=title, cax=cax)

        def anim_init():
            pass

        def animate(ifrm):
            # Clear the axes
            ax1.cla()
            # update mny, mxy
            dy  = y[ifrm][1 ] - y[ifrm][0]
            mny = y[ifrm][ 0] - 0.5*dy
            mxy = y[ifrm][-1] + 0.5*dy
            # update title
            ax1.set_title('t = %1.1f [Myr]' % tarr[ifrm])
            if ifrm != iani[0]:
                if not nolog:
                    data[ifrm] = get_log(data[ifrm])
                if smooth:
                    data[ifrm] = get_smooth(data[ifrm],dx,smooth)

            # make the plot
            im = ax1.imshow(data[ifrm], extent=[mnx,mxx,mny,mxy],
                        origin='lower',aspect=asp,
                        vmin = qmin, vmax = qmax, 
                        interpolation=interp)
            # Make sure the axes stay the same
            if quant == 'lv':
                ax1.set_ylim(mny0,mxy0)
            ax1.set_xlabel(xlab)
            ax1.set_ylabel(ylab)
            # Clear colorbar
            cax.cla()
            # Make new colorbar
            fig.colorbar(im,label=title,cax=cax)

            return #[im] 
        # do the animation
        ani = animation.FuncAnimation(fig, animate, range(iani[0],iani[1]+1), 
                                          init_func=anim_init,repeat=False)
        if save:
            print("[main]: Saving animation")
            #ani.save(quant+".gif",writer='imagemagick')
            ani.save(quant+".avi")
    else:
        # Get spacing
        dx = x[1] - x[0]
        dy = y[ifrm][1] - y[ifrm][0]
        #  get extent
        mnx = x[ 0] - 0.5*dx
        mxx = x[-1] + 0.5*dx
        mny = y[ifrm][ 0] - 0.5*dy
        mxy = y[ifrm][-1] + 0.5*dy

        if not nolog:
            data[ifrm] = get_log(data[ifrm])
        if smooth:
            data[ifrm] = get_smooth(data[ifrm],dx,smooth)

        im = ax1.imshow(data[ifrm], extent=[mnx,mxx,mny,mxy],
                        origin='lower',aspect=asp,
                        vmin = qmin,vmax = qmax,
                        interpolation = interp)
        ax1.set_xlim(mnx,mxx)
        ax1.set_ylim(mny,mxy)
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

