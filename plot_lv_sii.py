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
# scipy imports
from scipy.ndimage.filters import gaussian_filter
from scipy.interpolate import interp1d

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


    if quant == 'lv' or quant=='asym':
        t, x, y, lv, lvc = read_lv(file,prec, **kwargs)
        data = lv
    elif quant == 'lvc':
        t, x, y, lv, lvc = read_lv(file,prec, **kwargs)
        data = lvc
    elif quant == 'sii':
        t, x, y, data = read_sii(file,prec,**kwargs)
        if iline == 2:
            data = data[0]/data[1] # line ratio
            data[np.isnan(data)] = 0
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
    data[np.where(data <= 0.0)] = 1e-80
    data = np.log10(data)
    return data 

#================================
# get_title() function
def get_title(quant,nolog,**kwargs):
    lab = 'def'
    for key in kwargs:
        if key == 'iline':
            iline = kwargs[key] 
    if quant == 'lv' or quant == 'lvc':
        lab = 'Intensity of HI emission'
    elif quant == 'sii':
        if   iline == 0:
            lab = '6717 $\AA$ [SII] emission'
        elif iline == 1:
            lab = '6731 $\AA$ [SII] emission'
        elif iline == 2:
            lab = '[SII] line ratio, 6717/6731' 
        else: 
            print('[get_title]: iline not recognized \n')
            quit()
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
def init(quant, files, myfrms,  **kwargs):
    # myfrms are the frames of interest
    ctrl_files = 1
    # parse keyword args
    for key in kwargs:
        if key == 'ctrl_files':
            ctrl_files = kwargs[key]
            nctrl      = len(ctrl_files)

    # set length of simulation from myfrms
    #n    = len(files)
    n    = len(myfrms)
    # setup time array
    tarr = np.zeros(n)

    # setup lists to store data
    grdt  = []
    yvals = [] # if quant = lv, y changes w/ time
    
    i = 0
    for iff in myfrms:
        t, x, y, data = get_data(files[iff], quant=quant,**kwargs)
        # Set time dependent data 
        yvals.append(y)
        # asym measure
        if quant == 'asym':
            nrays = len(x)
            nrays2 = nrays/2 
            grdt.append( np.sum(data[:,0:nrays2]) / np.sum(data[:,nrays2:]) )
        else:
            grdt.append(data)
        tarr[i] = t
        i += 1    
    
    return tarr, x, yvals, grdt 

#==========================================================
#\func factors()
# given n, this returns the factors 
def factors(n): 
    return sorted(reduce(list.__add__,
                 ([i,n//i] for i in 
                 range(1,int(pow(n,0.5)+1)) if n%i == 0)))

#==========================================================
# \func get_levels() 
# gets the levels for the contour tracing of the cloud
# This returns levels that draw contours encompassing given
# percentages of the total intensity of emission   
def get_levels(img,pcts=np.array([0.9,0.5])):
    # img : array of cloud mass 
    # pcts: array of desired percentages (in descending order) 

    # define 1D array of values inside img 
    n    = 100
    t    = np.linspace(0, img.max(), n)  
    # Compute integral within each t 
    integral = (( img >= t[:, None, None]) * img).sum(axis=(1,2)) 
    # Interpolate on the integral 
    f        = interp1d(integral, t, fill_value='extrapolate',bounds_error=False)
    
    # As long as extrapolated values are negative, should be okay 
    levels = list(f(pcts))  

    return levels  


#==========================================================
#\func get_args()
# this function parses CMD line args
def get_args():
    # Read in system arguments
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter) 
    parser.add_argument("quant",type=str,
                        help="Plotting options:\n"
                             "  lv: HI emission f(l,v)\n"
                             " sii: [SII] emission f(x,y)\n")
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
    parser.add_argument("--ctrace",dest="ctrace",action='store_true',
                        default=False, help="Switch to overplot a contour for the cloud (l,v) emission\n")
    parser.add_argument("--vmnmx", dest="vmnmx",nargs=2,required=False,default=[-400,400],type=float,
                        help="For (l,v) diagrams, set the plotting range")
    parser.add_argument("--ipos", dest="ipos",type=int,default=0,required=False,
                        help="Observer position flag for (l,v) diagrams.\n"
                             "0: Position in disk @ Rsun\n"
                             "1: Position is from Andromeda @ angle 0  deg\n"
                             "2: Position is from Andromeda @ angle 45 deg\n"
                             "3: Position is from Andromeda @ angle 90 deg\n")
    parser.add_argument("--noplot", dest="noplot",action='store_true',
                        default=False, required=False, 
                        help="Switch to return only stitched together array\n"
                             "To be used if this file is imported from another file\n")
    return parser.parse_args() 
     
#===============================
# main()
def main(args):
    # Get input file
    athdir = os.getcwd()+'/'
    athin = [fname for fname in os.listdir(athdir) if fname.startswith("athinput")][0]
    # Read input file
    base, params = read_athinput.readath(athdir+athin) 
    # Set precision
    prec  = 32
    
    # parsing arguments            
    quant   = args.quant
    anim    = args.anim
    iani    = args.iani
    qminmax = args.qminmax 
    ifrm   = args.ifrm
    save   = args.save
    nolog  = args.nolog
    interp = args.interp
    smooth = args.smooth
    iline  = args.iline
    old    = args.old
    ctrace = args.ctrace
    vmnmx  = args.vmnmx
    ipos   = args.ipos
    noplot = args.noplot
    
    # Get panel flag 
    pflag = True if np.size(ifrm) > 1 else False 
    # Get qminmax flag
    qflag = True if np.size(qminmax) > 1 else False 
    if np.size(ifrm) == 1: ifrm = ifrm[0] 

    # Get lv files 
    if quant == 'lv' or quant == 'lvc':
        files = get_files(athdir,'id0','*.lv.*')
        if qflag:
            qmin = qminmax[0]
            qmax = qminmax[1]
        else:
            qmin = -5
            qmax =  5

    elif quant == 'sii':
        files = get_files(athdir,'id0','*.'+quant+'.*')
        # Change default range for sii files
        if qflag:
            qmin = qminmax[0]
            qmax = qminmax[1]
        else:
            qmin = -15
            qmax = -5
    elif quant == 'asym':
        files = get_files(athdir,'id0','*.lv.*')
        ifrm
        
    else: 
        print('[main]: quant not understood, aborting...\n')
        quit()

    # Sort them
    files.sort()
    # Get total number of files 
    n   = len(files)

    # modify iani from default value 
    if iani[1] == 0:
        iani[1] = n
    
    # Determine myframes
    if anim:
        myfrms = range(iani[0],iani[1])
    elif pflag: 
        myfrms = ifrm
    elif quant == 'asym':
        myfrms = range(0,n) 
    else:
        myfrms = [ifrm] 


    # Conflicting argument checks 
    #if anim and myfrms != 0:
    #    print("[main]: specifying a single frame while animating doesn't make sense")
    #    quit()
    if myfrms[-1] > (n-1):
        print("[main]: frame out of range of simulation!")
        quit()
    if smooth and quant =='lv':
        print("[main]: No need to smooth lv diagrams")
        quit()
    if smooth and quant =='lvc':
        print("[main]: No need to smooth lvc diagrams")
        quit()
    if not noplot and quant == 'asym':
        print("[main]: Asymmetry measure cannot be plotted from this file!")
        quit() 
    
    # Get the data
    tarr, x, y, data = init(quant,files,myfrms,prec=prec,iline=iline,old=old,ipos=ipos)
    if ctrace: 
        tarr, x, y, cloud = init('lvc',files,myfrms,prec=prec,iline=iline,old=old,ipos=ipos)
        # Set global values for drawing contours  
        alpha   = 1.0
        colors=['#40E0D0','#00C957']#,'#FF3030']
        pcts = np.array([0.9,0.5]) # percentages for contours 
        if pflag:
            lw = 1.
        else: 
            lw = 2.

    if noplot:
        if quant == 'asym':
            return tarr, data # (l,v) asymmetry measure (sum of intensity divided on both sides)
        elif ctrace:
            return tarr, x, y, data, cloud
        else: 
            return tarr, x, y, data 
        quit() 
        
        # Get labels 
    title = get_title(quant,nolog,iline=iline)
    xlab  = get_xlabel(quant)
    ylab  = get_ylabel(quant)
    asp   = get_aspect(quant)

    if anim:
        print("[main]: Animating from t = %1.1f [Myr]\n" 
              "                    to t = %1.1f [Myr]\n"
                    %( tarr[0], tarr[-1] ) )  
    if interp != 'None':
        print("[main]: Interpolation is on! Using %s to"
              " interpolate\n" %( interp ) )



    if not pflag:
        # Open figure
        fig = plt.figure(figsize=(7.0,5.5),facecolor='white')
        ax1 = fig.add_subplot(111, axisbg='black')
    
    # Handle animation
    if anim:
        # Get spacing
        dx = x[1] - x[0]
        dy = y[0][1] - y[0][0]
        #  get extent
        mnx = x[ 0] - 0.5*dx
        mxx = x[-1] + 0.5*dx
        mny = y[0][ 0] - 0.5*dy
        mxy = y[0][-1] + 0.5*dy
        
        # set mny, mxy for whole animation
        if quant == 'lv' or quant == 'lvc':
            mny0 = vmnmx[0] 
            mxy0 = vmnmx[1]
            # Also set the background color
            #ax1.set_facecolor('black') 
      
        if smooth:
            data[0] = get_smooth(data[0],dx,smooth)
       
        if not nolog:
            data[0]  = get_log(data[0])
            if ctrace:
            # Normalize cloud
                cloud[0] /= np.sum(cloud[0]) 
        
        im = ax1.imshow(data[0], extent=[mnx,mxx,mny,mxy],
                        origin='lower',aspect=asp,
                        vmin = qmin, vmax = qmax,
                        interpolation = interp)

        if ctrace:
            m = np.amax(cloud[0])
            if m != 0:
                ax1.contour(cloud[0],levels=get_levels(cloud[0],pcts),
                             extent=[mnx,mxx,mny,mxy],
                             colors=colors,origin='lower',alpha=alpha,linewidths=lw)
    
        ax1.set_xlim(mnx,mxx)
        if quant == 'lv':
            ax1.set_ylim(mny0,mxy0)
        ax1.set_xlabel(xlab)
        ax1.set_ylabel(ylab)
        ax1.set_title('t = %1.1f [Myr]' % tarr[0]) 

        div  = make_axes_locatable(ax1)
        cax  = div.append_axes('right', '5%', '5%')
        cbar = fig.colorbar(im,label=title, cax=cax)

        def anim_init():
            pass

        def animate(ifrm):
            # Clear the axes
            ax1.cla()
            # update title
            ax1.set_title('t = %1.1f [Myr]' % tarr[ifrm])
            
            # Update extent 
            dy  = y[ifrm][1 ] - y[ifrm][0]
            mny = y[ifrm][ 0] - 0.5*dy
            mxy = y[ifrm][-1] + 0.5*dy

            if ifrm != iani[0]:
                if smooth:
                    data[ifrm] = get_smooth(data[ifrm],dx,smooth)
                if not nolog:
                    data[ifrm]  = get_log(data[ifrm])
                    if ctrace:
                        # Normalize cloud
                        cloud[ifrm] /= np.sum(cloud[ifrm]) 
            # make the plot
            im = ax1.imshow(data[ifrm], extent=[mnx,mxx,mny,mxy],
                        origin='lower',aspect=asp,
                        vmin = qmin, vmax = qmax, 
                        interpolation=interp)

            if ctrace:
                m = np.amax(cloud[ifrm])
                if m != 0:
                    ax1.contour(cloud[ifrm],levels=get_levels(cloud[ifrm],pcts),
                                 extent=[mnx,mxx,mny,mxy],
                                 colors=colors,origin='lower',alpha=alpha,linewidths=lw)
            # Make sure the axes stay the same
            if quant == 'lv':
                ax1.set_ylim(mny0,mxy0)
            ax1.set_xlabel(xlab)
            ax1.set_ylabel(ylab)
            # Clear colorbar
            cax.cla()
            # Make new colorbar
            fig.colorbar(im,label=title,cax=cax)
            return  

        # do the animation
        ani = animation.FuncAnimation(fig, animate, range(len(myfrms)), 
                                          init_func=anim_init,repeat=False)
        if save:
            print("[main]: Saving animation")
            ani.save(quant+".gif",writer='imagemagick')
            #ani.save(quant+".avi")

    # Handle making panel plots 
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
        
        # define ratio
        ratio = float(nxp)/float(nyp) 
        fsize = (ratio*7., 7.)
        # Create figure object
        fig, axes = plt.subplots(nyp,nxp, sharex='col', sharey='row',facecolor='white',
                                 figsize=(ratio*7.,7.)) 
        fig.subplots_adjust(hspace=0.1, wspace=0.1)
        # define the colorbar
        fig.subplots_adjust(top=0.8) 
        cax = fig.add_axes([0.16,0.85,0.7,0.02])

        # Set spacing stuff
        dx  = x[1] - x[0]
        mnx = x[0 ] - 0.5*dx
        mxx = x[-1] + 0.5*dx
        if quant == 'lv' or quant == 'lvc':
            mny0 = vmnmx[0]
            mxy0 = vmnmx[1]


        # Now plot
        for (ax, iff) in zip(axes.flat, range(len(myfrms))):
            ax.set_axis_bgcolor('black')
            # Get spacing
            dy = y[iff][1] - y[0][0]
            #  get extent
            mny = y[iff][ 0] - 0.5*dy
            mxy = y[iff][-1] + 0.5*dy
            
            # Take log
            if not nolog:
                data[iff]  = get_log(data[iff]) 
                if ctrace:
                    m = np.amax(cloud[iff])
                    if m != 0:
                        # Normalize cloud
                        cloud[iff] /= np.sum(cloud[iff]) 


            im = ax.imshow(data[iff], extent=[mnx,mxx,mny,mxy],
                            origin='lower',aspect=asp,
                            vmin = qmin, vmax = qmax, 
                            interpolation=interp)
            if ctrace: 
                m = np.amax(cloud[iff])
                if m != 0:
                    ax.contour(cloud[iff],levels=get_levels(cloud[iff],pcts),
                                 extent=[mnx,mxx,mny,mxy],
                                 colors=colors,origin='lower',alpha=alpha,linewidths=lw)
            ax.set_xlim(mnx ,mxx )
            if quant == 'lv':
                ax.set_ylim(mny0,mxy0)
                mxy = mxy0 
            ax.text(0.9*mnx, 0.8*mxy, 't = %1.1f [Myr]' % (tarr[iff]),
                        bbox={'facecolor':'white', 'alpha':0.9, 'pad':5})
            if quant == 'lv' and ipos == 0:
                ax.set_xticks([-100,0,100])
            
        # define global labels
        fig.text(0.5, 0.04, xlab, ha='center')
        fig.text(0.03, 0.5, ylab, va='center', rotation='vertical') 

        # Set colorbar
        cb = fig.colorbar(im,cax=cax, orientation='horizontal') 
        cax.xaxis.set_ticks_position('bottom') 
        cb.ax.set_title(title) 

        
    # Handle plotting single frame 
    else:
        # Get spacing
        dx = x[1] - x[0]
        dy = y[0][1] - y[0][0]
        #  get extent
        mnx = x[ 0] - 0.5*dx
        mxx = x[-1] + 0.5*dx
        mny = y[0][ 0] - 0.5*dy
        mxy = y[0][-1] + 0.5*dy

        if quant == 'lv' or quant == 'lvc':
            mny0 = vmnmx[0]
            mxy0 = vmnmx[1]
        
        if smooth:
            data[0] = get_smooth(data[0],dx,smooth)
        if not nolog:
            data[0]  = get_log(data[0])
            if ctrace: 
                # Normalize cloud
                cloud[0] /= np.sum(cloud[0])

        im = ax1.imshow(data[0], extent=[mnx,mxx,mny,mxy],
                        origin='lower',aspect=asp,
                        vmin = qmin,vmax = qmax,
                        interpolation = interp)
        if ctrace:
            m = np.amax(cloud[0])
            if m != 0:
                ax1.contour(cloud[0],levels=get_levels(cloud[0],pcts),
                             extent=[mnx,mxx,mny,mxy],
                             colors=colors,origin='lower',alpha=alpha,linewidths=lw)

        ax1.set_xlim(mnx,mxx)
        ax1.set_ylim(mny0,mxy0)
        ax1.set_xlabel(xlab)
        ax1.set_ylabel(ylab)
        ax1.set_title('t = %1.1f [Myr]' % tarr[0]) 
        plt.colorbar(im,label=title)

        if save:
            print("[main]: Saving figure")
            plt.savefig(quant+str(ifrm)+".eps")
    
    if save:
        print("[main]: Program complete")
    else:
        plt.show()
    #================================================================================#
            
if __name__ == '__main__':
    # plot_lv_sii.py called from cmd line
    args = get_args()
    main(args) 

