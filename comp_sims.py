#!/usr/bin/env python
import numpy as np
import sys
import os
import subprocess as sbp
import argparse
from argparse import RawTextHelpFormatter
import matplotlib.pyplot as plt
import matplotlib as mpl
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

import plot_sims as ps
import plot_lv_sii as plv
import plot_otf as potf

#=====================================================
#
#  Code: comp_sims.py
#
#  Purpose:  Compares and plots simulations  
#
#  Keywords: python comp_sims.py -h   
#
#  Usage: python comp_sims.py   
#          
#  Author: John Dupuy 
#          UNC Chapel Hill
#  Date:    04/13/18 
#  Updated: 04/13/18 
#=====================================================

#=================================
# get_dirs() function
def get_dirs():
    # Purpose: get the directories of each simulation
    #          and store them in a list
    # input: 
    #         cwd: current working directory
    # output: 
    #         a list of strings that contain
    #         all the directories of simulations
    
    # Use the find function to get directory structure
    process = sbp.Popen(["find",os.getcwd(),"-type","d"], stdout=sbp.PIPE)
    stdout, stderr = process.communicate()

    dirlist = []
    
    print('[get_dirs]: Finding simulation directories...')
    # Loop through stdout and get sims
    for d in stdout.splitlines():
        # The slice is to avoid appending processor directories 
        if ('fac' in d) and ('id' not in d[-10:]):
            dirlist.append(d+'/')

    return dirlist 

def get_args():
    # Read in system arguments
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter) 
    # Arguments for (l,v) and (sim) plotting
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
    parser.add_argument("--comp", dest="comp",action='store_true',
                        default=False, required=False,
                        help="Comparison plot, currently need to specify whole path"
                             "for the comparison of interest.\n"
                             " ONLY FOR 1D PLOTS!\n")
    parser.add_argument("--vvec",dest="vvec",action='store_true',
                        default=False, required=False,
                        help="Overplot velocity vectors\n")
    parser.add_argument("--ms",dest="ms",action='store_true',
                        default=False, required=False,
                        help="Marker for the sun's position\n")
    parser.add_argument("--rtrace",dest="rtrace",action='store_true',
                        default=False, required=False,
                        help="Switch to trace out rays for (l,v) diagrams") 
    parser.add_argument("--noplot",dest="noplot",action='store_true',
                        default=True, required=False,
                        help="Switch to return only stitched together array\n"
                             "To be used if this file is imported from another file\n")
    parser.add_argument("--com", dest="com", action='store_true',
                    default=False, required=False,
                    help="Switch to compute cloud center of mass as a function of time, and\n"
                        "plot it.")
    parser.add_argument("--stat",dest="stat",
                        default='None', type=str,
                        help="Type of statistic to compute error bars in CoM measurement") 
    # Arguments for (l,v) and sii plotting
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
    parser.add_argument("--vmnmx", dest="vmnmx",nargs=2,required=False,default=[-400,400],type=float,
                        help="For (l,v) diagrams, set the plotting range")
    parser.add_argument("--ipos", dest="ipos",type=int,default=0,required=False,
                        help="Observer position flag for (l,v) diagrams.\n"
                             "0: Position in disk @ Rsun\n"
                             "1: Position is from Andromeda @ angle 0  deg\n"
                             "2: Position is from Andromeda @ angle 45 deg\n"
                             "3: Position is from Andromeda @ angle 90 deg\n")
    parser.add_argument("--sumv", dest="sumv",action='store_true',
                    default=False, required=False,
                    help="Sum over velocity bins to get I(l)")                 
    # Arguments for otf plotting 
    parser.add_argument("--rmnmx", dest="rmnmx",nargs=2,required=False,
                    default=[5000., 13000.],type=float,
                    help="Average A1 from rmn to rmx") 
    parser.add_argument("--iang", dest="iang",type=int, default=0,
                    required=False,
                    help="Angle to plot if using vrot:\n"
                         "  [0-(Nang-1)]: index that corresponds to [0-305 deg]\n"
                         " -1: All angles are plotted\n")
    parser.add_argument("--cmap", dest="cmap",action='store_true',
                        default=False,
                        help="Switch to create a colormap for A1 & A2")
    parser.add_argument("--wire", dest='wire',action='store_true',
                        default=False,
                        help="Switch to make a 3d wireframe plot")
    parser.add_argument("--fit", dest="fit",type=str,
                        default='None',
                        help="Switch to do fit for <A1> perturbations")
    parser.add_argument("--ctrl", dest="ctrl",type=str,
                        default='None',
                        help="String to compare measurement to control files:\n"
                             "  sub: Subtract control from sim and take abs value\n"
                             "  div: Divide sim by control")
    # Arguments for selecting simulations to plot
    parser.add_argument("--rpos",dest="rpos",type=str, default=['5500'],
                        required=False,nargs='+',
                        help="Position of cloud, can be any number of entries so\n"
                             "as to compare different position.\n"
                             " MUST correspond to directory however")
    parser.add_argument("--rc",dest="rc",type=str, default=['500'],
                        required=False,nargs='+',
                        help="Radius of cloud, can be any number of entries so\n"
                             "as to compare different position.\n"
                             " MUST correspond to directory however")
    parser.add_argument("--fac",dest="fac",type=str, default=['-1.0'],
                        required=False,nargs='+',
                        help="Factor for cloud, can be any number of entries so\n"
                             "as to compare different position.\n"
                             " MUST correspond to directory however")
    parser.add_argument("--ang",dest="ang",type=str, default=['0.0'],
                        required=False,nargs='+',
                        help="Angle for cloud, can be any number of entries so\n"
                             "as to compare different position.\n"
                             " MUST correspond to directory however")
    parser.add_argument("--mc",dest="mc",type=str, default=['1e7'],
                        required=False,nargs='+',
                        help="Mass for cloud, can be any number of entries so\n"
                             "as to compare different position.\n"
                             " MUST correspond to directory however")
    parser.add_argument("--te",dest="te",type=str,default=['325'],
                        required=False,nargs='+',
                        help="Ending time for cloud collision")

    return parser.parse_args() 

#\func factors()
# given n, this returns the factors 
def factors(n): 
    return sorted(reduce(list.__add__,
                 ([i,n//i] for i in 
                 range(1,int(pow(n,0.5)+1)) if n%i == 0)))

def get_sims(args, cwd = os.getcwd()):
    # Parse relevant args for simulation selection
    rpos     = args.rpos
    rhvc     = args.rc
    facvhvc  = args.fac
    ahvc     = args.ang
    mhvc     = args.mc
    thvce    = args.te

    # Construct simulation names 
    cwd += '/' 

    # Define a list to store directory names
    mysims = []

    # Now loop over parameters and construct directory names
    print("[get_sims]: Constructing relevant simulation directories...") 
    for m in mhvc:
        for te in thvce:
            for rc in rhvc:
                for r in rpos:
                    for a in ahvc:
                        for fac in facvhvc:
                            simname = cwd+'m'+m+'/te'+te+'/rc'+rc+'/r'+r+'/a'+a+'/fac'+fac+'/'
                            mysims.append(simname) 
    return mysims 

def get_data(args, sims):
    # Parse relevant arguments for the data
    quant = args.quant
    # Select proper plotting file based on quant

    # Using lv
    if (quant=='lv' or   quant=='lvc' or 
        quant=='sii' or quant=='asym'):

        if quant=='asym':
            lvflag = False
        elif args.sumv:
            lvflag = False 
        else:
            lvflag  = True
        otfflag = False
        reader  = plv
    # Using otf
    elif (quant=='<A1>' or quant=='<A2>' or
          quant=='LoR'  or quant=='RoL'  or
          quant=='mcR'  or quant=='mcL'):

        lvflag  = False
        otfflag = True
        reader  = potf
    # Using sims
    else: 
        lvflag  = False 
        otfflag = False
        reader  = ps

    # Create data list 
    data_list = []
    print('[get_data]: Reading simulations files...') 
    # Now loop over simulations and read data
    for s in sims:
        os.chdir(s)
        try:
            data_list.append(reader.main(args)) 
        except OSError:
            print('[get_data]: Sim files not present for %s!' %(s))

    return lvflag, otfflag, data_list  

def make_plots(args,lvflag,otfflag,data,sims):
    # Specify the control simulation directory
    ctrl_path = '/srv/scratch/jdupuy26/jid/jid_adiabatic'\
                '/longevity_jid/m0.0/te325/rc500/r10000/a0.0'\
                '/fac1.0'

    #-----------------------------------------#
    # First parse relevant plotting arguments
    vmnmx    = args.vmnmx
    interp   = args.interp
    qminmax  = args.qminmax
    mnx, mxx = args.mnmx
    ctrace   = args.ctrace
    qflag    = True if np.size(qminmax) > 1 else False 
    quant    = args.quant
    nolog    = args.nolog
    iline    = args.iline
    ipos     = args.ipos
    ifrm     = args.ifrm
    cart     = args.cart
    log      = args.log
    iunit    = args.units  
    com      = args.com
    sumv     = args.sumv
    mnx, mxx = args.mnmx
    mny, mxy = mnx, mxx

    # Parse relevant args for simulation selection
    rpos     = args.rpos
    rhvc     = args.rc
    facvhvc  = args.fac
    ahvc     = args.ang
    mhvc     = args.mc
    thvce    = args.te

    # Determine proper labels 
    if len(rpos) > 1:
        var  = rpos
        pstr  = '$r_{\\rm pos,0}$ = '
        punit = ' [pc]'  
    elif len(facvhvc) > 1:
        var  = facvhvc
        pstr  = '$f_c$ = '
        punit = ''
    elif len(rhvc) > 1:
        var = rhvc
        pstr  = '$r_{c}$ = '
        punit = ' [pc]'
    elif len(ahvc) > 1:
        var = ahvc
        pstr  = '$\\phi_{pos,0}$ = '
        punit = ' [rad]'
    elif len(mhvc) > 1:
        var = mhvc
        pstr  = '$M_{c}$ = '
        punit = ' [M$_{\\odot}$]'
    else:
        var = thvce
        pstr  = '$\\Delta t$ = '
        punit = ' [Myr]' 

    # ctrace stuff
    if ctrace:
        colors = ['#40E0D0','#00C957']  
        lw     = 1.
        alpha  = 1.
        
    #-----------------------------------------#

    #-----------------------------------------#
    # Create the proper figure
    if otfflag or quant == 'asym' or com:
        fig = plt.figure(figsize=(10.0,7.0),facecolor='white')
    else:   
        # Determine the size of the panel
        fact = factors(len(data))
        if len(fact) == 2:
            nxp = np.max(fact)
            nyp = np.min(fact) 
        else:
            nxp = np.max(fact[1:-1])
            nyp = np.min(fact[1:-1])
        ratio = float(nxp)/float(nyp)
        fsize = (ratio*7.,7.)
        # Create figure object
        fig, axes = plt.subplots(nyp,nxp,sharex='col', sharey='row',facecolor='white',
                                 figsize=fsize) 

        fig.subplots_adjust(hspace=0.1, wspace=0.1)
        # define the colorbar
        fig.subplots_adjust(top=0.8) 
        cax = fig.add_axes([0.16,0.85,0.7,0.02])
    #-----------------------------------------#

    
    #-----------------------------------------#
    # Handle making (l,v) plots first 
    if lvflag:
        if qflag:
            qmin = qminmax[0]
            qmax = qminmax[1]
        else:
            qmin = -5
            qmax =  5

        # Get labels 
        title = plv.get_title(quant,nolog,iline=iline)
        xlab  = plv.get_xlabel(quant)
        ylab  = plv.get_ylabel(quant)
        asp   = plv.get_aspect(quant) 

        # Animation functionality? How to implement, do this next week 
        for iff in range(len(ifrm)):
            # Now plot
            for (ax, dat, v) in zip(axes.flat, data, var):
                # Set background color for image
                ax.set_axis_bgcolor('black')
                # Parse data
                if ctrace:
                    tarr, x, y, imgs, imgc = dat  
                else: 
                    tarr, x, y, imgs       = dat


                # Set spacing stuff
                dx  = x[1 ] - x[0]
                mnx = x[0 ] - 0.5*dx
                mxx = x[-1] + 0.5*dx
                if quant == 'lv' or quant == 'lvc':
                    mny0 = vmnmx[0]
                    mxy0 = vmnmx[1]
                # Get spacing
                dy = y[iff][1] - y[iff][0]
                #  get extent
                mny = y[iff][ 0] - 0.5*dy
                mxy = y[iff][-1] + 0.5*dy

                # Take log
                if not nolog:
                    imgs[iff]  = plv.get_log(imgs[iff]) 
                    if ctrace:
                        m = np.amax(imgc[iff])
                        if m != 0:
                            # Normalize cloud
                            imgc[iff] /= np.sum(imgc[iff]) 
                
                im = ax.imshow(imgs[iff], extent=[mnx,mxx,mny,mxy],
                                origin='lower',aspect=asp,
                                vmin = qmin, vmax = qmax, 
                                interpolation=interp)
                if ctrace: 
                    if tarr[iff] > params[0].thvcs:
                        ax.contour(imgc[iff],levels=plv.get_levels(imgc[iff]),
                                     extent=[mnx,mxx,mny,mxy],
                                     colors=colors,origin='lower',alpha=alpha,linewidths=lw)
                ax.set_xlim(mnx ,mxx )
                if quant == 'lv':
                    ax.set_ylim(mny0,mxy0)
                    mxy = mxy0 
                ax.text(0.9*mnx, 0.8*mxy, '%s%s%s' % (pstr,v,punit),
                            bbox={'facecolor':'white', 'alpha':0.9, 'pad':5})
                if quant == 'lv' and ipos == 0:
                    ax.set_xticks([-100,0,100])
                
                # define global labels
                fig.text(0.5, 0.04, xlab, ha='center')
                fig.text(0.03, 0.5, ylab, va='center', rotation='vertical') 

            # Set colorbar
            cb = fig.colorbar(im,cax=cax, orientation='horizontal') 
            cax.xaxis.set_ticks_position('bottom') 
            cb.ax.set_title(title + ' @ t = %1.1f [Myr]' % (tarr[iff]) ) 

    # Now handle otf figures
    elif otfflag or quant=='asym' or com: 
        # Define labels 
        if quant == 'asym':
            ylab = '(l,v) asymmetry measure'
            step = 1
        elif com:
            ylab = '$R_{\\rm com}$ [kpc]'
            step = 1

        else:
            ylab = potf.get_ylabel(quant)
            step = 5
        xlab = 't [Myr]' 
        # Define cutoff (if plotting '<A1>') 
        cut  = 0.05
        # Define plotting markers
        #colors = ['b-','m-','k-','g-']
        mycmap = mpl.cm.get_cmap('viridis')
        colors = [ mycmap(x) for x in np.linspace(0.0, 0.8, 4)]
        lines  = ['-','-','-','-']
        markers= ['s','p','^','o']

        # Put labels on plot
        plt.xlabel(xlab,fontsize=16)
        plt.ylabel(ylab,fontsize=16)

        for (dat, v, c, l, m) in zip(data, var, colors, lines, markers): 
            if com:
                tarr, rlo, rhi, vals = dat
            else:
                tarr, vals = dat 
            plt.plot(tarr[::step], vals[::step], linestyle=l, color = c, label= '%s%s%s' % (pstr,v, punit),
                     linewidth=2., marker=m, markersize=7.)  

            # Error on CoM
            if com: 
                plt.fill_between(tarr[::step], rlo[::step], rhi[::step], alpha=0.2,antialiased=True, color=c)
        
        if quant == '<A1>':
            # Plot cutoff value
            plt.plot(tarr, cut*np.ones(len(tarr)), 'r--', lw=3.)
        #plt.ylim(0, 0.3)
        plt.legend(loc=2) 
    
    # Now handle plotting simulations 
    else: 
        # Get labels
        xlab, ylab, clab = ps.get_labels(quant,iunit,log) 
        # Base parameters off of control simulation 
        base, params = ps.get_athinput(ctrl_path)
        # get face centered meshgrid
        x1f, x2f           = ps.get_fc(cart,params) 
        # get cartesian version of face centered meshgrid
        if cart:
            x1cf, x2cf = x1f, x2f
        else:
            x1cf, x2cf = x1f*np.cos(x2f), x1f*np.sin(x2f) 
        # get cartesian version of 'cell centered' meshgrid
        x1c, x2c           = ps.get_fc2(cart,params)
        if cart:
            x1cc, x2cc = x1c, x2c
        else:
            x1cc, x2cc = x1c*np.cos(x2c),  x1c*np.sin(x2c)

        # Animation functionality? Need to implement later 
        for iff in range(len(ifrm)):
            for (ax, dat, v) in zip(axes.flat, data, var):
                # Parse data
                if ctrace:
                    tarr, x, y, imgs, imgc = dat  
                else: 
                    tarr, x, y, imgs       = dat

                # Select qmin/qmax
                if qflag:
                    qmin, qmax = qminmax[0], qminmax[1]
                else:
                    qmin, qmax = np.min(imgs), np.max(imgs)
                
                # Now plot
                im = ax.pcolorfast(x1cf,x2cf,imgs[iff],vmin=qmin,vmax=qmax) 
                if ctrace:
                    if tarr[iff] > params[0].thvcs:  
                        imc  = ax.contour(x1cc,x2cc,imgc[iff],levels=ps.get_levels(imgc[iff])
                                                             ,colors=colors,alpha=alpha,linewidths=lw)
                ax.set_xlim(mnx,mxx)
                ax.set_ylim(mnx,mxx)
                #ax.set_aspect('equal')  
                ax.text(0.9*mnx, 0.8*mxy, '%s%s%s' % (pstr,v,punit),
                            bbox={'facecolor':'white', 'alpha':0.9, 'pad':5})
                # This makes .eps files manageable 
                im.set_rasterized(True) 

            # define global labels
            fig.text(0.5, 0.04, xlab, ha='center')
            fig.text(0.03, 0.5, ylab, va='center', rotation='vertical') 

            # Set colorbar
            cb = fig.colorbar(im,cax=cax, orientation='horizontal') 
            cax.xaxis.set_ticks_position('bottom') 
            cb.ax.set_title(clab + ' @ t = %1.1f [Myr]' %(tarr[iff]) ) 

    return 

def main():
    cwd  = os.getcwd() 
    args = get_args() 
    
    # Now based off of simulation selection, get relevant directories 
    sims = get_sims(args, cwd)  
    # Now get the data for each simulation 
    lvflag, otfflag, data = get_data(args,sims)
    # Now construct the plots 
    make_plots(args,lvflag,otfflag,data,sims) 
    plt.show() 

    return 

main() 
