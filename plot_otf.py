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
# Import from correct directory
sys.path.insert(0,'/afs/cas.unc.edu/users/j/d/'
                  'jdupuy26/Johns_work/'
                  'misc_scripts/plotting_tools')
import read_otf as rotf
import read_athinput
from units_class import * 


#=====================================================
#
#  Code: plot_otf.py
#
#  Purpose: Plot files read in from otf dumps 
#
#  Keywords: python plot_otf.py -h  
#
#  Usage: python plot_otf.py quant (--ifrm or --anim) --iang  
#          
#  Author: John Dupuy 
#          UNC Chapel Hill
#  Date:   11/7/17
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
    
    t, mhvc, rhvc, rpos,\
    acc_rate, facvhvc, ahvc,\
    mcR, mcL,\
    r, ang, vrot = rotf.readbin(file,prec, **kwargs)

    if iunit == 1: 
        u = units_CGS()
    elif iunit == 2:
        u = units_SI()
    else: 
        u = units_COMP()

    # Do unit conversion
    t        *= u.myr
    mhvc     *= u.msun
    rhvc     *= u.pc
    rpos     *= u.pc
    acc_rate *= u.msun/u.myr
    mcR      *= u.msun
    mcL      *= u.msun
    r        *= u.pc
    ang      *= 1.     # [in rad]
    vrot     *= u.v

    return t, mhvc, rhvc, rpos,\
           acc_rate, facvhvc, ahvc,\
           mcR, mcL,\
           r, ang, vrot 

#================================
# get_ylabel() function
def get_ylabel(quant):
    lab = 'def'
    if quant == 'mcR':
        lab = 'mass in right circle [M$_{\\odot}$]'
    elif quant == 'mcL':
        lab = 'mass in left circle [M$_{\\odot}$]'
    elif quant == 'RoL':
        lab = 'mcR/mcL [unitless]'
    elif quant == 'LoR':
        lab = 'mcL/mcR [unitless]'
    else: 
        lab = 'v [pc/Myr]'
    return lab
        

#================================
# initialize() function
def init(quant, files, **kwargs):
    Nang    = 8
    nx1_dom = 256
    # parse keyword args
    for key in kwargs:
        if key == 'Nang':
            Nang    = kwargs[key]
        if key == 'nx1_dom':
            nx1_dom = kwargs[key]

    # set length of simulation from file list
    n    = len(files)
    # setup time array
    tarr = np.zeros(n)

    # choose quant
    if quant != "vrot":
        mcR_arr = np.zeros(n)
        mcL_arr = np.zeros(n)
        
        # Read in data
        for i in range(n):
            t, mhvc, rhvc, rpos,\
            acc_rate, facvhvc, ahvc,\
            mcR, mcL,\
            r, ang, vrot        = get_data(files[i], **kwargs)
            
            tarr[i]     = t
            mcR_arr[i]  = mcR
            mcL_arr[i]  = mcL

        # print HVC parameters 
        print("\n[init]: mhvc = %1.3e [M_sun]\n"
              "        rhvc = %1.3e [pc]\n"
              "        rpos = %1.3e [pc]\n"
              "    acc_rate = %1.3e [M_sun/Myr]\n"
              "     facvhvc = %1.3e [unitless] \n"
              "        ahvc = %1.3e [rad]\n" % (mhvc, rhvc, rpos, acc_rate,facvhvc,ahvc))
        # Return quantities
        if quant == "mcR":
            return tarr, mcR_arr
        elif quant == "mcL":
            return tarr, mcL_arr
        elif quant == "LoR":
            return tarr, mcL_arr/mcR_arr
        else:
            return tarr, mcR_arr/mcL_arr
    else: 
        vrot_arr = np.zeros((n,Nang,nx1_dom))
        for i in range(n):
            t, mhvc, rhvc, rpos,\
            acc_rate, facvhvc, ahvc,\
            mcR, mcL,\
            r, ang, vrot        = get_data(files[i], **kwargs)
            
            tarr[i]     = t
            vrot_arr[i] = vrot
        # print HVC parameters 
        print("\n[init]: mhvc = %1.3e [M_sun]\n"
              "        rhvc = %1.3e [pc]\n"
              "        rpos = %1.3e [pc]\n"
              "    acc_rate = %1.3e [M_sun/Myr]\n"
              "     facvhvc = %1.3e [unitless] \n"
              "        ahvc = %1.3e [rad]\n" % (mhvc, rhvc, rpos, acc_rate,facvhvc,ahvc))
        data = (r, ang, vrot_arr)
        return tarr, data

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

    # Set Nang and nx1_dom
    Nang, nx1_dom = int(params[0].Nang), int(params[0].nx1)

    # Get otf files 
    otf_files = get_files(athdir,'id0','*.otf.*')
    n   = len(otf_files)

    # Read in system arguments
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter) 
    parser.add_argument("quant",type=str,
                        help="Plotting options:\n"
                             "  mcR: mass w/in circleR\n"
                             "  mcL: mass w/in circleL\n"
                             "  RoL: mcR/mcL\n"
                             "  LoR: mcL/mcR\n"
                             "  vrot: rotation curve\n")
    parser.add_argument("--iang", dest="iang",type=int, default=None,
                        required=False,
                        help="Angle to plot if using vrot:\n"
                             "  [0-(Nang-1)]: index that corresponds to [0-305 deg]\n"
                             " -1: All angles are plotted\n")
    parser.add_argument("--anim", dest="anim",action='store_true',
                        default=False,
                        help="Switch to do animation\n")
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
    iang  = args.iang
    anim  = args.anim
    ifrm  = args.ifrm
    save  = args.save
    
    # Conflicting argument checks 
    if iang > Nang:
        print("[main]: iang must be less than Nang!")
        quit()
    if anim and args.quant != 'vrot':
        print("[main]: only able to animate vrot!")
        quit()
    if anim and ifrm != None:
        print("[main]: specifying a single frame while animating"
              "        doesn't make sense")
        quit()
    if ifrm > (n-1):
        print("[main]: frame out of range of simulation!")
        quit()
    
    # Get the data
    tarr, data = init(quant, otf_files, prec=32, 
                      nx1_dom=nx1_dom,Nang=Nang) 

    #===== PLOTTING =============================================
    # set colors for lines
    colormap = cm.gist_ncar
    if Nang != 4:
        c    = [colormap(i) for i in np.linspace(0,0.4,Nang)]
    else: c = ['k','b','g','r'] 
    
    # set ylabel
    ystr = get_ylabel(quant)
    # for the mass stuff
    if quant != 'vrot':
        # Do plotting
        plt.figure(figsize=(10,8))
        plt.plot(tarr, data)
        plt.xlabel('t [Myr]')
        plt.ylabel(ystr) 
        if save:
            print("[main]: Saving figure")
            plt.savefig(quant+".eps")

    # for rotation curves
    else:
        xstr = 'r [pc]'
        # unpack data
        r, ang, vrot = data[0], data[1], data[2] 
        # convert ang to degrees
        ang *= 180./np.pi

        # Handle plotting all angles
        if iang == -1:
            fig, axs = plt.subplots(nrows=Nang/2,ncols=2, figsize=(10,8),
                                    sharex=True)
            # make into 1D array
            axs = np.ravel(axs)
            # a single frame
            if ifrm is not None:
                t = tarr[ifrm]
                # loop over axes and angles
                for i in range(Nang):
                    axs[i].plot(r, vrot[ifrm,i],color=c[i],
                                  label='$\\theta$ = %1.1f$^\\circ$' % ang[i])
                    axs[i].legend(loc=4)
                # make labels 
                # xlabels
                axs[-2].set_xlabel(xstr)
                axs[-1].set_xlabel(xstr)
                # common ylabel
                fig.text(0.04, 0.5, ystr, va='center',rotation='vertical')
                # title
                plt.suptitle('t = %1.1f [Myr]' % tarr[ifrm], size='large')
                if save:
                    print("[main]: Saving figure")
                    plt.savefig(quant+str(ifrm)+".eps")
            
            # handle animation   
            elif anim:
                # create list of lines 
                lines = [axs[iang].plot(r, vrot[0,iang],color=c[iang],
                            label='$\\theta$ = %1.1f$^\\circ$' % ang[iang]) for iang in range(Nang)]
                
                # Make time independent labels
                for iang in range(Nang):
                    axs[iang].legend(loc=4)
                axs[-2].set_xlabel(xstr)
                axs[-1].set_xlabel(xstr)
                # common ylabel
                fig.text(0.04, 0.5, ystr, va='center',rotation='vertical')
                # initial title
                plt.suptitle('t = %1.1f [Myr]' % tarr[0], size='large')
                
                def animate(ifrm):
                    for iang in range(Nang):
                        lines[iang][0].set_ydata(vrot[ifrm,iang])
                    # update title 
                    plt.suptitle('t = %1.1f [Myr]' % tarr[ifrm], size='large')
                    return lines
                
                # do the animation
                ani = animation.FuncAnimation(fig, animate, range(n), repeat=False)
                if save:
                    print("[main]: Saving animation")
                    ani.save(quant+".mp4")
            else: 
                print("[main]: Not sure what to plot; must specify ifrm or anim. Type "
                            "'python plot_otf.py -h' for help.")
                quit()
        # handle plotting of single rotation curve 
        else: 
            # setup axes
            fig, ax = plt.subplots(figsize=(10,8))
            ax.set_xlabel(xstr) 
            ax.set_ylabel(ystr)

            # no animation
            if ifrm is not None:
                t = tarr[ifrm]
                ax.set_title('t = %1.1f [Myr]' % tarr[ifrm])
                ax.plot(r, vrot[ifrm, iang], color=c[iang], label='$\\theta$ = %1.1f$^\\circ$' % ang[iang])
                ax.legend(loc=4)
                if save:
                    print("[main]: Saving figure")
                    plt.savefig(quant+str(ifrm)+'_'+str(iang)+".eps")
            # handle animation
            elif anim: 
                line, = ax.plot(r, vrot[0,iang], color=c[iang], 
                        label='$\\theta$ = %1.1f$^\\circ$' % ang[iang])   
                ax.legend(loc=4)
                def animate(ifrm):
                    line.set_ydata(vrot[ifrm,iang])
                    ax.set_title('t = %1.1f [Myr]' % tarr[ifrm])
                    return line,
                ani = animation.FuncAnimation(fig, animate, range(n), repeat=False)
                if save:
                    print("[main]: Saving animation")
                    ani.save(quant+str(iang)+".mp4")
            else: 
                print("[main]: Not sure what to plot; must specify ifrm or anim. Type "
                            "'python plot_otf.py -h' for help.")
                quit()
    # show 
    if save:
        print("[main]: Program complete")
    else:
        plt.show()
    #=============================================================
                
main()
