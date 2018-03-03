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
from mpl_toolkits.mplot3d import Axes3D 
# scipy imports 
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from scipy.signal import argrelextrema
import scipy.interpolate as itp

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

from read_bin import read_otf as rotf
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
#  Usage: python plot_otf.py quant (--ifrm, --anim, --iani) --iang  
#          
#  Author: John Dupuy 
#          UNC Chapel Hill
#  Date:    11/07/17
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
    
    t, mhvc, rhvc, rpos,\
    acc_rate, facvhvc, ahvc,\
    mcR, mcL,\
    r, ang, vrot,\
    A1, A2           = rotf(file,prec, **kwargs)

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
           r, ang, vrot,\
           A1, A2 

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
    elif quant == 'vrot': 
        lab = 'v [pc/Myr]'
    elif quant == 'A1':
        lab = 'A$_1$ [unitless]'
    elif quant == 'A2':
        lab = 'A$_2$ [unitless]'
    elif quant == '<A1>':
        lab = '<A$_1$> [unitless]'
    elif quant == '<A2>':
        lab = '<A$_2$> [unitless]'
    return lab

#================================
# get_title() function
def get_title(quant,comp):
    lab = {'A1':{'None':'A1 [unitless]',
                 'sub' :'$|$A1 - A1$_{control}$$|$ [unitless]',
                 'div' :'A1/A1$_{control}$ [unitless]'},
           'A2':{'None':'A2 [unitless]',
                 'sub' :'$|$A2 - A2$_{control}$$|$ [unitless]',
                 'div' :'A2/A2$_{control}$ [unitless]'}
          }
    return lab[quant][comp]
    

#================================
# exp_func() function
def exp_func(x,a,b,c=0):
    return a*np.exp(b*x) + c
        
#===============================
def lin_func(x,m,b):
    return m*x + b

#===============================
def exp_fit(x,dat, c=0):
    dat -= c
    dat = np.log(dat)
    b, log_a = np.polyfit(x, dat, 1)
    a = np.exp(log_a)
    return a, b

#==============================
def lin_fit(x,dat):
    m, b = np.polyfit(x, dat, 1)
    return m, b
    
#================================
# initialize() function
def init(quant, files, **kwargs):
    Nang       = 8
    nx1_dom    = 256
    # parse keyword args
    for key in kwargs:
        if key == 'Nang':
            Nang    = kwargs[key]
        if key == 'nx1_dom':
            nx1_dom = kwargs[key]
        if key == 'rmnmx':
            rmnmx = kwargs[key]

    # set length of simulation from file list
    n    = len(files)
    # setup time array
    tarr = np.zeros(n)
    # Parse rmnmx
    rmn = rmnmx[0]
    rmx = rmnmx[1]

    # choose quant
    if (quant == 'mcR' or quant == 'mcL' or 
        quant == 'RoL' or quant == 'LoR'):
        mcR_arr = np.zeros(n)
        mcL_arr = np.zeros(n)
        
        # Read in data
        for i in range(n):
            t, mhvc, rhvc, rpos,\
            acc_rate, facvhvc, ahvc,\
            mcR, mcL,\
            r, ang, vrot,\
            A1, A2         = get_data(files[i], **kwargs)
            
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
        # Set data
        if quant == "mcR":
            data = mcR_arr
        elif quant == "mcL":
            data = mcL_arr
        elif quant == "LoR":
            data = mcL_arr/mcR_arr
        else:
            data = mcR_arr/mcL_arr

    elif quant == 'vrot': 
        vrot_arr = np.zeros((n,Nang,nx1_dom))
        for i in range(n):
            t, mhvc, rhvc, rpos,\
            acc_rate, facvhvc, ahvc,\
            mcR, mcL,\
            r, ang, vrot,\
            A1, A2        = get_data(files[i], **kwargs)
            
            tarr[i]     = t
            vrot_arr[i] = vrot
        # print HVC parameters 
        print("\n[init]: mhvc = %1.3e [M_sun]\n"
              "        rhvc = %1.3e [pc]\n"
              "        rpos = %1.3e [pc]\n"
              "    acc_rate = %1.3e [M_sun/Myr]\n"
              "     facvhvc = %1.3e [unitless] \n"
              "        ahvc = %1.3e [rad]" % (mhvc, rhvc, rpos, acc_rate,facvhvc,ahvc))
        data = (r, ang, vrot_arr)
    
    elif (quant == 'A1' or quant == 'A2' or
          quant == '<A1>' or quant == '<A2>'): 
        A1 = np.zeros((n,nx1_dom)) 
        A2 = np.zeros((n,nx1_dom))
        # Read in data
        for i in range(n):
            t, mhvc, rhvc, rpos,\
            acc_rate, facvhvc, ahvc,\
            mcR, mcL,\
            r, ang, vrot,\
            a1, a2         = get_data(files[i], **kwargs)
            
            tarr[i] = t
            A1[i]   = a1
            A2[i]   = a2  

        # print HVC parameters 
        print("\n[init]: mhvc = %1.3e [M_sun]\n"
              "        rhvc = %1.3e [pc]\n"
              "        rpos = %1.3e [pc]\n"
              "    acc_rate = %1.3e [M_sun/Myr]\n"
              "     facvhvc = %1.3e [unitless] \n"
              "        ahvc = %1.3e [rad]" % (mhvc, rhvc, rpos, acc_rate,facvhvc,ahvc))
        # Set data
        if   quant == 'A1':
            data = [r, A1]
        elif quant == 'A2':
            data = [r, A2]    
        elif quant == '<A1>':
            # Get r values of peaks in bar zone
            rbar_peaks  = r[np.argmax(A1[:,(0 < r) & ( r < 1200)],axis=1)]
            # Get r values of peaks in disc zone
            rdisc_peaks = r[np.argmax(A1[:,(rmn < r) & (r < rmx)],axis=1) + 197]

            diff = rdisc_peaks - rbar_peaks
            # Average 
            avg = np.mean(A1[:, (rmn<r) & (r < rmx)],axis=1)
            try:
                mxi = int(tarr[ diff > np.max(diff)/4.][0] - 250.)
            except IndexError:
                mxi = np.argmax(avg)  

            data = (avg, mxi)
        elif quant == '<A2>':
            data = np.mean(A2[:,(rmn < r) & (r < rmx)],axis=1)
    
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
    # Sort them
    otf_files.sort()
    
    n   = len(otf_files)

    # Read in system arguments
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter) 
    parser.add_argument("quant",type=str,
                        help="Plotting options:\n"
                             "  mcR: mass w/in circleR\n"
                             "  mcL: mass w/in circleL\n"
                             "  RoL: mcR/mcL\n"
                             "  LoR: mcL/mcR\n"
                             " vrot: rotation curve\n"
                             "   A1: Lopsidedness parameter for m=1\n"
                             "   A2: Lopsidedness parameter for m=2\n"
                             " <A1>: spatial average of A1\n"
                             " <A2>: spatial average of A2\n")
    parser.add_argument("--iang", dest="iang",type=int, default=0,
                        required=False,
                        help="Angle to plot if using vrot:\n"
                             "  [0-(Nang-1)]: index that corresponds to [0-305 deg]\n"
                             " -1: All angles are plotted\n")
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
    parser.add_argument("--rmnmx", dest="rmnmx",nargs=2,required=False,
                        default=[1200., 4000.],type=float,
                        help="Average A1 from rmn to rmx") 
    parser.add_argument("--save", dest="save",action='store_true',
                        default=False,
                        help="Switch to save anim or figure")
    parser.add_argument("--cmap", dest="cmap",action='store_true',
                        default=False,
                        help="Switch to create a colormap for A1 & A2")
    parser.add_argument("--wire", dest='wire',action='store_true',
                        default=False,
                        help="Switch to make a 3d wireframe plot")
    parser.add_argument("--fit", dest="fit",type=str,
                        default='None',
                        help="Switch to do fit for <A1> perturbations")
    parser.add_argument("--comp", dest="comp",type=str,
                        default='None',
                        help="String to compare measurement to control files:\n"
                             "  sub: Subtract control from sim and take abs value\n"
                             "  div: Divide sim by control")
    parser.add_argument("--interp", dest="interp",action='store_true',
                        default=False,
                        help="Switch to do interpolation on <A1>")
    parser.add_argument("--fmt", dest="fmt",type=str,default='eps',
                        help="Type of figure to save, only applies to single frames, default: eps") 

    # parsing arguments            
    args  = parser.parse_args()
    quant = args.quant
    iang  = args.iang
    anim  = args.anim
    iani  = args.iani
    ifrm  = args.ifrm
    save  = args.save
    rmnmx = args.rmnmx
    cmap  = args.cmap
    wire  = args.wire
    fit   = args.fit
    comp  = args.comp
    interp= args.interp
    fmt   = args.fmt

    if   fit == 'exp':
        fit_func = exp_fit
        m_func   = exp_func
    elif fit == 'lin':
        fit_func = lin_fit
        m_func   = lin_func

    # Conflicting argument checks 
    if iang > Nang:
        print("[main]: iang must be less than Nang!")
        quit()
    if anim and (quant == 'mcR' or quant == 'mcL' or 
                 quant == 'LoR' or quant == 'RoL' or
                 quant == '<A1>' or quant == '<A2>'):
        print("[main]: unable to animate these parameters")
        quit()
    if anim and ifrm != None:
        print("[main]: specifying a single frame while animating doesn't make sense")
        quit()
    if cmap and (quant != 'A1' and quant != 'A2'):
        print("[main]: colormaps must be plotted using quant A1 or A2")
        quit()
    if ifrm > (n-1):
        print("[main]: frame out of range of simulation!")
        quit()

    # Get the data
    if quant == '<A1>':
        tarr, (data, mxi) = init(quant, otf_files, prec=32,
                               nx1_dom=nx1_dom, Nang=Nang, rmnmx = rmnmx)
    else:
        tarr, data = init(quant, otf_files, prec=32, 
                          nx1_dom=nx1_dom,Nang=Nang,rmnmx=rmnmx) 
    
    print("       thvcs = %1.3e [Myr]\n"
          "       thvce = %1.3e [Myr]\n" % (params[0].thvcs, params[0].thvce) ) 

    # Check if comparing to control simulation
    if comp != 'None':
        ctrl_path  = "/srv/analysis/jdupuy26/longevity_study/"+\
                     "m0.0/te265/rc400/r1100/a0.0/fac1.0/"
        ctrl_files = get_files(ctrl_path,'id0','*.otf.*')
        ctrl_files.sort()
        nctrl      = len(ctrl_files)
        if nctrl != n:
            print("\n[main]: WARNING - no. of control files != no. of sim files\n"
                  "          nctrl = %i, n = %i" %(nctrl, n))
        print("\n[main]: Using control files located at:\n %s\n" % ctrl_path)
        
        # Get the control data
        tarr, dcon = init(quant, ctrl_files, prec=32,
                          nx1_dom=nx1_dom,Nang=Nang,rmnmx=rmnmx)
        
        
        # Print out some useful information
        # 1) MSE, mean squared error
        mse  = np.sum((data[1] - dcon[1])**2.)
        mse /= float(data[1].shape[0] * data[1].shape[1])

        print("\n[main]: MSE between control and sim: %1.3e \n" % (mse) ) 
 
        # Decide what to do with the control data
        if comp == 'sub':
            data[1] -= dcon[1]
            data[1]  = abs(data[1])
        elif comp == 'div':
            data[1] /= dcon[1]
        else:
            print('[main]: Comparison %s not recognized, exiting!' %(comp))
    
    
    if anim:
        print("[main]: Animating from t = %1.1f [Myr]\n" 
              "                    to t = %1.1f [Myr]\n"
                    %( tarr[iani[0]], tarr[iani[1]] ) )  

    if (quant == '<A1>' or quant == '<A2>'):
        print("[main]: Mean value of " + quant +
              " between %1.0f [pc] and %1.0f [pc] "
              "for whole simulation: %1.4f" 
              " \n" % (rmnmx[0],rmnmx[1],np.mean(data)) ) 

    #===== PLOTTING =============================================
    # set colors for lines
    colormap = cm.gist_ncar
    if Nang != 4:
        c    = [colormap(i) for i in np.linspace(0,0.4,Nang)]
    else: c = ['k','b','g','r'] 
    
    # set ylabel
    ystr = get_ylabel(quant)
    # for the mass/average A{1,2} stuff
    if (quant == 'mcR' or quant == 'mcL' or 
        quant == 'LoR' or quant == 'RoL' or
        quant == '<A1>' or quant == '<A2>'):

        # Define cutoff
        cut = 0.05 

        # Do plotting
        fig = plt.figure(figsize=(7.0,5.0),facecolor='white')
        plt.plot(tarr, data,'b-',marker='.')
        plt.plot(tarr, np.ones(len(data))*cut, 'r--') 
        plt.xlabel('t [Myr]')
        plt.ylabel(ystr) 

        # Get time above cutoff <A1> 
        try:
            ts = tarr[data > cut][-1] - tarr[data > cut][0]
        except IndexError:
            ts = 0.0 
        # Put on figure
        plt.text(0.7,0.82,'$\\tau$ = %1.1f [Myr]' %(ts), 
                 transform=fig.transFigure) 
        
        # Fit perturbation to exponential
        if quant == '<A1>' and fit != 'None':
            
            popt = fit_func(tarr[mxi:], data[mxi:]) 

            mx = np.max(data)
            ts_half = tarr[data > mx/2.0][-1] - tarr[data > mx/2.0][0]

            plt.plot(tarr,m_func(tarr,*popt),'k--',label=fit+' fit')
            print('[main]: Half life for <A1> perturbation: %1.2f' %( -1.0/popt[1] * np.log(2) ) )
            print('[main]: Time scale based on     1/2 max: %1.2f' %( ts_half ) )

        if quant == '<A1>' and interp:
            f = UnivariateSpline(tarr,data,k=2,s=1.0)
            #f = interp1d(tarr[::2], data[::2], kind='linear',fill_value='extrapolate')
            tnew = np.arange(0,1000,1) 
            tcut = tnew[np.where(f(tnew) < data[0])]
            
            if len(tcut)==0:
                f = UnivariateSpline(tarr,data,k=1,s=1.0)
                tcut = tnew[np.where(f(tnew) < data[0])] 
            #print(tcut[0]) 
            plt.plot(tarr,f(tarr),'r--') 
            plt.plot(tarr,data[0]*np.ones(len(data)),'k-')
       # if save:
       #     print("[main]: Saving figure")
       #     plt.savefig(quant+".eps")

    # for rotation curves
    elif (quant == 'vrot'):
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
                for iang in range(Nang):
                    axs[iang].plot(r, vrot[ifrm,iang],color=c[iang],
                                  label='$\\theta$ = %1.1f$^\\circ$' % ang[iang])
                    axs[iang].legend(loc=4,prop={'size': 10})
                    # Set the plotting range 
                    axs[iang].set_ylim(0.9*np.min(vrot[iani[0]:iani[1],:]),
                                       1.1*np.max(vrot[iani[0]:iani[1],:]) )

                # make labels 
                # xlabels
                axs[-2].set_xlabel(xstr)
                axs[-1].set_xlabel(xstr)
                # common ylabel
                fig.text(0.04, 0.5, ystr, va='center',rotation='vertical')
                # title
                plt.suptitle('t = %1.1f [Myr]' % tarr[ifrm], size='large')
               # if save:
               #     print("[main]: Saving figure")
               #     plt.savefig(quant+str(ifrm)+".eps")
            
            # handle animation   
            elif anim:
                # create list of lines 
                lines = [axs[iang].plot(r, vrot[0,iang],color=c[iang],
                            label='$\\theta$ = %1.1f$^\\circ$' % ang[iang]) for iang in range(Nang)]
                
                # Make time independent labels
                for iang in range(Nang):
                    axs[iang].legend(loc=2, prop={'size': 10})
                    # Set the plotting range 
                    axs[iang].set_ylim(0.9*np.min(vrot[iani[0]:iani[1],:]),
                                       1.1*np.max(vrot[iani[0]:iani[1],:]) )

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
                ani = animation.FuncAnimation(fig, animate, range(iani[0],iani[1]+1), repeat=False)
               # if save:
               #     print("[main]: Saving animation")
               #     ani.save(quant+".mp4")
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
                #if save:
                #    print("[main]: Saving figure")
                #    plt.savefig(quant+str(ifrm)+'_'+str(iang)+".eps")
            # handle animation
            elif anim: 
                line, = ax.plot(r, vrot[0,iang], color=c[iang], 
                        label='$\\theta$ = %1.1f$^\\circ$' % ang[iang])   
                ax.legend(loc=4)
                # Set the plotting range for the whole animation 
                ax.set_ylim(0.9*np.min(vrot[iani[0]:iani[1],iang]),
                            1.1*np.max(vrot[iani[0]:iani[1],iang]) )
                def animate(ifrm):
                    line.set_ydata(vrot[ifrm,iang])
                    ax.set_title('t = %1.1f [Myr]' % tarr[ifrm])
                    return line,
                ani = animation.FuncAnimation(fig, animate, range(iani[0],iani[1]+1), repeat=False)
                #if save:
                #    print("[main]: Saving animation")
                #    ani.save(quant+str(iang)+".mp4")
            else: 
                print("[main]: Not sure what to plot; must specify ifrm or anim. Type "
                            "'python plot_otf.py -h' for help.")
                quit()
    # Handle plotting of lopsidedness parameters
    elif quant == 'A1' or quant == 'A2':
        xstr = 'r [pc]'
        # setup axes
        fig, ax = plt.subplots(figsize=(10,8),facecolor='white')
        
        # unpack data
        r, A = data[0], data[1] 

                # handle plotting a colormap
        if (cmap) or (wire):
                        # 'Cell centered' time 
            dt = tarr[1] - tarr[0]
            tc = np.linspace(tarr[0]-0.5*dt, tarr[-1]+0.5*dt,len(tarr)) 
            # 'Cell centered' radii
            lr = np.log(r)
            dx = lr[1] - lr[0]
            rc = np.exp(np.linspace(lr[0]-0.5*dx, lr[-1]+0.5*dx,len(r)))

            # Create the meshgrid
            R, T = np.meshgrid(rc,tc)
            
            if cmap:
                # Make the colormap
                im = ax.pcolormesh(R,T,A,cmap='magma')
                cbar = plt.colorbar(im,pad=0.1,label=get_title(quant,comp)) 
            else: 
                ax  = fig.add_subplot(111,projection='3d')
                surf = ax.plot_surface(R,T,A,rstride=3,cstride=3,cmap='rainbow',shade=False)
                ax.set_zlabel(get_title(quant,comp))
            
            
            ax.set_ylabel('t [Myr]')
            ax.set_xlabel('r [pc]')
            ax.set_xlim(rc[0],rc[-1])
            ax.set_ylim(tc[0],tc[-1])

        else:
            ax.set_xlabel(xstr) 
            ax.set_ylabel(ystr)
            
            # no animation
            if ifrm is not None:
                t = tarr[ifrm]
                ax.set_title('t = %1.1f [Myr]' % tarr[ifrm])
                ax.plot(r, A[ifrm])
                #if save:
                #    print("[main]: Saving figure")
                #    plt.savefig(quant+str(ifrm)+'_'+str(iang)+".eps")
            # handle animation
            elif anim: 
                line, = ax.plot(r, A[0]) 
                # Set the plotting range for the whole animation 
                ax.set_ylim(0.9*np.min(A[iani[0]:iani[1]]),
                            1.1*np.max(A[iani[0]:iani[1]]))
                def animate(ifrm):
                    line.set_ydata(A[ifrm])
                    ax.set_title('t = %1.1f [Myr]' % tarr[ifrm])
                    return line,
                ani = animation.FuncAnimation(fig, animate, range(iani[0],iani[1]+1), repeat=False)
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
        mydir  = '/srv/analysis/jdupuy26/figures/'
        myname = os.getcwd().split('longevity_study/',1)[1].replace('/','_')
        plt.savefig(mydir+myname+'_'+base+'_'+quant+'.'+fmt,format=fmt,bbox_inches='tight')
    else:
        plt.show()
    #=============================================================
                
main()
