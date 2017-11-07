#!/usr/bin/python
import numpy as np
import sys 
import os
import fnmatch as fnm
# matplotlib imports
from matplotlib import pyplot as plt
import matplotlib.animation as animation
# Import from correct directory
sys.path.insert(0,'/afs/cas.unc.edu/users/j/d/'
                  'jdupuy26/Johns_work/'
                  'analysis_files/plotting_tools')
import read_otf as rotf
import read_athinput
from units_class import * 


#=====================================================
#
#  Code: plot_otf.py
#
#  Purpose: Plot files read in from otf dumps 
#
#  Keywords: ???
#
#  Usage: ??? 
#          
#  Author: John Dupuy 
#          UNC Chapel Hill
#  Date:   11/7/17
#=====================================================

#============FUNCTIONS===========
                                   
def get_files(dirs,proc,flag):
    # input: ALL TYPE str
    #       dirs: directory that contains athinput
    #       proc: processor no.
    #       flag: unique flag for your file 
    # get the files you want to do analysis on 
    return [dirs+proc+'/'+fname for 
            fname in os.listdir(dirs+proc+'/')
            if fnm.fnmatch(fname,flag)] 

def get_data(file,prec,**kwargs):
    # Set units 
    iunit = 0
    for key in kwargs: 
        if key == 'iunit':
            iunit = kwargs[key] 
    
    t, mhvc, rhvc, rpos,\
    acc_rate, mcR, mcL,\
    r, ang, vrot = rotf.readbin(file,prec,**kwargs)

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
           acc_rate, mcR, mcL,\
           r, ang, vrot 

def main():
    # Get input file
    athdir = os.getcwd()+'/'
    athin = [fname for fname in os.listdir(athdir) if fname.startswith("athinput")][0]
    # Read input file
    base, params = read_athinput.readath(athdir+athin) 
    # Set precision
    prec  = 32

    # Get otf files 
    otf_files = get_files(athdir,'id0','*.otf.*')

    # Create arrays to store t-dependent data
    n = len(otf_files)
    tarr     = np.zeros(n)
    mcR_arr  = np.zeros(n)
    mcL_arr  = np.zeros(n)
    vrot_arr = np.zeros((n,8,256))

    # Read in data
    for i in range(n):
        t, mhvc, rhvc, rpos,\
        acc_rate, mcR, mcL,\
        r, ang, vrot        = get_data(otf_files[i],prec)
        
        tarr[i]     = t
        mcR_arr[i]  = mcR
        mcL_arr[i]  = mcL
        vrot_arr[i] = vrot 
   
    # convert angles to deg
    ang *= 180./np.pi

    iang = 7
    mycol = ['b.','g.','r.','c.','m.','y.','k.','b-']
    
    # Animate the rotation curve
    fig, ax = plt.subplots()
    line, = ax.plot(r,vrot_arr[0,0],mycol[iang], 
                        label='$\\theta$ = %1.2f$^\\circ$' % ang[iang])
    ax.set_xlabel('r [pc]')
    ax.set_ylabel('v [pc/Myr]') 
    ax.legend(loc=4)
    def animate(i):
        line.set_ydata(vrot_arr[i,iang]) 
        ax.set_title('t = %1.3f [Myr]' % tarr[i])
        return line,
    ani = animation.FuncAnimation(fig, animate, range(n))
    ani.save('my_rot.mp4')
    
    #plt.show() 

main()

