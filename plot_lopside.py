#!/usr/bin/python
import numpy as np
import sys
import os
import subprocess as sbp
import matplotlib.pyplot as plt
import argparse
from argparse import RawTextHelpFormatter

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
from read_bin import read_otf as rotf
import read_athinput

#=====================================================
#
#  Code: plot_lopside.py
#
#  Purpose: Reads in otf files from hvc simulations
#           and does data analysis 
#
#  Keywords: python plot_lopside -h   
#
#  Usage: python plot_lopside.py quant (--nrand, --stat)  
#          
#  Author: John Dupuy 
#          UNC Chapel Hill
#  Date:    01/16/18 
#  Updated: 01/16/18 
#=====================================================


#=======FUNCTIONS=================

#=================================
# get_dirs() function
def get_dirs(cwd):
    # Purpose: get the directories of each simulation
    #          and store them in a list
    # input: 
    #         cwd: current working directory
    # output: 
    #         a list of strings that contain
    #         all the directories of simulations
    
    # Use the find function to get directory structure
    process = sbp.Popen(["find",cwd,"-type","d"], stdout=sbp.PIPE)
    stdout, stderr = process.communicate()

    # Create empty lists for each mass range
    m0   = []
    m1e5 = []
    m1e6 = []
    
    print('\n[get_dirs]: Finding simulation directories...')
    # Loop through stdout and get sims
    for d in stdout.splitlines():
        if 'fac' in d:
            if 'id' not in d:
                if 'm0' in d:
                    m0.append(d+'/')
                if 'm1e5' in d:
                    m1e5.append(d+'/')
                if 'm1e6' in d:
                    m1e6.append(d+'/')

    return m0, m1e5, m1e6

#==================================
# get_files() function
def get_files(dirs,proc,flag):
    # input: ALL TYPE str
    #       dirs: directory that contains athinput
    #       proc: processor no.
    #       flag: unique flag for your file 
    # output: a list of lists that for each simulation, contains 
    #         the paths of ever file containing 'flag' 
    # get the files you want to do analysis on 
    
    m = '???'
    if 'm0' in dirs[0]:
        m = 'm0'
    if 'm1e5' in dirs[0]:
        m = 'm1e5'
    if 'm1e6' in dirs[0]:
        m = 'm1e6'  

    print("[get_files]: Getting files containing '%s' for each %4s simulation..." % (flag,m) )  
    # Create empty list to store file names for 
    # each simulation
    files = []
    for d in dirs:
        try:
            files.append(sorted([d+proc+'/'+fname for
                          fname in os.listdir(d+proc+'/')
                          if flag in fname]))
        except OSError: 
            print('[get_files]: WARNING, sim files not present for %s ' % (d) )
    return files  

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



#==================================
# get_stats() function
def get_stats(files, quant, nrand=10, stat='rand_chc',**kwargs):
    # input: 
    #       files: output of get_files()
    #       quant: quantity of interest
    #       nrand: number of random samples
    #       
    # output: 
    #        numpy array(s) of desired data
    Nang       = 8
    nx1_dom    = 256
    ctrl_files = 1
    rmn = 10.
    rmx = 5000.
    # parse keyword args
    for key in kwargs:
        if key == 'Nang':
            Nang    = kwargs[key]
        if key == 'nx1_dom':
            nx1_dom  = kwargs[key]
        if key == 'ctrl_files':
            ctrl_files = kwargs[key]
            nctrl      = len(ctrl_files)
        if key == 'rmn':
            rmn = kwargs[key]
        if key == 'rmx':
            rmx = kwargs[key]

    # Total number of simulations
    nsim = len(files)
    
    # Create numpy array to store data
    mydata = np.zeros(nsim)

    # Begin main loop over simulations 
    for i in range(nsim):
        # Get length of each simulation
        n    = len(files[i])
        # time varying data array
        tdat = np.zeros(n) 
        for j in range(n):
            t, mhvc, rhvc, rpos,\
            acc_rate, facvhvc, ahvc,\
            mcR, mcL,\
            r, ang, vrot,\
            A1, A2         = get_data(files[i][j], **kwargs)

            if quant == 'A1':
                tdat[j] = np.mean(A1[(rmn < r) & (r < rmx)])
            elif quant == 'A2':
                tdat[j] = np.mean(A2[(rmn < r) & (r < rmx)]) 
            elif quant == 'LoR':
                tdat[j] = mcL/mcR
            elif quant == 'RoL':
                tdat[j] = mcR/mcL
            else: 
                print('[get_stats]: quant not understood, exiting...')
                quit()
                
        # Now do the statistical analysis
        if stat == 'rand_chc':   
            mydata[i] = np.mean(np.random.choice(tdat,nrand))
        elif stat == 't_mean':
            mydata[i] = np.mean(tdat)
        else: 
            print('[get_stats]: stat not understood, exiting...')
            quit() 

    return mydata

# main() function 
def main():
    # First find the directories of the simulations
    m0d, m1e5d, m1e6d = get_dirs(os.getcwd())
    # Set the processor and unique flag
    proc = 'id0'
    flag = 'otf'
    # Get global parameters from the m0 simulation
    athin = [fname for fname in os.listdir(m0d[0]) if 'athinput' in fname]
    prec  = 32
        # Read input file
    base, params = read_athinput.readath(m0d[0] + athin[0])
    nx1_dom = int(params[0].nx1)


    # Get all the paths for the simulation files     
    m0f   = get_files(m0d  ,proc,flag)
    m1e5f = get_files(m1e5d,proc,flag)
    m1e6f = get_files(m1e6d,proc,flag)

    # Now parse input 
        # Read in system arguments
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter) 
    parser.add_argument("quant",type=str,
                        help="Quantity options:\n"
                             "  RoL: mcR/mcL\n"
                             "  LoR: mcL/mcR\n"
                             "   A1: Lopsidedness parameter for m=1\n"
                             "   A2: Lopsidedness parameter for m=2")
    parser.add_argument("--stat", dest="stat",type=str, default='rand_chc',
                        required=False,
                        help="Statistics to carry out on data\n"
                             "  rand_chc: randomly selects nrand elements from t-dependent data\n"
                             "  t_mean: takes time average of t-dependent data")
    parser.add_argument("--nrand", default=10,required=False,
                        type=int,
                        help="number of points of quant(t) to randomly select\n")
    parser.add_argument("--rmnmx",type=float,dest='rmnmx',nargs=2,
                        required=False,default=(params[0].x1min,params[0].x1max),
                        help="Min/max radii for computing <A1> and <A2>\n")
    parser.add_argument("--save", dest="save",action='store_true',
                        default=False,
                        help="Switch to save figure")
    parser.add_argument("--cut", dest="cut",default=0.05,type=float,
                        help="Cutoff value for determining lopsidedness")

    # parsing arguments
    args = parser.parse_args()
    quant    = args.quant
    stat     = args.stat
    nrand    = args.nrand
    rmn, rmx = args.rmnmx
    save     = args.save
    cut      = args.cut

    # Now read in data
    dat0   = get_stats(m0f  , quant, nrand=nrand, stat=stat, nx1_dom=nx1_dom, rmn=rmn, rmx=rmx)  
    dat1e5 = get_stats(m1e5f, quant, nrand=nrand, stat=stat, nx1_dom=nx1_dom, rmn=rmn, rmx=rmx)
    dat1e6 = get_stats(m1e6f, quant, nrand=nrand, stat=stat, nx1_dom=nx1_dom, rmn=rmn, rmx=rmx)


    # print out simulation info for simulations that are above cutoff 
    m1e5d_above_cut = np.array(m1e5d)[dat1e5 > cut]  
    m1e6d_above_cut = np.array(m1e6d)[dat1e6 > cut] 

    print('[main]: m1e5, N_lop/N_sim = %3d/%3d' % (len(m1e5d_above_cut), len(m1e5d)) )
    print('[main]: m1e6, N_lop/N_sim = %3d/%3d' % (len(m1e6d_above_cut), len(m1e6d)) )
     
    plt.figure()
    plt.plot(len(dat1e5)/2, dat0,'bo',label='control sim')
    plt.plot(dat1e5,'r.',label='m1e5 sims') 
    plt.plot(dat1e6,'g.',label='m1e6 sims')
    plt.xlabel('Simulation number')
    plt.ylabel('Stat: %s of %s' % (stat,quant))
    plt.legend()
    if save:
        plt.savefig(quant+stat+'.eps') 
    else:
        plt.show()  


    return

main()
