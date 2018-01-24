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

#=======DICTIONARIES==============

#=================================
# sym dictionary
sym = {"m1e6":"M$_c$=10$^6$ M$_{\odot}$", "m1e5": "M$_c$=10$^5$ M$_{\odot}$",
       "a":"$\phi_{c,0}$=", "fac":"f$_{c}$=", "te":"t$_e$=",
       "r":"r$_{pos,0}$=", "rc":"r$_c$=", "f":"$f$", "A1":"<A1>",
       "A2":"<A2>", "RoL": "M$_{cr,R}$/M$_{cr,L}$", "LoR":"M$_{cr,L}$/M$_{cr,R}$"} 
 

#=======FUNCTIONS=================

#=================================
# get_dirs() function
def get_dirs(cwd, inr3000):
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
        if inr3000:
            if 'fac' in d:
                if 'id' not in d:
                    if 'm0' in d:
                        m0.append(d+'/')
                    if 'm1e5' in d:
                        m1e5.append(d+'/')
                    if 'm1e6' in d:
                        m1e6.append(d+'/')
        else:
            if 'r3000' not in d:
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
    for i in range(len(dirs)):
        try:
            sfiles = sorted([dirs[i]+proc+'/'+fname for
                          fname in os.listdir(dirs[i]+proc+'/')
                          if flag in fname]) 
            if len(sfiles) != 51:
                print("[get_files]: WARNING, exluding sim that did not complete: %s" %dirs[i])
                dirs[i] = 'removed'  # remove this simulation from the directory list
            else: 
                files.append(sfiles)  

        except OSError: 
            print('[get_files]: WARNING, sim files not present for %s ' % (dirs[i]) )
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
# get_mask() function
def get_mask(dirs, flag):
    # Given a unique flag, returns a boolean
    # array, used to plot simulations as 
    # functions of mass parameters
    return np.array([flag in x for x in dirs])


#==================================
# get_sortdat() function
def get_sortdat(dat, err, dirs, param,inr3000=False):
    # Given a certain parameter
    # Separates data according to that 
    # parameter

    dats = []
    errs = []
    
    # Decide what param we are plotting
    if inr3000: 
        select = {'rc':[200,400,600],'te':[265,275,285],
                  'r':[500,1000,2000,3000],'fac':[0.0,1.0,2.0],
                  'a':[0.0,1.5708]}
    else: 
        select = {'rc':[200,400,600],'te':[265,275,285],
                  'r':[500,1000,2000],'fac':[0.0,1.0,2.0],
                  'a':[0.0,1.5708]}
    try:
        params = select[param]
    except KeyError:
        print('[get_sortdat]: ERROR, param %s not understood, exiting...\n' % param)
        quit()

    # reconstruct dirs so that it gets rid of removed simulations/entries
    dirs = [x for x in dirs if x != 'removed']

    for p in params:
        dats.append(dat[get_mask(dirs,param+str(p))])
        errs.append(err[get_mask(dirs,param+str(p))])
    
    return params, dats, errs  

#==================================
# get_cmf() function
def get_cmf(dat, err, x): 
    # Given data, and a x-axis value, x,
    # the cumulative function represents
    # the fraction of data that has value 
    # that is greater than x.
    # Also computes the error based off,
    # low and high estimates.   
    
    f   = np.zeros(len(x))
    flo = np.zeros(len(x))
    fhi = np.zeros(len(x))

    for i in range(0,len(x)):
        flo[i] = len(dat[(dat - err) > x[i]])
        f[i]   = len(dat[dat > x[i]])
        fhi[i] = len(dat[(dat + err) > x[i]])
    
    flo /= len(dat)
    f   /= len(dat)
    fhi /= len(dat)
    lo_e = f - flo
    hi_e = fhi - f

    return f, (lo_e, hi_e)  

        

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
        if key == 'ntrials':
            ntrials = kwargs[key]

    # Total number of simulations
    nsim = len(files)
    
    # Create numpy array to store data
    mydata = np.zeros((ntrials,nsim))

    # Begin main loop over simulations 
    for i in range(nsim):
        # Get length of each simulation
        n    = len(files[i])
        
        # time varying data array
        tdat = np.zeros(n) 
        
        # start time loop
        for j in range(n):
            t, mhvc, rhvc, rpos,\
            acc_rate, facvhvc, ahvc,\
            mcR, mcL,\
            r, ang, vrot,\
            A1, A2         = get_data(files[i][j], **kwargs)

            t, mhvc, rhvc, rpos,\
            acc_rate, facvhvc, ahvc,\
            mcR, mcL,\
            r, ang, vrot,\
            A1c, A2c       = get_data(ctrl_files[j], **kwargs)

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

        for k in range(ntrials):
            if stat == 'rand_chc':   
                mydata[k,i] = np.mean(np.random.choice(tdat,nrand))
            elif stat == 't_mean':
                mydata[k,i] = np.mean(tdat)
            elif stat == 'max':
                mydata[k,i] = np.max(tdat)
            else: 
                print('[get_stats]: stat not understood, exiting...')
                quit() 
    # Now take get the mean and std dev from mydata
    mean = np.mean(mydata,axis=0)
    err  = np.sqrt(np.var(mydata,axis=0))

    return mean, err


# main() function 
def main():
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
                        required=False,default=(10.,5000.),
                        help="Min/max radii for computing <A1> and <A2>\n")
    parser.add_argument("--save", dest="save",action='store_true',
                        default=False,
                        help="Switch to save figure")
    parser.add_argument("--cut", dest="cut",default=0.05,type=float,
                        help="Cutoff value for determining lopsidedness")
    parser.add_argument("--param", dest='param',type=str, default='number',
                        required=False,
                        help='Parameter to plot on x-axis')
    parser.add_argument("--cmf", dest='cmf',action='store_true',
                        default=False,
                        help="Switch to get the cumulative function f")
    parser.add_argument("--ntrials",type=int,dest='ntrials',required=False,
                        default=1,
                        help="Number of times to do the random choice with nrand."
                             "If ntrials>1, then errorbars on <A1> represent the"
                             "std dev of the measurement, and the measurement itself"
                             "is the mean.\n") 
    parser.add_argument("--inr3000", dest='inr3000',action='store_true',
                        default=False,
                        help="Switch to include r3000 sims, default is False.")

    # parsing arguments
    args = parser.parse_args()
    quant    = args.quant
    stat     = args.stat
    nrand    = args.nrand
    rmn, rmx = args.rmnmx
    save     = args.save
    cut      = args.cut
    param    = args.param
    cmf      = args.cmf
    ntrials  = args.ntrials
    inr3000  = args.inr3000


    # First find the directories of the simulations
    m0d, m1e5d, m1e6d = get_dirs(os.getcwd(),inr3000)
    
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

    # error checking for conflicting arguments
    if cmf and param == 'number':
        print('[main]: Cumulative function cannot be plotted with number!')
        quit()

    if stat != 'rand_chc' and ntrials != 1:
        print('[main]: Using more than one trial does not make sense for stat: %s!' %(stat))
        quit()

    # Now read in data
    print("[get_stats]: Starting main loop to fill data array.\n")

    dat0  , err0   = get_stats(m0f  , quant, nrand=nrand, stat=stat, nx1_dom=nx1_dom, 
                                     rmn=rmn, rmx=rmx, ntrials=ntrials)  
    dat1e5, err1e5 = get_stats(m1e5f, quant, nrand=nrand, stat=stat, nx1_dom=nx1_dom, 
                                     rmn=rmn, rmx=rmx, ntrials=ntrials)
    dat1e6, err1e6 = get_stats(m1e6f, quant, nrand=nrand, stat=stat, nx1_dom=nx1_dom, 
                                     rmn=rmn, rmx=rmx, ntrials=ntrials)


    # print out simulation info for simulations that are above cutoff 
    '''
    m1e5d_above_cut = np.array(m1e5d)[dat1e5 > cut]  
    m1e6d_above_cut = np.array(m1e6d)[dat1e6 > cut] 

    print('[main]: m1e5, N_lop/N_sim = %3d/%3d' % (len(m1e5d_above_cut), len(m1e5d)) )
    print('[main]: m1e6, N_lop/N_sim = %3d/%3d' % (len(m1e6d_above_cut), len(m1e6d)) )
    '''


    print('[main]: Max error in 1e5, %e' %(np.max(err1e5)) ) 
    print('[main]: Max error in 1e6, %e\n' %(np.max(err1e6)) ) 

    if param == 'number':
        xv1e5 = np.arange(len(dat1e5))
        xv1e6 = np.arange(len(dat1e6))
        
        plt.figure()
        plt.errorbar(len(dat1e5)/2,dat0,yerr=err0,fmt='bo',label='control sim')
        plt.errorbar(xv1e5,dat1e5,yerr=err1e5,fmt='r.',label=sym['m1e5']) 
        plt.errorbar(xv1e6,dat1e6,yerr=err1e6,fmt='g.',label=sym['m1e6'])
        plt.xlabel('Simulation number')
        plt.ylabel('Stat: %s of %s' % (stat,quant))
        plt.legend()
        if save:
            plt.savefig(quant+stat+'.eps') 
        else:
            plt.show()  
    else:
        params, dats1e5, errs1e5 = get_sortdat(dat1e5, err1e5, m1e5d, param,inr3000)
        params, dats1e6, errs1e6 = get_sortdat(dat1e6, err1e6, m1e6d, param,inr3000)
        
        plt.figure()
        # General plotting variables 
        lines = ['-','--','-.',':']
        # Define cmf variables 
        if cmf: 
            amin  = 0.8*np.min(dat1e5) 
            amax  = 1.2*np.max(dat1e6)
            da    = 5*(amax-amin)/len(dat1e6)
            
        # Plot m1e6
        for (p,d,e,l) in zip(params,dats1e6,errs1e6,lines):
            if cmf:
                avals = np.arange(amin,amax,da)
                fvals, er = get_cmf(d,e,avals)
                plt.plot(avals, fvals, 'b'+l+'x',
                         label=sym['m1e6']+', '+sym[param]+str(p))
            else:
                plt.plot(p*np.ones(len(d)), d, 'b.',label=sym['m1e6'])
        
        # Plot m1e5 
        for (p,d,e,l) in zip(params,dats1e5,errs1e5,lines):
            if cmf:
                avals = np.arange(amin,amax,da)
                fvals, er = get_cmf(d,e,avals)
                plt.plot(avals, fvals, 'g'+l+'+',
                         label=sym['m1e5']+', '+sym[param]+str(p))
            else:
                plt.plot(p*np.ones(len(d)), d, 'g+',label=sym['m1e5'])
        
        if cmf:
            plt.xlim(avals[0], avals[-1])
            plt.ylim(-0.1,1.1)
            plt.xlabel(sym[quant])
            plt.ylabel('Cumulative function, '+sym['f'])
        else:
            # Get spacing for plotting stuff
            dxp = (params[-1] - params[0])/len(params)
            
            plt.plot(0.5*dxp, dat0, 'bo', label='Control sim')
            plt.xlim(params[0]-0.1*dxp, params[-1]+0.1*dxp)
            plt.xlabel('Parameter: %s' %(sym[param]) ) 
            plt.ylabel('Stat: %s of %s' % (stat,sym[quant]))
        plt.legend() 
        if save:
            plt.savefig(quant+stat+'.eps')
        else:
            plt.show() 
        
    return

main()
