#!/usr/bin/python
import numpy as np
import sys
import os
import subprocess as sbp
# matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D 
# scipy
from scipy.interpolate import UnivariateSpline as uspl
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
from read_bin import read_otf as rotf
import read_athinput

#=====================================================
#
#  Code: plot_mse.py
#
#  Purpose: Reads in otf files from hvc simulations
#           and does computes the mse error between 
#           the simulation and the control 
#
#  Keywords: python plot_lopside -h   
#
#  Usage: python plot_lopside.py quant   
#          
#  Author: John Dupuy 
#          UNC Chapel Hill
#  Date:    01/24/18 
#  Updated: 02/15/18 
#=====================================================

#=======DICTIONARIES==============

#=================================
# sym dictionary
sym = {"m1e6":"M$_c$=10$^6$ M$_{\odot}$", "m1e5": "M$_c$=10$^5$ M$_{\odot}$",
       "m3e6":"M$_c$=3 $\\times$ 10$^6$ M$_{\odot}$", 
       "m6e6":"M$_c$=6 $\\times$ 10$^6$ M$_{\odot}$",
       "m5e6":"M$_c$=5 $\\times$ 10$^6$ M$_{\odot}$",
       "m1e7":"M$_c$=10$^7$ M$_{\odot}$",
       "m5e7":"M$_c$=5 $\\times$ 10$^7$ M$_{\odot}$",
       "m1e8":"M$_c$=10$^8$ M$_{\odot}$",  
       "m0"  :"M$_c$=0 M$_{\odot}$",
       "a":"$\phi_{c,0}$", "fac":"f$_{c}$", "te":"t$_e$",
       "r":"r$_{pos,0}$", "rc":"r$_c$", "f":"$f$", "A1":"<A1>",
       "A2":"<A2>", "RoL": "M$_{cr,R}$/M$_{cr,L}$", "LoR":"M$_{cr,L}$/M$_{cr,R}$"} 

# units dictionary
units = {"a":" [rad]","rc":" [pc]","fac":" [unitless]","r":" [pc]",
         "te":" [Myr]"} 
 

#=======FUNCTIONS=================

#=================================
# get_dirs() function
def get_dirs(cwd, inr3000,inm1e8):
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

    # find mass directories 
    masses = os.walk('.').next()[1]
    if not inm1e8:
        masses.remove('m1e8') 
    # create empty dictionary of emtpy lists
        # for each mass  
    mdict  = {m: [] for m in masses}
    
    print('\n[get_dirs]: Finding simulation directories...')
    # Loop through stdout and get sims
    for d in stdout.splitlines():
        if inr3000:
            if ('fac' in d) and ('id' not in d):
                for m in masses:
                    if m in d:
                        mdict[m].append(d+'/') 
        else:
            if ('r3000' not in d) and ('fac' in d) and ('id' not in d):
                for m in masses:
                    if m in d:
                        mdict[m].append(d+'/')
    return mdict 

#==================================
# get_files() function
def get_files(dirs,proc,flag, num_files):
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
    if 'm3e6' in dirs[0]:
        m = 'm3e6'
    if 'm5e6' in dirs[0]:
        m = 'm5e6'
    if 'm6e6' in dirs[0]:
        m = 'm6e6' 
    if 'm1e7' in dirs[0]:
        m = 'm1e7'
    if 'm5e7' in dirs[0]:
        m = 'm5e7'
    if 'm1e8' in dirs[0]:
        m = 'm1e8'

    print("[get_files]: Getting files containing '%s' for each %4s simulation..." % (flag,m) )  
    # Create empty list to store file names for 
    # each simulation
    files = []
    for i in range(len(dirs)):
        try:
            sfiles = sorted([dirs[i]+proc+'/'+fname for
                          fname in os.listdir(dirs[i]+proc+'/')
                          if flag in fname]) 
            if len(sfiles) < 0.6*num_files:
                print("[get_files]: WARNING, excluding sim that did not complete: %s" %dirs[i])
                dirs[i] = 'removed'  # remove this simulation from the directory list
            else: 
                files.append(sfiles)  

        except OSError: 
            print('[get_files]: WARNING, excluding sim that did not start   : %s ' % (dirs[i]) )
            dirs[i] = 'removed' # remove this simulation from directory list 
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
# get_groups() function
def get_groups(dirs, param, inr3000=False):
    # Given a unique param, returns a list of 
    # groups that pertain to that param 

    print("[get_groups]: Separating files into groups based on param: %s" % (param) ) 

    if inr3000: 
        select = {'rc':['100','500'],'te':['325','345'],
                  'r':['1000','2000','3000'],'fac':['-1.0','0.0','1.0'],
                  'a':['0.0','1.5708']}
    else: 
        select = {'rc':['500'],'te':['325','345'],
                  'r':['1000','2000','3000'],'fac':['-1.0','0.0','1.0'],
                  'a':['0.0','1.5708']}
    
    mygroups = []

    # This code is super ugly, but it works and is fast enough 
    if param == 'fac':
        for te in select['te']: 
            for rc in select['rc']:
                for r in select['r']:
                    for a in select['a']:
                        flag = 'te'+te+'/rc'+rc+'/r'+r+'/a'+a
                        group = [x for x in dirs if flag in x]
                        if group:
                            mygroups.append(group) 
    elif param == 'a':
        for te in select['te']:
            for rc in select['rc']:
                for r in select['r']:
                    for fac in select['fac']:
                        group = []
                        for a in select['a']:
                            flag = 'te'+te+'/rc'+rc+'/r'+r+'/a'+a+'/fac'+fac 
                            for x in dirs:
                               if flag in x:
                                   group.append(x)  
                        if group:
                            mygroups.append(group) 
    elif param == 'r':
        for te in select['te']:
            for rc in select['rc']:
                for a in select['a']:
                    for fac in select['fac']:
                        group = []
                        for r in select['r']:
                            flag = 'te'+te+'/rc'+rc+'/r'+r+'/a'+a+'/fac'+fac 
                            for x in dirs:
                               if flag in x:
                                   group.append(x)  
                        if group:
                            mygroups.append(group) 
    elif param == 'rc':
        for te in select['te']:
            for r in select['r']:
                for a in select['a']:
                    for fac in select['fac']:
                        group = []
                        for rc in select['rc']:
                            flag = 'te'+te+'/rc'+rc+'/r'+r+'/a'+a+'/fac'+fac 
                            for x in dirs:
                               if flag in x:
                                   group.append(x)  
                        if group:
                            mygroups.append(group)
    elif param == 'te':   
        for a in select['a']:
            for rc in select['rc']:
                for r in select['r']:
                    for fac in select['fac']:
                        group = []
                        for te in select['te']:
                            flag = 'te'+te+'/rc'+rc+'/r'+r+'/a'+a+'/fac'+fac 
                            for x in dirs:
                               if flag in x:
                                   group.append(x)  
                        if group:
                            mygroups.append(group)  
    else:
        print('[get_groups]: ERROR, param %s not understood, exiting... \n' % param)
        quit()
    
    # Return the list of params as well
    params = select[param]

    return mygroups, [float(x) for x in params]  

#==================================
# get_groups3d() function
def get_groups3d(dirs, inr3000=False):
    # Given a unique param, returns a list of 
    # groups that pertain to that param 

    if inr3000: 
        select = {'rc':['500'],'te':['325'],
                  'r':['1000','2000','3000'],'fac':['-1.0','0.0','1.0'],
                  'a':['0.0','1.5708']}
    else: 
        select = {'rc':['500'],'te':['325','345'],
                  'r':['1000','2000','3000'],'fac':['-1.0','0.0','1.0'],
                  'a':['0.0','1.5708']}
    
    mygroups = []  # dictionary sorted using   

    afix = '0.0'

    # sort simulations according to r, only for a=0.0 
    for te in select['te']:
        for rc in select['rc']:
            for fac in select['fac']:
                group = []
                for r in select['r']:
                    flag = 'te'+te+'/rc'+rc+'/r'+r+'/a'+afix+'/fac'+fac 
                    for x in dirs:
                       if flag in x:
                           group.append(x)  
                if group:
                    mygroups.append(group) 

    # Return the list of params as well
    params = select['r']

    return mygroups, [float(x) for x in params]  


#==================================
# get_sortdat() function
def get_sortdat(dirs, param, inr3000=False):
    # Given a certain parameter
    # Separates data according to that 
    # parameter

    groups = []
    
    # Decide what param we are plotting
    if inr3000: 
        select = {'rc':[200,400,500,600,800],'te':[265,275,285,325,345],
                  'r':[500,1000,2000,3000],'fac':[-1.0,0.0,1.0,2.0],
                  'a':[0.0,1.5708]}
    else: 
        select = {'rc':[200,400,500,600,800],'te':[265,275,285,325,345],
                  'r':[500,1000,2000],'fac':[-1.0,0.0,1.0,2.0],
                  'a':[0.0,1.5708]}
    try:
        params = select[param]
    except KeyError:
        print('[get_sortdat]: ERROR, param %s not understood, exiting...\n' % param)
        quit()

    for p in params:
        groups.append([x for x in dirs if param+str(p) in x])
    
    return groups, params

#==================================
# get_cmf() function
def get_cmf(dat, x): 
    # Given data, and a x-axis value, x,
    # the cumulative function represents
    # the fraction of data that has value 
    # that is greater than x.
    
    f   = np.zeros(len(x))
    for i in range(0,len(x)):
        f[i]   = len(dat[dat > x[i]])
    
    f   /= len(dat)

    return f 

#==================================
# exp_fit() function
def exp_fit(x, dat, c=0):
    # fits data to decaying exponential
    dat -= c
    dat  = np.log(dat)
    b, log_a = np.polyfit(x,dat,1)
    a = np.exp(log_a)
    return a, b

#==================================
# get_stats() function
def get_stats(files, quant, stat, **kwargs):
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
            ctrl_files = kwargs[key][0]
            nctrl      = len(ctrl_files)
        if key == 'rmn':
            rmn = kwargs[key]
        if key == 'rmx':
            rmx = kwargs[key]
        if key == 'ntrials':
            ntrials = kwargs[key]

    # Total number of simulations
    nsim = len(files)
    tarr = np.zeros(nsim) 
    
    # Create numpy array to store data
    mydata = []
    # Begin main loop over simulations 
    for i in range(nsim):
        # Get length of each simulation
        n    = len(files[i])
        tarr = np.zeros(n)
        
        # time varying data array
        tdats = np.zeros((n,nx1_dom)) # for the simulation
        tdatc = np.zeros((n,nx1_dom)) # for the control
        
        # Error check
        if (nctrl != n):
            print("[get_stats]: WARNING! no. of ctrl files does not equal no. of sim files!")
        #    quit()

        # start time loop
        for j in range(n):
            t, mhvc, rhvc, rpos,\
            acc_rate, facvhvc, ahvc,\
            mcR, mcL,\
            r, ang, vrot,\
            A1, A2         = get_data(files[i][j], **kwargs)
           
            # No need to get control files for timescale measures 
            if 'ts' not in stat:
                t, mhvc, rhvc, rpos,\
                acc_rate, facvhvc, ahvc,\
                mcR, mcL,\
                r, ang, vrot,\
                A1c, A2c       = get_data(ctrl_files[j], **kwargs)

            if quant == 'A1':
                tdats[j] = A1
                if 'ts' not in stat:
                    tdatc[j] = A1c
            elif quant == 'A2':
                tdats[j] = A2
                if 'ts' not in stat:
                    tdatc[j] = A2c 
            else: 
                print('[get_stats]: quant not understood, exiting...')
                quit()
            # Set time array
            tarr[j] = t
                
        # Now do the statistical analysis
        if stat == 'MSE':   
            mydata.append(np.sum((tdats - tdatc)**2.0)/float(tdats.size)) 
        elif stat == 'RMS':
            mydata.append(np.sqrt(np.sum((tdats - tdatc)**2.0)/float(tdats.size)))
        elif stat == 'cmf':
            mydata.append( np.mean(tdats[:, (rmn < r) & (r < rmx)], axis=1) )

        # time scale measures 
        elif stat == 'ts_cut':
            # time scale measured by time above cutoff
            cut = 0.05
            avg = np.mean(tdats[:, (rmn < r) & (r < rmx)], axis=1) 
            mx  = np.max(avg) 
            try:
                ts_cut = tarr[avg > cut][-1] - tarr[avg > cut][0]
            except IndexError:
                ts_cut = 0
            mydata.append( (mx, ts_cut) ) 

        elif stat == 'ts_half':
            # time scale measured by time it takes to go to half the max
            avg  = np.mean(tdats[:, (rmn < r) & (r < rmx)], axis=1) 
            mx   = np.max(avg)
            ts_half = tarr[avg > mx/2.0][-1] - tarr[avg > mx/2.0][0]
            mydata.append( (mx, ts_half) ) 

        elif stat == 'ts_lin':
            # time scale measured by linear fit to data
            ex  = 5 # exclude first 4 points for the linfit 
            avg = np.mean(tdats[:, (rmn < r) & (r < rmx)], axis=1) 
            tarr -= 250.

            m,b = np.polyfit(tarr[ex:],avg[ex:],1)
            if abs(m) > 1e-5:
                mydata.append( ( np.mean(avg),(0.05 - b)/m ) ) 
            else: 
                mydata.append( (np.mean(avg), 0) )

        elif stat == 'ts_spl':
            # time scale measured by spline fit to data
            # if quadratic spline fails, fall back to linear spline
            avg = np.mean(tdats[:, (rmn < r) & (r < rmx)], axis=1)
            avgc= np.mean(tdatc[:, (rmn < r) & (r < rmx)], axis=1) # for the control
            tarr -= 250.
            # start at 50 to avoid immediately falling under the value
            tnew = np.arange(50,10000,1)
            # Try quadratic spline
            f    = uspl(tarr,avg,k=2,s=1.0)
            tcut = tnew[np.where(f(tnew) < np.mean(avgc))]

            # if tcut is empty try linear
            if len(tcut) == 0:
                f = uspl(tarr,avg,k=1,s=1.0)
                tcut = tnew[np.where(f(tnew) < np.mean(avgc))]

            #if tcut[0] == 5:
            #    print(files[i][0])

            # define the timescale
            try: 
                ts = tcut[0] - tarr[1] 
                if ts == 0:
                    print(files[i][0]) 
            except IndexError:
                print('[get_stats]: Unable to get timescale for sim %s'
                      ' setting timescale to 0' % (files[i][0]))
                ts = 0
            mydata.append( (np.max(avg), ts) ) 
        
        elif stat == 'ts_exp':
            # timescale measured by fitting data to a exponential
            # NOTE: this is really only a good fit if we do the 
                #   averaging b/w 1200 - 4000 pc (i.e. disc region)
            avg = np.mean(tdats[:, (rmn < r) & (r < rmx)], axis=1)
            # Assume exponential form after the separation of disc and bar perts 
            mxi  = np.argmax(avg)
            a, b = exp_fit(tarr[mxi:], avg[mxi:])

            # Define the time scale as the half-life
            ts = -(1.0/b) * np.log(2) 

            mydata.append( (np.max(avg), ts) ) 
                
        else: 
            print('[get_stats]: stat not understood, exiting...')
            quit() 

    return tarr, mydata

# main() function 
def main():
    # Now parse input 
    
    # Read in system arguments
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter) 
    parser.add_argument("quant",type=str,
                        help="Quantity options:\n"
                             "   A1: Lopsidedness parameter for m=1\n"
                             "   A2: Lopsidedness parameter for m=2")
    parser.add_argument("--stat", dest="stat",type=str, default='ts_cut',
                        required=False,
                        help="Statistics to carry out on data\n"
                             "  MSE: Mean square error b/w sim and control\n"
                             "  RMS: Root mean square = SQRT(MSE)\n"
                             "  cmf: Plot cumulative function\n")
    parser.add_argument("--rmnmx",type=float,dest='rmnmx',nargs=2,
                        required=False,default=(1200.,4000.),
                        help="Min/max radii for computing <A1> and <A2>\n")
    parser.add_argument("--save", dest="save",action='store_true',
                        default=False,
                        help="Switch to save figure")
    parser.add_argument("--param", dest='param',type=str, default='number',
                        required=False,
                        help='Parameter to plot on x-axis')
    parser.add_argument("--ifrm", dest='ifrm',type=int,
                        default=0,
                        help='Frame of simulation to plot, if doing cmf analysis')
    parser.add_argument("--iani", dest="iani",nargs=2,required=False,
                        default=[0,100],type=int,
                        help="Animate from frame iani[0] to iani[1]\n")
    parser.add_argument("--anim", dest='anim',action='store_true',
                        default=False,
                        help="Switch for animation, if doing cmf analysis")
    parser.add_argument("--inr3000", dest='inr3000',action='store_true',
                        default=False,
                        help="Switch to include r3000 sims, default is False.")
    parser.add_argument("--inm1e8", dest='inm1e8',action='store_true',
                        default=False, 
                        help='Switch to include m1e8 sims, default is False')
    parser.add_argument("--3d", dest="threed", action='store_true',
                        default=False, 
                        help='Switch to make 3d plots showing parameter effects') 

    # parsing arguments
    args = parser.parse_args()
    quant    = args.quant
    stat     = args.stat
    rmn, rmx = args.rmnmx
    save     = args.save
    param    = args.param
    inr3000  = args.inr3000
    ifrm     = args.ifrm
    iani     = args.iani
    anim     = args.anim
    inm1e8   = args.inm1e8 
    threed   = args.threed 
    
    # First find the directories of the simulations
    mdict = get_dirs(os.getcwd(),inr3000,inm1e8)   # dictionary of {mass: mdsim}
    mlist = sorted(os.walk('.').next()[1])  # list of mass directory names
    mlist = [m for m in mlist if 'restart' not in m] # if restart is in this directory 
    if not inm1e8:
        mlist.remove('m1e8') 
    
    # Set the processor and unique flag
    proc = 'id0'
    flag = 'otf'
    
    # Get global parameters from the m0 simulation
    athin =      [fname for fname in os.listdir(mdict[mlist[0]][0]) if 'athinput' in fname]
    n     =  len([fname for fname in os.listdir(mdict[mlist[0]][0]+'/'+proc) if flag in fname])
    prec  = 32

    
    # Read input file
    base, params = read_athinput.readath(mdict[mlist[0]][0] + athin[0])
    nx1_dom = int(params[0].nx1)

    # Get all the paths for all the simulation files     
        # again, this is a list that has the form 
        # mlist 
    af = []
    for m in mlist:
        af.append(get_files(mdict[m], proc, flag, num_files=n)) 


    # Make sure ifrm and iani[1] are not out of range 
    n = len(af[0][0]) 
    if (ifrm > n-1): 
        ifrm = n-1 
    if (iani[1] > n-1):
        iani[1] = n-1 

    # now do the data analysis 
    dats = []
    # exclude the control sim
    for f in af:
        tarr, dat = get_stats(f, quant, stat, ctrl_files = af[0], rmn=rmn, rmx=rmx)
        dats.append(dat) 

    # create dictionaries of results  
    #    first get rid of simulations that didn't finish
    #    the directories of the sims that didn't finish 
    #    are marked by 'removed', c.f. get_files 
    m_dirs_trim = []
    dict_list   = [] 
    for m in mlist:
        m_dirs_trim.append([x for x in mdict[m] if x != 'removed'])
        # key - sim directory: value - mse
    for (md,d) in zip(m_dirs_trim,dats):
        dict_list.append(dict(zip(md,d))) 

    # Stuff for plotting 
    lines  = ['-','--','-.',':']
    colors = ['r','g','b','m','k','c','r']  # colors for each mass 
    markrs = ['s','p','^','o','*','d','<']
    masses = mlist 



    # plotting section 
    if threed: 
        fig = plt.figure(facecolor='white',figsize=(10.,5.))
        #ax  = fig.add_subplot(111, projection='3d') 
        ax  = Axes3D(fig) 
    else:
        fig = plt.figure(facecolor='white',figsize=(10.,5.))
    if stat == 'cmf': 
        #   set parameters for the avals 
        amin  = 0.8*np.min(dats[ 0])
        amax  = 1.2*np.max(dats[1])
        da    = 2*(amax-amin)/len(dats[1])
        avals = np.arange(amin,amax,da)

        
        # sort the data based on param
        groups = []
        for mdt in m_dirs_trim:
            group, params = get_sortdat(mdt, param, inr3000)
            groups.append(group)

        if anim:
            line1 = []
            # Now loop over masses and plot
            for c,mgroup,md,mstr,m in zip(colors,groups,dict_list,masses,markrs):
                for (l,p,g) in zip(lines,params,mgroup):
                    dat   = np.array([md[d][ifrm] for d in g])
                    fvals = get_cmf(dat,avals)
                    line1.append(plt.plot(avals, fvals,c+l+m,
                             label=sym[mstr]+', '+sym[param]+'='+str(p)))
            
            def animate(ifrm):
                # Now loop over masses and plot
                i = 0
                for c,mgroup,md in zip(colors,groups,dict_list):
                    for (l,p,g) in zip(lines,params,mgroup):
                        dat   = np.array([md[d][ifrm] for d in g])
                        fvals = get_cmf(dat,avals)
                        line1[i][0].set_ydata(fvals)
                        i += 1 
                plt.title('t = %1.1f [Myr]' % tarr[ifrm]) 
                return  

            ani = animation.FuncAnimation(fig, animate, range(iani[0],iani[1]), repeat=False)

        else: 
             # Now loop over masses and plot
            for c,mgroup,md,mstr,m in zip(colors,groups,dict_list,masses,markrs):
                for (l,p,g) in zip(lines,params,mgroup):
                    dat   = np.array([md[d][ifrm] for d in g])
                    fvals = get_cmf(dat,avals)
                    plt.plot(avals, fvals,c+l+m,
                             label=sym[mstr]+', '+sym[param]+'='+str(p))

            
            plt.title('t = %1.1f [Myr]' % tarr[ifrm]) 
        
        # These are the same whether doing animation or not 
        plt.xlabel('<'+quant+'>')
        plt.ylabel('Cumulative function, $f$')
        plt.ylim(-0.1,1.1)
        plt.xlim(amin,amax)
        plt.legend()
        

    else: 
        if param == 'number':
            for mdt,md,c,mstr,m in zip(m_dirs_trim,dict_list,colors,masses,markrs):
                # Exclude control simulation
                if mstr != 'm0':
                    xvals = range(len(mdt))
                    for x,d in zip(xvals,mdt):
                        if x == 0:
                            if 'ts' in stat:
                                plt.plot(md[d][0],md[d][1],c+m,label=sym[mstr])
                            else:
                                plt.semilogy(x,md[d],c+m,label=sym[mstr])
                        else:
                            if 'ts' in stat:
                                plt.plot(md[d][0],md[d][1],c+m)
                            else:
                                plt.semilogy(x,md[d],c+m)
            
            plt.ylabel('$\\tau$ [Myr]')
            plt.ylim(-10,600) 
            if 'ts' in stat:
                plt.xlabel('max(<A$_1$>) [unitless]') 
            else:
                plt.xlabel('Simulation number') 
            plt.legend(bbox_to_anchor=(0.,1.02,1.,0.102),loc=3,ncol=4,mode='expand',borderaxespad=0.)
        
        else: 
            # Handle 3D plots 
            if threed: 
                groups = [] 
                for mdt in m_dirs_trim:
                    group, paramr = get_groups3d(mdt, inr3000)
                    groups.append(group)

                for mgroup,mdt,md,c,mstr,l,m in zip(groups,m_dirs_trim,dict_list,colors,masses,lines,markrs):
                    i = 0 # easy way to handle legends  
                    for g in mgroup:
                        if 'ts' in stat:
                            vals = [md[d][1] for d in g]
                        
                        # Get correct value of paramf 
                        if 'fac-1.0' in g[0]:
                            paramf = -1.0
                        elif 'fac0.0' in g[0]:
                            paramf = 0.0
                        else:
                            paramf = 1.0 
                        paramf *= np.ones(len(paramr)) 

                        if i == 0:
                            ax.plot(paramr,paramf,vals,c+m+l,alpha=1.0,label=sym[mstr],linewidth=2.5)
                        else:
                            ax.plot(paramr,paramf,vals,c+m+l,alpha=1.0,linewidth=2.5)
                        i += 1
                # Make plot labels 
                ax.set_xlabel('\n'+sym['r']+' '+units['r'],linespacing=2.2)
                ax.set_ylabel('\n'+sym['fac']+' '+units['fac'],linespacing=2.2)
                ax.zaxis.set_rotate_label(False) 
                ax.set_zlabel('$\\tau$ [Myr]', rotation=90.) 
                # Change default view 
                ax.view_init(elev=30., azim=225.) 
                # set axes limits 
                ax.set_xlim(800.,3200.)
                ax.set_ylim(-1.2,1.2) 
                ax.set_zlim(-10,600.)
                ax.legend() 

            else: 
                # sort the data based on param
                groups = []
                for mdt in m_dirs_trim:
                    group, params = get_groups(mdt, param, inr3000)
                    groups.append(group)
                    plt.ylim(-10,600) 
                
                for mgroup,mdt,md,c,mstr,l,m in zip(groups,m_dirs_trim,dict_list,colors,masses,lines,markrs):
                    i = 0 # easy way to handle legends  
                    for g in mgroup:
                        if 'ts' in stat:
                            vals = [md[d][1] for d in g]
                        if len(vals) == len(params):
                            if i == 0:
                                plt.plot(params,vals,c+m+l,alpha=1.0,label=sym[mstr],linewidth=2.5)
                            else:
                                plt.plot(params,vals,c+m+l,alpha=1.0,linewidth=2.5)
                            i += 1



                plt.xlabel(sym[param]+' '+units[param])
                plt.ylabel('$\\tau$ [Myr]' )
                plt.ylim(-10,600) 
                dx = 0.25*(max(params) - min(params))/( float(len(params)) )  
                plt.xlim(min(params)-dx, max(params)+dx)
                plt.legend(bbox_to_anchor=(0.,1.02,1.,0.102),loc=3,ncol=4,mode='expand',borderaxespad=0.) 
    
    
    if save:
        if anim:
            ani.save(quant+stat+param+'.mp4')
        else:
            plt.savefig(quant+'_'+stat+'_'+param+'.eps')
            plt.savefig(quant+'_'+stat+'_'+param+'.png') 
    else:
        plt.show() 
    
main()

