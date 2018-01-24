#!/usr/bin/python
import numpy as np
import sys
import os
import subprocess as sbp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
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
#  Updated: 01/24/18 
#=====================================================

#=======DICTIONARIES==============

#=================================
# sym dictionary
sym = {"m1e6":"M$_c$=10$^6$ M$_{\odot}$", "m1e5": "M$_c$=10$^5$ M$_{\odot}$",
       "a":"$\phi_{c,0}$", "fac":"f$_{c}$", "te":"t$_e$",
       "r":"r$_{pos,0}$", "rc":"r$_c$", "f":"$f$", "A1":"<A1>",
       "A2":"<A2>", "RoL": "M$_{cr,R}$/M$_{cr,L}$", "LoR":"M$_{cr,L}$/M$_{cr,R}$"} 

# units dictionary
units = {"a":"[rad]","rc":"[pc]","fac":"[unitless]","r":"[pc]",
         "te":"[Myr]"} 
 

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
# get_groups() function
def get_groups(dirs, param, inr3000=False):
    # Given a unique param, returns a list of 
    # groups that pertain to that param 

    print("[get_groups]: Separating files into groups based on param: %s" % (param) ) 

    if inr3000: 
        select = {'rc':['200','400','600'],'te':['265','275','285'],
                  'r':['500','1000','2000','3000'],'fac':['0.0','1.0','2.0'],
                  'a':['0.0','1.5708']}
    else: 
        select = {'rc':['200','400','600'],'te':['265','275','285'],
                  'r':['500','1000','2000'],'fac':['0.0','1.0','2.0'],
                  'a':['0.0','1.5708']}
    
    mygroups = []

    # This code is super ugly, but it works and is fast enough 
    if param == 'fac':
        for te in select['te']: 
            for rc in select['rc']:
                for r in select['r']:
                    for a in select['a']:
                        flag = 'te'+te+'/rc'+rc+'/r'+r+'/a'+a
                        mygroups.append([x for x in dirs if flag in x]) 
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
                        mygroups.append(group)  
    else:
        print('[get_groups]: ERROR, param %s not understood, exiting... \n' % param)
        quit()
    
    # Return the list of params as well
    params = select[param]

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
            print(ctrl_files) 
            print("[get_stat]: Error! no. of ctrl files does not equal no. of sim files!")
            quit()

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
                tdats[j] = A1
                tdatc[j] = A1c
            elif quant == 'A2':
                tdats[j] = A2
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
    parser.add_argument("--stat", dest="stat",type=str, default='MSE',
                        required=False,
                        help="Statistics to carry out on data\n"
                             "  MSE: Mean square error b/w sim and control\n"
                             "  RMS: Root mean square = SQRT(MSE)\n"
                             "  cmf: Plot cumulative function\n")
    parser.add_argument("--rmnmx",type=float,dest='rmnmx',nargs=2,
                        required=False,default=(10.,5000.),
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
    m0f   = get_files(m0d  ,proc,flag) # ctrl_files 
    m1e5f = get_files(m1e5d,proc,flag)
    m1e6f = get_files(m1e6d,proc,flag)

    # Make sure ifrm and iani[1] are not out of range 
    n = len(m0f[0]) 
    if (ifrm > n-1): 
        ifrm = n-1 
    if (iani[1] > n-1):
        iani[1] = n-1 

    # now do the data analysis 
    tarr, dat1e5 = get_stats(m1e5f, quant, stat, ctrl_files = m0f,rmn=rmn,rmx=rmx)
    tarr, dat1e6 = get_stats(m1e6f, quant, stat, ctrl_files = m0f,rmn=rmn,rmx=rmx) 
   
    # create dictionaries of results  
    #    first get rid of simulations that didn't finish
    #    the directories of the sims that didn't finish 
    #    are marked by 'removed', c.f. get_files 
    m1e5d_trim = [x for x in m1e5d if x != 'removed']
    m1e6d_trim = [x for x in m1e6d if x != 'removed']
        # key - sim directory: value - mse 
    dict1e5 = dict(zip(m1e5d_trim,dat1e5)) 
    dict1e6 = dict(zip(m1e6d_trim,dat1e6))

    
    # plotting section 
    fig = plt.figure(facecolor='white')
    if stat == 'cmf': 
        #   set parameters for the avals 
        amin  = 0.8*np.min(dat1e5)
        amax  = 1.2*np.max(dat1e6)
        da    = 2*(amax-amin)/len(dat1e6)
        avals = np.arange(amin,amax,da)

        lines = ['-','--','-.',':']

        # sort the data based on param
        groups1e5, params = get_sortdat(m1e5d_trim, param, inr3000) 
        groups1e6, params = get_sortdat(m1e6d_trim, param, inr3000)

        if anim:
            line1 = []
            for (l,p,g) in zip(lines,params,groups1e5):
                dat   = np.array([dict1e5[d][ifrm] for d in g])
                fvals = get_cmf(dat,avals)
                line1.append(plt.plot(avals, fvals,'r'+l+'x',
                         label=sym['m1e5']+', '+sym[param]+'='+str(p)))
            line2 = [] 
            for (l,p,g) in zip(lines,params,groups1e6):
                dat   = np.array([dict1e6[d][ifrm] for d in g])
                fvals = get_cmf(dat,avals)
                line2.append(plt.plot(avals, fvals,'g'+l+'+',
                             label=sym['m1e6']+', '+sym[param]+'='+str(p)))
            def animate(ifrm):
                iline = range(len(line1))
                for (il,l,p,g) in zip(iline,lines,params,groups1e5):
                    dat   = np.array([dict1e5[d][ifrm] for d in g])
                    fvals = get_cmf(dat,avals)
                    line1[il][0].set_ydata(fvals)
                for (il,l,p,g) in zip(iline,lines,params,groups1e6):
                    dat   = np.array([dict1e6[d][ifrm] for d in g])
                    fvals = get_cmf(dat,avals)
                    line2[il][0].set_ydata(fvals)
                plt.title('t = %1.1f [Myr]' % tarr[ifrm]) 
                return  

            ani = animation.FuncAnimation(fig, animate, range(iani[0],iani[1]), repeat=False)

        else: 
            for (l,p,g) in zip(lines,params,groups1e5):
                dat   = np.array([dict1e5[d][ifrm] for d in g])
                fvals = get_cmf(dat,avals)
                plt.plot(avals, fvals,'r'+l+'x',
                         label=sym['m1e5']+', '+sym[param]+'='+str(p))
            
            for (l,p,g) in zip(lines,params,groups1e6):
                dat   = np.array([dict1e6[d][ifrm] for d in g])
                fvals = get_cmf(dat,avals)
                plt.plot(avals, fvals,'g'+l+'+',
                         label=sym['m1e6']+', '+sym[param]+'='+str(p))

            plt.title('t = %1.1f [Myr]' % tarr[ifrm]) 
        
        # These are the same whether doing animation or not 
        plt.xlabel('<'+quant+'>')
        plt.ylabel('Cumulative function, $f$')
        plt.ylim(-0.1,1.1)
        plt.xlim(amin,amax)
        plt.legend()
        

    else: 
        if param == 'number':
            x = range(len(m1e5d_trim))
            for x,d in zip(x,m1e5d_trim):
                plt.semilogy(x,dict1e5[d],'rx')
            x = range(len(m1e6d_trim))
            for x,d in zip(x,m1e6d_trim):
                plt.semilogy(x,dict1e6[d],'g+')
            plt.ylabel(stat+'('+quant+', '+quant+'$_{control}$)')
            plt.xlabel('Simulation number') 
        
        else: 
            groups1e5, params = get_groups(m1e5d_trim, param, inr3000)
            groups1e6, params = get_groups(m1e6d_trim, param, inr3000)

            for g in groups1e5:
                vals = [dict1e5[d] for d in g]
                if len(vals) == len(params):
                    plt.semilogy(params, vals,'r-') 
            for g in groups1e6:
                vals = [dict1e6[d] for d in g]
                if len(vals) == len(params):
                    plt.semilogy(params, vals,'g-') 
            
            plt.xlabel(sym[param]+' '+units[param])
            plt.ylabel(stat+'('+quant+', '+quant+'$_{control}$)')
    
    
    if save:
        plt.savefig(quant+stat+param+'.eps') 
    else:
        plt.show() 
    
main()

