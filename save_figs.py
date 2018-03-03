#!/usr/bin/env python 
import numpy as np
import sys
import os
import subprocess as sbp 
import shutil
import argparse
from argparse import RawTextHelpFormatter


#=====================================================
#
#  Code: save_figs.py
#
#  Purpose:  Loops through simulation directories
#            create and images   
#
#  Keywords: python save_figs.py -h   
#
#  Usage: python save_figs.py   
#
#  WARNING: THIS MUST BE RUN FROM SIMULATION DIRECTORY 
#          
#  Author: John Dupuy 
#          UNC Chapel Hill
#  Date:    03/01/18 
#  Updated: 03/01/18 
#=====================================================


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

    # find mass directories 
    masses = os.walk('.').next()[1]
    # create empty dictionary of emtpy lists
        # for each mass  
    mdict  = {m: [] for m in masses}

    dirlist = []
    
    print('\n[get_dirs]: Finding simulation directories...')
    # Loop through stdout and get sims
    for d in stdout.splitlines():
        if ('fac' in d) and ('id' not in d):
            for m in masses:
                if m in d:
                    dirlist.append(d+'/')

    return dirlist 


def main():
    
    # First get simulation directories 
    dirs = get_dirs(os.getcwd()) 
    # Define the plotting path
    plotpath = '/afs/cas.unc.edu/users/j/d/jdupuy26/'+\
               'Johns_work/misc_scripts/plotting_tools/'

    # Read in system arguments 
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument('ptype',type=str,
                        help="Type of plot to save, options:\n"
                             "      animd: animation of density\n"
                             "     paneld: panel plot of density\n"
                             "         A1: average quantity of A1")
    parser.add_argument('--fmt',type=str,default='eps',
                        help="Type of figure to save, default: eps") 
    
    # parsing arguments
    args  = parser.parse_args()
    ptype = args.ptype
    fmt   = args.fmt

    # Get plot file  
    if ptype == 'animd' or ptype == 'paneld':
        pfile = 'plot_sims.py'
        if ptype == 'animd':
            pcall = 'python plot_sims.py d --anim --save --ctrace --qminmax '
        else: 
            pcall = 'python plot_sims.py d --ifrm 2 5 10 20 30 40 --save --fmt '+ fmt +' --ctrace --qminmax ' 
    elif ptype == 'A1':
        pfile = 'plot_otf.py'
        pcall = "python plot_otf.py '<A1>' --save --fmt " + fmt  
    else: 
        print('[main]: ptype %s not understood, exiting...\n' %(ptype) ) 
        quit() 
    # Define qmnmx 
    qmnmx = ''

    # Now loop through each directory and create figure 
    for d in dirs:
        os.chdir(d) # change to sim directory
        
        # Copy the python script to sim directory  
        shutil.copy2(plotpath+pfile, d) 

        # Depending on ptype and mass choose qminmax  

        if ptype != 'A1':
            if 'm5e7' in d or 'm1e8' in d:
                qmnmx = '0 40'
            else: 
                qmnmx = '0 20'

        # Create the figure
        sbp.call(pcall+qmnmx,shell=True) 

        # Delete the python script
        sbp.call('rm '+pfile,shell=True)
         
    return 
main() 
