#!/usr/bin/env python
import numpy as np
import sys
import os
import subprocess as sbp
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

import plot_sims as ps
import plot_lv_sii as plv

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
                        default=False, required=False,
                        help="Switch to return only stitched together array\n"
                             "To be used if this file is imported from another file\n")
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


    return parser.parse_args() 

def main():
    cwd  = os.getcwd() 
    dirs = get_dirs()
    args = get_args() 


    for d in dirs:
        os.chdir(d)
        try:  
            x = plv.main(args) 
        except OSError:
            print('Sim files not present for %s!' %(d)) 


    return 

main() 
