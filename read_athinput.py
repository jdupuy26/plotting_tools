#!/usr/bin/python
import numpy as np
import sys


#==========================================
"""
 Code: read_athinput.py

 Purpose: Read in parameters file from athinput.----
          and store vital simulation mesh information in 
          a more compact structure. Also pull the base
          file name for the files to be read in


 Keywords: quant [string] - quantity to be

 Warning: 

      
 paramsrams[i] = { 
                   level        : resolution of level
                   nx1          : number of x1 zones
                   nx2          : number of x2 zones
                   nx3          : number of x2 zones
                   x1min        : minimum x1 
                   x1max        : maximum x1
                   ilog         : switch for log grid 
                   x1rat        : ratio of each cell to the next cell in x1
                   x2min        : minimum x2
                   x2max        : maximum x2
                   x3min        : minimum x3
                   x3max        : maximum x3
                   ngridx1      : number of proc. divisions in x1
                   ngridx2      : number of proc. divisions in x2
                   ngridx3      : number of proc. divisions in x3
                   idisp        : i displacement of level  (only SMR)
                   jdisp        : j displacement of level  (only SMR)
                   kdisp        : k displacement of level  (only SMR)
                   omg          : bar rotation [1/Myr] 
                   gam          : Gamma value (C_P/C_V)
                   mu           : mean molecular weight
                   scaleh       : scale height [pc] 
                   Nang         : no. of angles for analysis 
                   thvcs        : start time of hvc accretion [Myr]
                   thvce        : end time of hvc accretion [Myr] 
                 }
        JLD edited: (9/3/17) 
                    get x1rat and omg from athinput 
                    (10/11/17)
                    get mu, gam, and scaleh from athinput 
                    (11/8/17)
                    get Nang for no. of angles of analysis
                    (11/14/17)
                    get thvcs, thvce

"""
#==========================================
              
def readath(fl):
    #create structure class
    class athparams:
        def __init__(self):
            self.level=0
            self.nx1=0
            self.nx2=0
            self.nx3=0
            self.x1min=0
            self.x1max=0
            self.ilog=0
            self.x1rat=0
            self.x2min=0
            self.x2max=0
            self.x3min=0
            self.x3max=0
            self.ngridx1=0
            self.ngridx2=0
            self.ngridx3=0
            self.idisp=0
            self.jdisp=0
            self.kdisp=0
            self.omg=0
            self.gam=0
            self.mu=0
            self.scaleh=0
            self.Nang=0
            self.thvcs=0
            self.thvce=0
        def __repr__(self):
            return "<athparams level:%s nx1:%s nx2=%s nx3=%s x1min=%s x1max=%s ilog=%s x1rat=%s x2min=%s x2max=%s x3min=%s x3max=%s ngridx1=%s ngridx2=%s ngridx3=%s idisp=%s jdisp=%s kdisp=%s gamma=%s iso_csound=%s omg=%s \
                     gam=%s mu=%s scaleh=%s Nang=%s thvcs=%s thvce=%s>" % (self.level, self.nx1,self.nx2,self.nx3,self.x1min,self.x1max,self.ilog,self.x1rat,self.x2min,self.x2max,
                                                                           self.x3min,self.x3max,self.ngridx1,self.ngridx2,self.ngridx3,self.idisp,self.jdisp,self.kdisp,
                                                                           self.omg,self.gam,self.mu,self.scaleh,self.Nang,self.thvcs,self.thvce)


    #Read in athinput file into array of lines
    try:
       file = open(fl)
       content = file.readlines()
    except:
        print 'read_athinput: failed to read file input'
        quit()


    #strip the whitespace
    for i in range(0,len(content)-1):
        content[i] = content[i].strip()

    jobind = content.index('<job>');


    #find number of domains and problem ID
    num=0
    flag=0
    jobind+=1
    while (flag == 0 and jobind < len(content)):
        splitstring = content[jobind].split()
        if splitstring[0]  == 'problem_id':
            base = splitstring[2]
    
        elif splitstring[0] == 'num_domains':
            num = int(splitstring[2])
            flag=1
        jobind +=1


    if num == 0: 
        print 'read_athinput: NO DOMAIN KEYWORD FOUND'
        quit()
    #Create sturcture for holding data for all domains
    params = [athparams()]
    for x in range(num-1):
        params.append(athparams())
    #populate structure with parameter file values
    check=0
    for i in range(1,num+1):
        ind = content.index('<domain'+str(i)+'>')
        while (check ==0):
            ind +=1 
            if ind >  len(content)-1:
                break
            
            line = content[ind]
            
            if line == '<domain'+str(i+1)+'>':
                check=1
            #if line == '<problem>':
            #    check=1

            splitline = content[ind].split()
        
            if len(splitline) > 0:
                if splitline[0] == 'level':
                    params[i-1].level   = float(splitline[2])
                if splitline[0] == 'Nx1':
                    params[i-1].nx1     = float(splitline[2])
                if splitline[0] == 'Nx2':
                    params[i-1].nx2     = float(splitline[2])
                if splitline[0] == 'Nx3':
                    params[i-1].nx3     = float(splitline[2])
                if splitline[0] == 'x1min':
                    params[i-1].x1min   = float(splitline[2])
                if splitline[0] == 'x1max':
                    params[i-1].x1max   = float(splitline[2])
                if splitline[0] == 'ilog':
                    params[i-1].ilog    = float(splitline[2])
                if splitline[0] == 'x1rat':
                    params[i-1].x1rat   = float(splitline[2])
                if splitline[0] == 'x2min':
                    params[i-1].x2min   = float(splitline[2])
                if splitline[0] == 'x2max':
                    params[i-1].x2max   = float(splitline[2])
                if splitline[0] == 'x3min':
                    params[i-1].x3min   = float(splitline[2])
                if splitline[0] == 'x3max':
                    params[i-1].x3max   = float(splitline[2])
                if splitline[0] == 'NGrid_x1':
                    params[i-1].ngridx1 = float(splitline[2])
                if splitline[0] == 'NGrid_x2':
                    params[i-1].ngridx2 = float(splitline[2])
                if splitline[0] == 'NGrid_x3':
                    params[i-1].ngridx3 = float(splitline[2])
                if splitline[0] == 'iDisp':
                    params[i-1].idisp   = float(splitline[2])
                if splitline[0] == 'jDisp':
                    params[i-1].jdisp   = float(splitline[2])
                if splitline[0] == 'kDisp':
                    params[i-1].kdisp   = float(splitline[2])
                if splitline[0] == 'omg':
                    params[i-1].omg     = float(splitline[2])
                if splitline[0] == 'gamma':
                    params[i-1].gam     = float(splitline[2])
                if splitline[0] == 'mu':
                    params[i-1].mu      = float(splitline[2])
                if splitline[0] == 'scaleh':
                    params[i-1].scaleh  = float(splitline[2])
                if splitline[0] == 'Nang':
                    params[i-1].Nang    = float(splitline[2])
                if splitline[0] == 'thvcs':
                    params[i-1].thvcs   = float(splitline[2])
                if splitline[0] == 'thvce':
                    params[i-1].thvce   = float(splitline[2])

        check=0 

    return base,params
