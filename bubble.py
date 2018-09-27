import numpy as np
import matplotlib.pyplot as plt 
from scipy import ndimage 

# the following file assumes units: [pc, Myr, & solar masses] throughout

class Bubble(object):
    """ 
    Bubble class for a single HI bubble/shell 

    Attributes:
        x1min: minimum grid value in x1-dir
        x1max: maximum grid value in x1-dir
        x2min: ''
        x2max: '' 

        pos: position of bubble center in polar coords
        size: size of the bubble (radius)   
        d: density of bubble  (np array of size (nx2, nx1))
        v: velocity of bubble (np array of size (nx2, nx1))
    """

    def __init__(self, params=(1.e3,1.5e4,0.0,2.*np.pi)):
        """
        Initialize attributes of grid 
        """
        ( self.x1min, self.x1max, 
          self.x2min, self.x2max )  = params


        # initialize bubble position & size
        self.pos = np.array([ np.random.uniform(self.x1min,self.x1max), 
                              np.random.uniform(self.x2min,self.x2max) ])
        sizescale = 0.175e3 # [pc]
        #sizescale = 1.0e3
        self.size = np.random.exponential(size=1,scale=sizescale)[0] 

        # initialize velocity of bubble 
        self.vel  = np.random.uniform(10.,30.) # [pc/Myr]  

        # initialize rotation rate 
        self.omg  = 0.0337 # [1/Myr] 

    def apply_bubble(self, quant, data, X, Y, mode='default', bub_dens=None, fTeq=None):
        """
        Takes in a single snapshot of 
        simulation data and places the bubble on the
        data. 

        Inputs:
            self: bubble object
            quant: str corresponding to type of data (e.g. 'd', 'v', etc.)
            data: array of grid data
            X: array of cell-centered X values (size (nx2,nx1))
            Y: array of cell-centered Y values (size (nx2,nx1))
            mode: mode for applying bubble
            bub_dens: array of bub_dens data (to be used to get Teq)
            fTeq: 1D interpolating function for equilibrium temp calculation.
                  this is computed in plot_sims from the file mean_Teq.dat.
                  It is an average over the metallicity gradient so that 
                  Teq = Teq(d) only.   

        Output:
            bubble_data: data after applying bubble 
        """ 
        # Preliminaries
        xpos, ypos = (self.pos[0]*np.cos(self.pos[1]),  # x,y position of bubble
                      self.pos[0]*np.sin(self.pos[1]))

        # get mask to apply bubbles 
        mask = ((X - xpos)**2. + (Y - ypos)**2.) < self.size**2.
        bubble_data = data.copy()

        # quantity specific things
        if mode == 'default':
            if quant == 'd':
                bubble_data[:,mask] -= data[:,mask] - 1.0 # add small background density  

            elif quant == 'v1':
                bubble_data[:,mask] -= data[:,mask] - 1.0 # add small velocity background
                
            elif quant == 'v2':
                bubble_data[:,mask] -= data[:,mask] - 1.0 # add small velocity background  

            elif quant == 'T':
                bubble_data[:,mask] -= data[:,mask] - 1e7 # add high temperature 

            else:
                raise Exception('[UnknownQuant]: do not know how to apply {} to data'.format(quant))

        # high density shell w/ expansion velocity  
        elif mode == 'shell':
            # step 1: convert mask to binary  
            bindata = mask.astype(int) 
            # step 2: compute boundary mask
            bmsk   = ndimage.laplace(bindata) > 0 
            # step 3: get size of boundary mask
            bmsk_size = float(bubble_data[:,bmsk].size)

            if quant == 'd':
                bubble_data[:,mask] -= data[:,mask] - 1.0
                # place a high density shell at the boundary  
                if bmsk_size > 0: 
                    bubble_data[:,bmsk]  = np.sum(data[:,mask])/bmsk_size  

            elif quant == 'v1' or quant == 'v2':
                bubble_data[:,mask] -= data[:,mask] - 1.0

                # now compute the expansion velocity for the boundary
                    # step 1: get coords of boundary in sim frame
                xbnds = X[bmsk]
                ybnds = Y[bmsk] 
                rbnds = np.sqrt( xbnds**2. + ybnds**2. )
                    # step 2: get coords of boundary in bub frame
                xpbnds = xbnds - xpos
                ypbnds = ybnds - ypos    
                    # step 3: get gamma angle 
                gam = np.arctan2(ypbnds, xpbnds)
                    # step 4: compute vx, vy 
                vx = self.vel*np.cos(gam)
                vy = self.vel*np.sin(gam) 
                    # step 5: convert to R, \phi vel
                vr = (xbnds/rbnds)*vx + (ybnds/rbnds)*vy
                vp = (xbnds/rbnds)*vy - (ybnds/rbnds)*vx
                    # step 6: set velocity @ the boundary 
                if quant == 'v1':
                    bubble_data[:,bmsk] = vr
                else:
                    bubble_data[:,bmsk] = vp - self.omg*rbnds # convert to vel in rotating frame  

            elif quant == 'T':
                bubble_data[:,mask] -= data[:,mask] - 1e7

                # now compute equilibrium temperature for the boundary (requires bub_dens)  
                if bub_dens is None: 
                    bub_dens = self.apply_bubble('d', data, X, Y, mode=mode)  

                if bmsk_size > 0: 
                    bubble_data[:,bmsk]  = fTeq(bub_dens[:,bmsk])  
            
            else:
                raise Exception('[UnknownQuant]: do not know how to apply {} to data'.format(quant))

        else:
            raise Exception('[UnknownMode]: mode {} not understood'.format(mode)) 
        return bubble_data 
