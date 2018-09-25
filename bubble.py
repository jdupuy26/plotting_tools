import numpy as np
import matplotlib.pyplot as plt 
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
        sizescale = 0.25e3 # [pc]
        self.size = np.random.exponential(size=1,scale=sizescale)[0] 

        # initialize velocity of bubble 
        self.vel  = np.random.uniform(4.,30.) # pc/Myr  

    def apply_bubble(self, quant, data, X, Y):
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

        return bubble_data 
