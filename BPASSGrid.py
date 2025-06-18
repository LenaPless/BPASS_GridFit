
from astropy.table import Table, vstack

import pimodels as pi
from pimodels import PIModel as PI
from pifit import MultiNestFit as mf

import os
import numpy as np


class Grid(PI.PIModel):
    
    """
    Class to create pimodeling objects for a configured BPASS grid
    This class inherits from PI.PIModel and sets up the BPASS model
    """

    def __init__(self, name='bpass', pars={['CO', 'LOGZ', 'LOGU', 'XI', 'NH']}, lines=['OIII4959', 'OIII5007', 'HB', 'HA', 'OII3727'], fix=None, path=None, recalculate=False):
        
        """
        Initialize the BPASS model object
        
        name - name of the model
        pars - list of parameters to probe in the grid, default is ['CO', 'LOGZ', 'LOGU', 'XI', 'NH']
        lines - list of emission lines to include in the model, default is ['OIII4959', 'OIII5007', 'HB', 'HA', 'OII3727']
        fix - dictionary with parameters to fix in the grid
        path - path to the BPASS grid file, if None, it will use the default path
        recalculate - boolean, if True, it will recalculate the abundances from Gutkin et al. (2016)
        """
        
        PI.PIModel.__init__(self, name)

        self.recalculate = recalculate
        self.abundance_setup(recalculate=self.recalculate)
        self.path = path
        
        self.pars = pars
        self.fix = fix
        
        self.read_grid()
    
    def dirs(self):
        
        """
        Get the directory paths for the BPASS model
        """
        
        if self.path is not None:
            ROOT = os.path.dirname(self.path)
        else:
            ROOT = os.getcwd()

        return {"ROOT": ROOT}
    
    def abundance_setup(self, recalculate=False):
        
        """
        Set up abundances using Gutkin et al. (2016)
        """

        if self.recalculate:
            dirs = self.dirs()
            try:
                t = Table().read(dirs['ROOT']+'gutkin_abun.dat', format='ascii', header_start=0)
            except FileNotFoundError:
                raise FileNotFoundError("Gutkin abundance file not found. Please locate the file named 'gutkin_abun.dat' and move it to the current directory.")

            names = t['element']

            ninH = t['abundance']
            X, Y, Z = (0.70919985,  0.27556015,  0.01524000)
                
            absun = pi.AbundanceSet.AbundanceSet(names, ninH, X, Y, Z)

            X, Y, Z = absun.calculate_XYZ()
            self.Zsun = absun.Z
            self.Xsun = absun.X
            self.Ysun = absun.Y
        
        else:
            self.Zsun = 0.01524
            self.Xsun = 0.70919985
            self.Ysun = 0.27556015
    
    def read_grid(self):
        
        """
        Read the BPASS grid and set up the model
        """
        
        # read in the grid
        
        if self.path is not None:
            try:
                self.grid = Table.read(self.path)
            except FileNotFoundError:
                raise FileNotFoundError("BPASS grid file not found. Please provide a valid path.")
        else:
            try:
                self.grid = Table.read(self.dirs()['ROOT'] + 'BPASS_grid.fits')
            except FileNotFoundError:
                raise FileNotFoundError("BPASS grid file not found. Please provide a valid path.")
        
        for p in self.pars:
            self.parameters[p] = self.grid[p]
        self.parameter_order = self.pars
        
        for l in self.lines:
            self.add_predicted(l, self.grid[l], scale=True)
        
        self.parameters['LOGZ_ZSUN'] = self.parameters['LOGZ']-np.log10(self.Zsun)   
        
        if self.fix is not None:
            self.subset_grid(self.fix)