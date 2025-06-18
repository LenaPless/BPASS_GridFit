import os
import warnings as w

from pifit import MultiNestFit as mf

class Fit():
    
    """
    Class to run a fit using the MultiNestFit class from pifit.
    This class requires a BPASS_GridFit model grid and target objects with fluxes and errors.
    """
    
    def __init__(self, target, model, lines=None, path=None, version=1):

        """
        Initialize the fit object.
        
        target - Target object with fluxes and errors
        model - BPASS_GridFit model object
        lines - list of emission lines to fit, if None, it will use the lines present in the target object
        path - path to save the results, if None, it will use the default path './results/'
        version - version of the fit, default is 1
        """
        
        self.target = target
        self.model = model
        self.path = path
        self.id = target.id
        self.v = version

        if self.path is not None:
            try:
                os.path.exists(self.path)
            except:
                w.warn(f"Path {self.path} does not exist. Directory will be created.")
        else:
            try:
                os.path.exists('./results/')
            except:
                w.warn("Default path './results/' does not exist. Directory will be created.")
        
        self.target.read_flux()
        self.target.lines_to_fit()
        
        #only fit lines that are marked as True in the target's fit_lines and thus have a SNR above the threshold (DEFAULT: 3)
        
        self.lines = [x for x, y in zip(lines, self.target.fit_lines) if y == 'True']
        
        self.run_fit()
    
    def run_fit(self):
        
        """
        Run the fit.
        """
        
        self.fit = mf.MultiNestFit(self.model)
        
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        
        self.results = self.fit.fit(lines=self.lines, fluxes=self.target.flux, dfluxes=self.target.error, fit_dust=False, basename=self.path+f'{self.id:.0f}/v{self.v}_{self.id:.0f}')
        self.summary = self.results.summarised_results()
        
        return self.summary
        
    def show_corner(self):
        
        """
        Show the corner plot of the fit.
        """
        
        self.results.show_triangle()