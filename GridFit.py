import os
import warnings as w

import difflib
import csv

import numpy as np

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
            self.path = './results/'
            try:
                os.path.exists('./results/')
            except:
                w.warn("Default path './results/' does not exist. Directory will be created.")
        
        #only fit lines that are marked as True in the target's fit_lines and thus have a SNR above the threshold (DEFAULT: 3)
        
        self.lines = [x for x, y in zip(lines, self.target.fit_lines) if y == True]
        
        self.results = None
        self.summary = None
    
    def translate_line_labels(self):
        """
        Translate the line labels to the format used by the model.
        This is necessary because the model uses different labels for some lines.
        """
        
        for i, line in enumerate(self.lines):
            if line not in self.model.lines:
                try: 
                    for k in self.model.lines:
                        if difflib.SequenceMatcher(None, k.lower(), line.lower()).ratio() >= 0.4:
                            self.lines[i] = k
                            w.warn(f"Found similiar line label {k}. Using {k} instead of {line}. If you did not mean to do this, please check your input.")
                except ValueError:
                    raise ValueError(f"Line {line} not found in the Model Grid.")
    
    def run_fit(self):
        
        """
        Run the fit.
        """

        if self.lines != []:
            self.fit = mf.MultiNestFit(self.model)
        
            if not os.path.exists(self.path):
                os.makedirs(self.path)

            self.results = self.fit.fit(lines=self.lines, fluxes=np.array(self.target.flux), dfluxes=np.array(self.target.error), fit_dust=False, basename=self.path+f'{self.id:.0f}/v{self.v}_{self.id:.0f}')
            self.summary = self.results.summarised_results()
            
            return self.summary

        else:
            w.warn("No lines to fit. Consider lowering the SNR threshold in the target object or providing a list of lines to fit.")
            return None
        
    def show_corner(self):
        
        """
        Show the corner plot of the fit.
        """
        if self.results is None:
            w.warn("No results to show.")
            return
        else:
            self.results.show_triangle()
    
    def save_results(self):
        
        """
        Save the results of the fit to a file.
        The results will be saved in the path specified during initialization.
        """
        
        if self.results is None:
            w.warn("No results to save. Please run the fit first.")
            return
        
        self.summary['target_id'] = self.id
        self.summary['version'] = self.v
        self.summary['lines'] = ', '.join(self.lines)
        self.summary['snr'] = self.target.snr
        self.summary['flux'] = ', '.join([str(x) for x in self.target.flux])
        self.summary['error'] = ', '.join([str(x) for x in self.target.error])
        
        with open(self.path+f'{self.id:.0f}/v{self.v}_{self.id:.0f}_summary.csv', 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=self.summary.keys())
            writer.writeheader()
            writer.writerow(self.summary)