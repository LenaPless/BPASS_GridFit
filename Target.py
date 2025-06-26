import pandas as pd
import numpy as np

import difflib

import warnings as w

class Target():
    
    """
    Class to create a target object for a BPASS grid fit.
    This class holds the target ID, data frame with fluxes, and the lines to fit.
    It also provides methods to read the fluxes and errors for the target, and to determine which lines to fit based on the signal-to-noise ratio.
    """

    def __init__(self, id, df=None, lines=None, snr=3):
        
        """
        Initialize the Target object.
        id - target ID
        df - pandas DataFrame with fluxes and errors for the target
        lines - list of emission lines to fit, if None, it will use the lines from the DataFrame
        snr - signal-to-noise ratio threshold for fitting lines, default is 3
        """
        
        self.id = id
        self.df = df
        self.lines = lines
        self.snr = snr
        self.fit_lines = [False]*len(self.lines) if lines is not None else None
        
        self.read_flux()

    def lines_to_fit(self):
        
        """
        Return the lines to fit for the target based on the signal-to-noise ratio.
        """

        for i, l in enumerate(self.lines):
            if self.flux[i]/self.error[i] > self.snr:
                self.fit_lines[i] = True
    
    def filter_flux(self):
        
        """
        Filter the fluxes and errors for the target based on the lines to fit.
        This will return only the fluxes and errors for the lines that are marked as True in fit_lines.
        """
        
        if self.fit_lines is None:
            w.warn("No lines to fit. Please check the input DataFrame and lines.")
        else:
            self.flux = [self.flux[i] for i, fit in enumerate(self.fit_lines) if fit]
            self.error = [self.error[i] for i, fit in enumerate(self.fit_lines) if fit]


    def read_flux(self):
        
        f = []
        err = []
        
        for l in self.lines:
            """
            Read the fluxes and errors for the target from the DataFrame.
            l - emission line to read
            """

            if l not in self.df.columns:
                try:
                    for k in self.df.keys():
                        if difflib.SequenceMatcher(None, k, l).ratio() > 0.8:
                            l = k
                            w.warn(f"Found similiar line label {k}. Using {k} instead of {l}. If you did not mean to do this, please check your input.")
                except ValueError:
                    raise ValueError(f"Line {l} not found in the DataFrame columns.")

            if '372' in l and ('O2' in l or 'o2' in l or 'OII' in l or 'oII' in l):
                
                w.warn(f"Line {l} is a blend of OII lines commonly unresolved. The fluxes will be summed for the OII doublet at 3726 and 3729 A. If you did not mean to do this, please check your input.")

                f.append(self.df['o2_3726'].loc[self.df['ID']==self.id].values[0]+self.df['o2_3729'].loc[self.df['ID']==self.id].values[0])
                e = np.sqrt(self.df['o2_3726_err_new'].loc[self.df['ID']==self.id].values[0]**2+self.df['o2_3729_err_new'].loc[self.df['ID']==self.id].values[0]**2)
                err.append(e)

            elif '671' in l and ('S2' in l or 's2' in l or 'SII' in l or 'sII' in l):
                
                w.warn(f"Line {l} is a blend of SII lines commonly unresolved. The fluxes will be summed for the SII doublet at 6718 and 6733 A. If you did not mean to do this, please check your input.")

                f.append(self.df['s2_6718'].loc[self.df['ID']==self.id].values[0]+self.df['s2_6733'].loc[self.df['ID']==self.id].values[0])
                e = np.sqrt(self.df['s2_6718_err_new'].loc[self.df['ID']==self.id].values[0]**2+self.df['s2_6733_err_new'].loc[self.df['ID']==self.id].values[0]**2)
                err.append(e)
                
            else:
                f.append(self.df[l].loc[self.df['ID']==self.id].values[0])
                err.append(self.df[f'{l}_err_new'].loc[self.df['ID']==self.id].values[0])

        self.flux = f
        self.error = err

        self.lines_to_fit()
        self.filter_flux()
        
        return self.flux, self.error
