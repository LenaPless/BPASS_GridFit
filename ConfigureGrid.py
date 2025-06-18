
import h5py

import pandas as pd
import numpy as np
from scipy import integrate as integr

from astropy.table import Table, vstack


class ConfigureGrid():
    
    "Class to configure a BPASS grid hdf5 file for use with pimodeling by Brinchmann"
    
    def __init__(self, grid, params={'Z':[], 'U':[], 'xsi':[], 'nH':[], 'age':[], 'CO':[]}, columns=[]):
        
        """
        BPASS grid should probe different parameters:
        
        metallicities - Z (not in terms of solar metallicity!)
        ionizaition parameters log(U) - U
        dust-to-metal ratio ξ - xsi
        density log(n_H) - nH
        age in Myr - age
        CO fraction - CO
        
        grid - path to the BPASS grid hdf5 file
        params - dictionary with the parameters to probe
        columns - names of the columns in the grid including line names
        """
        
        self.grid = grid
        self.Z = params['Z']
        self.U = params['U']
        self.xsi = params['xsi']
        self.nH = params['nH']
        self.age = params['age']
        self.CO = params['CO']
        self.ages = [f'{y:.2e}' for y in self.age]
        self.columns = columns

    def get_data(self, g1, g2, g3, g4, g5):

        """
        Get data from the BPASS grid for a given CO, logU, xsi, nH, and age.
        """
        
        #g1 - CO fraction
        #g2 - log(U)
        #g3 - dust-to-metal ratio ξ
        #g4 - density log(n_H)
        #g5 - age
        
        g1 = 'CO_'+str(g1).replace('0.', '')
        g2 = str(g2)
        g3 = str(g3)
        g4 = str(g4)
        
        with h5py.File(self.grid, 'r') as f:
            name = list(f[g1][g2][g3][g4][g5].keys())
            data = f[g1][g2][g3][g4][g5][name[0]][:]
        return data
    
    def load_ages(self, CO, logU, xsi, nH):
        
        """
        load data for all ages for a given CO, logU, xsi, nH
        """
        
        d = []
        for a in self.ages:
            data = self.get_data(CO, logU, xsi, nH, a)
            d.append(data)
        self.d = d
        return d
    
    def make_df(self, d):
        
        """
        make a dataframe from the loaded data
        """
        
        c = self.columns.copy()
        c.remove('n_H')
        c.extend(['age'])
        df_tot = pd.DataFrame(columns=c)
        for i in range(len(d)):
            df = pd.DataFrame(d[i], columns=self.columns)
            df['age'] = self.age[i]
            df.drop(columns=['n_H'], inplace=True)
            df_tot = pd.concat([df_tot, df], ignore_index=True)

        return df_tot

    def integrate_age(self, df, CO, logU, xsi, nH):
        
        """
        integrate the dataframe over age for a given CO, logU, xsi, nH
        create a new dataframe with the integrated values
        """
    
        df_integrated = pd.DataFrame(columns=df.keys().drop('age'))
        df_integrated['Z'] = np.unique(df['Z'])
        df_integrated['CO'] = [CO]*len(df_integrated)
        df_integrated['LOGU'] = [logU]*len(df_integrated)
        df_integrated['XI'] = [xsi]*len(df_integrated)
        df_integrated['NH'] = [nH]*len(df_integrated)
        
        #integrate over age
        for k in np.unique(df['Z']):
            for l in df_integrated.keys():
                if l not in ['Z', 'CO', 'LOGU', 'XI', 'NH']:
                    df_integrated[l].loc[df_integrated['Z']==k] = integr.trapezoid(df[l].loc[df['Z']==k].to_numpy(), df['age'].loc[df['Z']==k].to_numpy())
        df_integrated.astype(np.float64)
        
        return df_integrated
    
    def load_grid(self):
        
        """
        Load the grid data from the HDF5 file.
        """
        
        c = self.columns.copy()
        c = c.extend(['CO', 'LOGU', 'XI', 'NH'])
        df = pd.DataFrame(columns=c)

        for co in self.CO:
            for logu in self.U:
                for x in self.xsi:
                    for n in self.nH:
                        d = self.load_ages(co, logu, x, n)
                        df_d = self.make_df(d)
                        df_int = self.integrate_age(df_d, co, logu, x, n)
                        df = pd.concat([df, df_int], ignore_index=True)
                        
        df['s2_6716'] = df['s2_6716_6731']/3.96
        df['s2_6731'] = df['s2_6716']*2.96
        
        df = df.astype(np.float64)
        
        self.df = df

        return self.df

    def grid_setup(self):
        """
        Set up the grid for pifmodel and pifit
        """
        
        tables = []

        for j in range(len(self.Z)):
            tab = self.df[self.df['Z'] == self.Z[j]]
            tab.rename(columns={'n_H':'NH', 'o1_6300':'OI6300', 'o2_3727':'OII3727', 'o3_1666':'OIII1666', 'o3_4959':'OIII4959', 'o3_5007':'OIII5007', 's2_6716_6731':'SII6716', 'n2_6548':'NII6548', 'n2_6584':'NII6584', 'n5_1243':'NII1240', 'he1_3965':'HeI3965', 'he1_4471':'HeI4471', 'he1_5876':'HeI5876', 'he1_6678':'HeI6678', 'he2_1640':'HeII1640', 'he2_4686':'HeII4686', 'c3_1909':'CIII1909', 'c4_1551':'CIV1551', 'si3_1892':'SiIII1892', 'he1_1083':'HeI1083', 'ha':'HA', 'hb':'HB', 'Z':'ZMET', 'Age':'AGE', 's2_6716':'SII6717', 's2_6731':'SII6731'}, inplace=True)
            tab['LOGZ']=np.log10(tab['ZMET'])
            tab = Table.from_pandas(tab)
            
            tables.append(tab)
        
        self.table = vstack(tables)

    def save_grid(self, filename='grid.fits'):
        """
        Save the grid to a FITS file.
        """
        self.grid_setup()
        self.table.write(filename, format='fits', overwrite=True)
        print(f"Grid saved to {filename}")