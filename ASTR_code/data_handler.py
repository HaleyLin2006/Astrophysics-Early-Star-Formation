import numpy as np
import os

class data():
    def __init__(self):
        return
    
    def save(self, param=[], radius=[], r_grid=[], density=[], v_r=[], v_theta=[], v_phi=[]):
        directory = self._get_directory(param)

        if not os.path.exists(directory):
            print('Create new file...')
            np.savez(directory)

        np.savez(directory, radius=radius, r_grid=r_grid, density=density, v_r=v_r, v_theta=v_theta, v_phi=v_phi)

    def load(self, param=[]):
        directory = self._get_directory(param)

        if not os.path.exists(directory):
            print("Data does not exsit")
            return

        file = np.load(directory)
        
        radius  = file['radius']
        r_grid  = file['r_grid']
        density = file['density']
        v_r     = file['v_r']
        v_theta = file['v_theta']
        v_phi   = file['v_phi']

        return radius, r_grid, density, v_r, v_theta, v_phi
    
    def delete(self, param=[]):
        directory = self._get_directory(param)

        if not os.path.exists(directory):
            print("Data does not exsit")
            return
        
        print('Deleting file...')
        os.remove(directory)

    def _get_directory(self, param=[]):
        directory = 'data/'
        for i in param:
            directory += f'_{i}'
        
        return directory+'.npz'