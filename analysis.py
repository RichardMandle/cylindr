import numpy as np
from numba import jit, njit, prange
from tqdm import tqdm
from fast_histogram import histogram2d # trial

class CDFAnalysis:
    def __init__(self, config, traj_processor):
        self.config = config
        self.traj_processor = traj_processor
        self.dist_length = None
        self.dist_radius = None
        self.cdf = None
        self.volume = None

    def calculate_distances(self):
        print('Calculating inter-com distances...')
        com = self.traj_processor.com
        self.dist_length = np.zeros((com.shape[0], com.shape[1], com.shape[1]), dtype=np.float32)
        self.dist_radius = np.zeros((com.shape[0], com.shape[1], com.shape[1]), dtype=np.float32)

        if self.config.args.selection_mode == 'default':
            for n in tqdm(range(com.shape[0]), unit=' timestep'):
                self.dist_length[n, :, :] = distance_matrix(com[n, :, 0])
                self.dist_radius[n, :, :] = np.sqrt(distance_matrix(com[n, :, 1])**2 + distance_matrix(com[n, :, 2])**2)
            
        elif self.config.args.selection_mode == 'hybrid':
            sel = self.traj_processor.sel
            for n in tqdm(range(com.shape[0]), unit=' timestep'):
                self.dist_length[n, :, :] = dual_matrix_distance(sel[n, :, 0], com[n, :, 0])
                self.dist_radius[n, :, :] = np.sqrt(dual_matrix_distance(sel[n, :, 1], com[n, :, 1])**2 + dual_matrix_distance(sel[n, :, 2], com[n, :, 2])**2)
        
        # set diagonals to large values suppresses the self interaction potential swamping the CDF.
        indices = np.arange(com.shape[1]) # done out of loop for speed (its maybe 20% faster)
        self.dist_length[:, indices, indices] = 10000
        self.dist_radius[:, indices, indices] = 10000
        
        if not self.config.args.symmetric:
            #use np.triu to set the lower part of the com matrix to some big number so its not head-tail symmetric (or not necesarily!)
            print(f'\n !! Ignoring diagonal symmetry in COM matrix !! \n') # temp debug while testing
            triu_indices = np.triu_indices(com.shape[1], k=1)
            self.dist_length[:, triu_indices[1], triu_indices[0]] = 10000 
            self.dist_radius[:, triu_indices[1], triu_indices[0]] = 10000
        
    def cylindrical_pcf(self):
        # make cdf (= a cylindrcal pair correlation function)
        length_cut = self.config.args.cutoff_length
        radius_cut = self.config.args.cutoff_radius
        res = self.config.args.res

        bins = ((2 * length_cut * res) - 1, (radius_cut * res) - 1)
        self.cdf = np.zeros([self.dist_length.shape[0], bins[0], bins[1]], dtype=np.float32)

        print('Performing CDF Analysis...')
        for n in tqdm(range(self.dist_length.shape[0]), unit=' timestep'):
            # using fast_histogram.histogram2d gives a big speedup over np.histogram2d or even a custom numba-jit function
            self.cdf[n, :, :] = histogram2d(
                self.dist_length[n, :, :].ravel(), self.dist_radius[n, :, :].ravel(),
                bins=[(2 * length_cut * res) - 1, (radius_cut * res) - 1],
                range=[(-length_cut, length_cut), (0, radius_cut)]
            )

        self.volume = np.zeros((bins[0], bins[1]))
        dL = 2 * length_cut / bins[0]
        for n in range(bins[0]):
            for m in range(bins[1]):
                r_in = m / res
                r_out = (m + 1) / res
                self.volume[n, m] = np.pi * (r_out**2 - r_in**2) * dL

    def save_additional_data(self, filename):
        # save additional data to read with lookup.py
        np.savez_compressed(filename, dist_length=self.dist_length, dist_radius=self.dist_radius)
        print(f"Saved additional data to {filename}") 
        
# custom numba jitted functions:
@jit(nopython=True, parallel=True)
def distance_matrix(a):
    # return the distance matrix for all centres-of-mass pairs
    x = np.ascontiguousarray(a).reshape((len(a), 1))    # use ascontiguiusarray rather than np.reshape for numba jit compatibility
    dist_matrix = x - x.T
    return dist_matrix

@jit(nopython=True, parallel=True)
def dual_matrix_distance(a, b):
    # return the distance matrix for all possible permutations of a/b 
    x = np.zeros((len(b), len(a)), dtype=np.float32)
    for n in prange(len(b)):
        x[n, :] = b.T - a[n]
    return x
