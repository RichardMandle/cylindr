import mdtraj as md
import numpy as np
from tqdm import tqdm

# TrajectoryProcessor class for... processing trajectories with mdtraj
class TrajectoryProcessor:
    def __init__(self, config):
        self.config = config
        self.traj = None
        self.com = None
        self.sel = None
        self.load_trajectory()
        self.calculate_com()

    def load_trajectory(self):
        print("Loading trajectory...")
        self.traj = md.load_trr(self.config.args.trajectory, top=self.config.args.topology)
        self.traj = self.traj[self.config.args.first_frame:self.config.args.end_frame]

    def calculate_com(self):
        print('Calculating centres of mass...')
        self.com = np.zeros((self.traj.n_frames, self.traj.n_residues, 3))
        for n in tqdm(range(self.traj.n_residues), unit=' Molecules'):
            if self.config.args.selection == 'com':
                self.com[:, n, :] = md.compute_center_of_mass(self.traj, select='resid ' + str(n)) * 10
            elif self.config.args.selection == 'name':
                self.com[:, n, :] = md.compute_center_of_mass(self.traj, select='name ' + self.config.args.selection_name + ' and resid ' + str(n)) * 10
            elif self.config.args.selection == 'element':
                self.com[:, n, :] = md.compute_center_of_mass(self.traj, select='element ' + self.config.args.selection_element + ' and resid ' + str(n)) * 10

        if self.config.args.selection_mode == 'hybrid':
            self.sel = np.zeros((self.traj.n_frames, self.traj.n_residues, 3))
            for n in tqdm(range(self.traj.n_residues), unit=' Selections'):
                if self.config.args.selection == 'com':
                    self.sel = self.com
                elif self.config.args.selection == 'name':
                    self.sel[:, n, :] = md.compute_center_of_mass(self.traj, select='name ' + self.config.args.selection_name + ' and resid ' + str(n)) * 10
                elif self.config.args.selection == 'element':
                    self.sel[:, n, :] = md.compute_center_of_mass(self.traj, select='element ' + self.config.args.selection_element + ' and resid ' + str(n)) * 10
