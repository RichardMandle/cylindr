import mdtraj as md
import numpy as np
from tqdm import tqdm

def rotation_matrix_from_vectors(vec1, vec2):
    """Return rotation matrix that rotates vec1 onto vec2."""
    a = np.array(vec1, dtype=float).reshape(3)
    b = np.array(vec2, dtype=float).reshape(3)

    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na == 0.0 or nb == 0.0:
        return np.eye(3)
    a /= na
    b /= nb

    c = float(np.dot(a, b))

    # Parallel: no rotation
    if np.isclose(c, 1.0):
        return np.eye(3)

    # Anti-parallel: 180 rotation about any axis
    if np.isclose(c, -1.0):
        if abs(a[0]) < 0.9:
            perp = np.array([1.0, 0.0, 0.0])
        else:
            perp = np.array([0.0, 1.0, 0.0])
        v = np.cross(a, perp)
        v /= np.linalg.norm(v)
        x, y, z = v
        K = np.array([[0, -z,  y],
                      [z,  0, -x],
                      [-y, x,  0]])
        return np.eye(3) + 2.0 * (K @ K)

    # General case
    v = np.cross(a, b)
    s = np.linalg.norm(v)         
    v /= s                         
    x, y, z = v
    K = np.array([[0, -z,  y],
                  [z,  0, -x],
                  [-y, x,  0]])

    # rodrigues' rotation formula with unit axis:
    return np.eye(3) + K * s + K @ K * (1.0 - c)
  
class TrajectoryProcessor:
    '''
    TrajectoryProcessor class for... processing trajectories with mdtraj
    '''
    def __init__(self, config):
        self.config = config
        self.traj = None
        self.com = None
        self.sel = None

        self.load_trajectory()
        self.calculate_com()
        self.apply_orientation()  # handle different orientation
        
        # we should only bother to build sel if we are in hybrid mode 
        if self.config.args.selection_mode == 'hybrid':
            self.calculate_sel()
            
            
    def load_trajectory(self):
        print("Loading trajectory...")
        self.traj = md.load_trr(self.config.args.trajectory, top=self.config.args.topology)
        self.traj = self.traj[self.config.args.first_frame:self.config.args.end_frame]
        
    def apply_orientation(self):
        '''
        Here we apply the reorientation of the trajectory
        '''
        args = self.config.args
        ori = args.orientation.lower()

        if ori == 'none':
            print('No rotational alignment requested.')
            return

        if self.com is None:
            raise RuntimeError('COM must be computed before orientation.')

        traj = self.traj
        n_frames = traj.n_frames
        n_res = traj.n_residues
        n_atoms = traj.n_atoms

        mol_size = n_atoms // n_res
        chain_indices = [[i + k for k in range(mol_size)]
                         for i in range(0, n_atoms, mol_size)]

        print('Computing global nematic director...')
        directors = md.compute_directors(traj, chain_indices)  
        D = directors.mean(axis=(0, 1))
        
        if np.linalg.norm(D) == 0.0:
            print('Warning: nematic director magnitude is zero; skipping rotation.')
            return
            
        D /= np.linalg.norm(D)

        tmp = np.array([0.0, 0.0, 1.0])
        
        if abs(np.dot(tmp, D)) > 0.9:
            tmp = np.array([0.0, 1.0, 0.0])
            
        e1 = D
        e2 = np.cross(e1, tmp); e2 /= np.linalg.norm(e2)
        e3 = np.cross(e1, e2)

        # here we are deciding which physical axis to align with +x
        if ori == 'nem':
            src = e1
            label = 'nematic director'
        elif ori == 'perp1':
            src = e2
            label = 'perpendicular direction perp1'
        elif ori == 'perp2':
            src = e3
            label = 'perpendicular direction perp2'
        elif ori == 'user':
            vec = np.array(args.orientation_vector.split(','), dtype=float)
            if np.allclose(vec, 0.0):
                print('User orientation vector is zero; skipping rotation.')
                return
            src = vec / np.linalg.norm(vec)
            label = f'user vector {src}'
        else:
            print(f'Unknown orientation "{ori}", skipping rotation.')
            return

        print(f'Rotating so that {label} aligns with +x axis.')

        # we'll do a single global rotation so we align with [1,0,0] which is the length of our cylinder.
        R = rotation_matrix_from_vectors(src, np.array([1.0, 0.0, 0.0]))

        self.com = np.einsum('ij,taj->tai', R, self.com)
        if self.sel is not None:
            self.sel = np.einsum('ij,taj->tai', R, self.sel)
            
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
    
    def calculate_sel(self):
        """
        build the self.sel for hybrid mode:
        distances from selected sites (name/element) to COMs.
        """
        args = self.config.args
        traj = self.traj
        top = traj.topology

        print('Identifying selection coordinates for hybrid COM/selection CDF...')
        self.sel = np.zeros_like(self.com)

        for r in top.residues:
            ridx = r.index

            if args.selection == 'com':
                # hybrid+com is effectively default mode;
                self.sel[:, ridx, :] = self.com[:, ridx, :]

            elif args.selection == 'name':
                names = [s.strip() for s in args.selection_name.split(',') if s.strip()]
                atom_idx = [a.index for a in r.atoms if a.name in names]

                if len(atom_idx) == 0:
                    raise ValueError(f"Residue {ridx}: no atoms match -selname {args.selection_name}")

                sub = traj.atom_slice(atom_idx)
                com_sel = md.compute_center_of_mass(sub) * 10.0  # A
                self.sel[:, ridx, :] = com_sel

            elif args.selection == 'element':
                elems = [s.strip() for s in args.selection_element.split(',') if s.strip()]
                atom_idx = [a.index for a in r.atoms
                            if (a.element is not None and a.element.symbol in elems)]

                if len(atom_idx) == 0:
                    raise ValueError(f"Residue {ridx}: no atoms match -selelement {args.selection_element}")

                sub = traj.atom_slice(atom_idx)
                com_sel = md.compute_center_of_mass(sub) * 10.0  # A
                self.sel[:, ridx, :] = com_sel

            else:
                raise ValueError('Hybrid mode requires -sel com|name|element.')
        print('Selection COMs computed for hybrid mode.')