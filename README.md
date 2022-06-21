# _cylindr_
_cylindr_ computes a cylindrical pair correlation function (CDF) from an MD simulation, this is a useful tool for identifying and visualising pairing modes in (polar) anisotropic fluids.

_cylindr_ looks for neighbouring molecules (or particles) in cylindrical shells by their centre-of-mass. Following this, neighbour distances are binned into a 2d histogram which is normalised against the volume of the cylindrical shells and also the average number of neighbours. This returns a the CDF plot, as shown below.  Trajectory and topology reading is performed using _MDtraj_, enabling the use of _cylindr_ with a number of different MD engines (tested with Gromacs 2019.2).

The default behaviour is for the length of the cylinder to be oriented along the vector describing the average orientation of the molecules/particles (the director), however the user can supply any 3-vector (_via_ the -ori user and -vec command line option) for orientation of the cylinder length. The cut off length (-L) and radius (-R) can be varied, as can the resolution (-res).

_cylindr_ optionally saves intermolecular distances for each frame in the trajectory as a compressed .npz file; a companion tool - _lookup.py_ is provided for accessing this data. The user supplies a lower and upper limits for length (-L) and radius (-R) cut offs to _lookup.py_, which returns the frame number and residue indices of any pairs of molecules that satisfy these distance constraints. This permits visualisation of pairs, but also lifetime studies and so on.

# Requirements
* MDtraj (tested on 1.98)
* numpy (tested using 1.22)
* Matplotlib (tested using 3.5.1)
* tqdm
* pickle

* some MD data (tested using Gromacs .trr and .gro files)

# Returns
* A cylindrical pair correlation function of the whole MD trajectory, e.g.:
![github](https://user-images.githubusercontent.com/101199234/174793536-b49d7abf-e66e-4a04-8ee7-b7b1ad93aa3a.png)


# Input Options
* -traj (required): trajectory file to be read by _cylindr_
* -top (required): topology information for -traj
* -O (required): output file name

* -mode (optional): choose between inter-COM distances (default) or COM-selection distances (-mode=hybrid); hybrid requires a valid selection for SEL to be specified with -sel flag.
* -sel (optional): selection mode; choose between molecular com (centre-of-mass), atom name (name) or element type (element)
* -selname (reqd. if -sel=name): used when -sel=name; pass atom names as string (e.g. Co1, Mn0 etc.)
* -selelement (reqd. if -sel=element): used when -sel=element; pass element names for selections as a string (e.g. H, C, Li).

* -L (optional; default=40 Å): cylindrical shell length cutoff in Å.
* -R (optional; default=40 Å): cylindrical shell radial cutoff in Å.
* -b (optional): first frame of -traj to read, default is 0.
* -e (optional): final frame of -traj to read, default is _end_.
* -res (optional; default = 4): spatial resolutuon in steps-per-angstrom.
* -ori (optional): orientation of cylinder length along nematic director (nem, default); perpendicular to director (perp1, perp2) or along user supplied vector (user, requires -vec)
* -vec (optional): 3vector for cylinder orientation; default is to orient along simulation principal orientation axis.

* -save (optional; default = no): save intermolecular distances as a .npz file for later reading by _lookup.py_
* -log (optional; default = no): Use log scale in CDF plot.
* -min (optional; default = 0): Specify minimum value used in CDF plot.
*  -max (optional): Specify maximum value used in CDF plot; default is to determine automatically.


# Usage Examples
_cylindr_
* python cylindr.py -traj traj.trr -top conf.gro                        - vanilla; computes the CDF with all default options
* python cylindr.py -traj traj.trr -top conf.gro -res 2                 - uses a resolution of 2 shells per angstrom (default = 4)
* python cylindr.py -O foo-traj traj.trr -top conf.gro -R 60            - saves output.npz as 'foo' and uses a radial cutoff of 60 angstroms.
* python cylindr.py -traj traj.trr -top conf.gro -L 20 -b 10000         - uses a length cutoff of 20 ansgtroms and starts from frame 10000
* python cylindr.py -traj traj.trr -top conf.gro -e 2500 -ori user -vec 0,1,1     - stop at frame 2500 and orient the cylinder length along 0,1,1
* python cylindr.py -traj traj.trr -top conf.gro -sel name -selname F1  - Only atoms whose name is F1 are used in the CDF analysis
* python cylindr.py -sel element -selelement H                          - Selects by element, in this case hydrogen
* python cylindr.py -mode hybrid -sel name -selname O6                  - hybrid mode - computes distance between COM of each molecule in turn and all atoms of name "O6". can use -sel element and -selelement in hybrid mode also.

_plotting_ and _output_
* python cylindr.py -traj traj.trr -top conf.gro -min 1 -max 5          - -min and -max set the minimum and maximum values of the colormap in the CDF plot
* python cylindr.py -traj traj.xtc -top conf.gro -log yes               - use a log normalised colormap in the CDF plot
* python cylindr.py -traj dump.traj -top conf.gro -save yes             - -save yes dumps distance matrix data as a .npz file, can be read by _lookup.py_ to return pair indicies over a specified length range for visualisation/lifetime anlysis
 
_lookup_
* python lookup.py -i foo -L 20,21.2 -R 6,6.5                           - return pair indicies in the length range 20 to 21 angstrom, radius range 6 to 6.5 angstrom.
* python lookup.py -i foo -L -7,-6 -R 5,6, -o foo2                      - return pair indicies in the length range -7 to -6 angstrom, radius range 5 to 6 angstrom, save output as foo2
# Questions, problems?
Make a Github issue :). Please be as clear and descriptive as possible, and please do feel free to reach out in person: (r[DOT]mandle[AT]leeds[DOT]ac[DOT]uk)
