# ! /usr/bin/env python

# cylindr - main code
# Dr. R. J. Mandle - University of Leeds, 2022


import mdtraj as md
import numpy as np
import matplotlib as mplot
from numpy import savetxt
import time
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import colors
from matplotlib.colors import LogNorm
import os
import platform
import argparse
from tqdm import tqdm
import pickle

## IDEA 1
# Would be nice to be able to do a COM-vs-atom map; i.e. make a CDF of the COM-to-atom.name distances.
# this probably needs a different version of distance_matrix that can handle two arrays...
# something like the distance between COM[:,n,:] and ATOM[:,:,:], but excluding ATOM[:,n,:]... how?
#
# IDEA 2
# Would be nice to be able to produce a map of average orientation vs distance(s); 
# like G1(r) and G2(r), would be like g1(r,h), g2(r,h) I guess...
def initialize():
    parser = argparse.ArgumentParser(description='Calculate Cylindrical Distribution Function (CDF) From an MD Trajectory')

    # input and output file names
    parser.add_argument('-top', '--topology', default='', type=str, help='Input topolgy filename')
    parser.add_argument('-traj', '--trajectory', default='', type=str, help='Input trajectory filename')
    parser.add_argument('-o', '--output_name', default='', type=str, help='name used for output CDF plots')

    # atom/molecule/residue selection arguments
    parser.add_argument('-mode','--selection_mode', default='default', type=str, help='default is to look at distances between centres of mass (COM). "-mode hybrid" computes distances between the COM and neighbouring atoms, specified according to -sel and -selname/-selelement')
    parser.add_argument('-sel','--selection', default='com', type=str, help='Selection type; defaults to centre of mass (COM). selecting "name" or "element" requires additional -selname or -selelement input.')
    parser.add_argument('-selname','--selection_name', default='', type=str, help='Used when -sel is "name"; pass atom name(s) as string, e.g. C1, O4, N9 etc.')
    parser.add_argument('-selelement','--selection_element', default='', type=str, help='Used when -sel is "element"; pass element name(s) as string, e.g. H, C, O, N etc.')

    # parameters for calculating the CDF
    parser.add_argument('-L', '--Cutoff_L', default=40, type=int, help='cylindrical shell length cutoff in angstroms')
    parser.add_argument('-R', '--Cutoff_R', default=15, type=int, help='cylindrical shell radial cutoff in angstroms')
    parser.add_argument('-b', '--first_frame', default=0, type=int, help='frame to start at')
    parser.add_argument('-e', '--end_frame', default=-1, type=int, help='frame to end at')
    parser.add_argument('-res', '--res', default=4, type=int, help='spatial resolution of integration grid in points per Angstrom')
    parser.add_argument('-ori', '--orientation', default='nem', type=str, help='specify orientation of cylinder length: "nem" = length along nematic director (default); "perp1" and "perp2" = perpendicular to nematic director (i.e. if nem=x, perp1=y and perp2=z); "user" = user supplied vector, please use the "-vec" flag and provide a 3vector (e.g. "1,0,0" = X axis)')
    parser.add_argument('-vec', '--orientation_vector', default='0,0,0', type=str, help='vector describing the orientation of the cylinder length; i.e. 1,0,0 = X, 0,1,0 = Y, 0,0,1= Z. Defaults to being along the nematic director')

    # output options: save data, plot formatting
    parser.add_argument('-save', '--save', default='no', type=str, help='save intermolecular distances as a compressed .npz file for reading by lookup.py. Default is no, because the files are big, but specify "yes" if you want them')
    parser.add_argument('-log', '--log_scale', default='no', type=str, help='plot CDF using a log scale')
    parser.add_argument('-min', '--mini', default ='0', type=float, help='specify the minimium value to use in CDF plot')
    parser.add_argument('-max', '--maxi', default ='0', type=float, help='specify the maximum value to use in CDF plot')

    return parser


def rotation_matrix_from_vectors(vec1, vec2):

    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))

    return rotation_matrix


def distance_matrix(a):

    x = np.reshape(a, (len(a), 1))

    return x - x.transpose()
    

def dual_matrix_distance(a, b):

    x = np.zeros([np.size(b),np.size(a)])

    for n in range(0,np.size(b)):
        x[n,:] = b.transpose() - a[n]
    
    # setting the diagonal elements of to a large value is a quick and dirty way to remove the intermolecular COM-SEL peak. 
    
    np.fill_diagonal(x,10000)   
    return x


def cylindrical_PCF(L_cut, R_cut, Res,Ldist,Rdist): # main CDF engine
    
    Bins = [np.array(range(- L_cut * Res, L_cut * Res)) / Res, np.array(range(R_cut * Res)) / Res]
    CDF = np.zeros([np.size(Ldist, 0), (2 * L_cut * Res) - 1, (R_cut * Res) - 1])

    print('Performing CDF Analysis...')

    for n in tqdm(range(0, np.size(Ldist, 0)), unit=' timestep'):
    
        CDF[n, :, :], XEdges, YEdges = np.histogram2d(np.ravel(Ldist[n, :, :]),\
        np.ravel(Rdist[n, :, :]), bins=Bins, range=[[-L_cut, L_cut], [0, R_cut/2]])

    print('Computing cylindrical shell volumes')

    Vol = np.zeros((np.size(XEdges) - 1, np.size(YEdges) - 1))  # This is the shell volume we normalise against:

    for n in tqdm(range(0, np.size(XEdges) -1), unit=' shells'):

        for m in range(0, np.size(YEdges) - 1):
        
            Vol[n, m] = ((2 * np.pi * (YEdges[m + 1]) ** 2) - (2 * np.pi * (YEdges[m - 1]) ** 2)) * np.abs(XEdges[2] - XEdges[1])

    return CDF, Vol


args = initialize().parse_args()

print("running cylindr.py on", platform.system(), platform.release(), platform.version())
print("Please cite the paper: 10.26434/chemrxiv-2022-z9rxk")

# load supplied trajectory (-traj) and topology (-top) with mdtraj
print("Loading trajectory...")

traj = md.load_trr(args.trajectory, top=args.topology)  # load trajectory (.trr) and topology (.gro)
traj = traj[args.first_frame:args.end_frame]            # cut traj down if user supplied -b/-e arguments
Atoms = traj.n_atoms                                    # Atom Count
NMols = traj.n_residues                                 # Molecule/residue count

chain_indices = [[n+x for x in range(int(Atoms/NMols))] for n in range(0, Atoms, (int(Atoms/NMols)))]

# Calculate the cetre-of-mass (COM) for each molecule, then reorient as required

print('Calculating centres of mass...')
COM = np.zeros((traj.n_frames, traj.n_residues, 3))

for n in tqdm(range(0, NMols), unit=' Molecules'):

    if args.selection == 'com':
        COM[:, n, :] = md.compute_center_of_mass(traj, select='resid ' + str(n)) * 10  # compute centres of mass, scale from nm to angstrom.

    if args.selection == 'name':        
        COM[:, n, :] = md.compute_center_of_mass(traj, select='name '+ args.selection_name + ' and resid ' + str(n)) * 10  # compute centres of mass, scale from nm to angstrom.

    if args.selection == 'element':        
        COM[:, n, :] = md.compute_center_of_mass(traj, select='element '+ args.selection_element + ' and resid ' + str(n)) * 10  # compute centres of mass, scale from nm to angstrom.

# ... and if we are using hybrid mode, we need to compute another set of selections:

if args.selection_mode == 'hybrid':
    if args.selection == 'com': # throw a warning if -sel=com & -mode=hybrid, because this is the same as mode=default...    
        print('user called -mode "hybrid" but also -sel "com"; this is the same as -mode "default"')
        print('\n check this is what you intended')
    
    print('\n Identifying selection coordinates for COM/' + str(args.selection) +'  CDF analysis...')
    SEL = np.zeros((traj.n_frames, traj.n_residues, 3))

    for n in tqdm(range(0, NMols), unit=' Selections'):
    
        if args.selection == 'com':
            SEL = COM
        
        if args.selection == 'name':
            SEL[:, n, :] = md.compute_center_of_mass(traj, select='name ' + args.selection_name + ' and resid ' + str(n)) * 10          # compute centres of mass, scale from nm to angstrom.

        if args.selection == 'element':
            SEL[:, n, :] = md.compute_center_of_mass(traj, select='element '+ args.selection_element + ' and resid ' + str(n)) * 10     # compute centres of mass, scale from nm to angstrom.

# rotate COM-map - default is to orient nematic director along X

if args.orientation != 'user':
    if args.orientation == 'nem':
        NewDirector = [1, 0, 0]
        
    if args.orientation == 'perp1':
        NewDirector = [0, 1, 0]   
        
    if args.orientation == 'perp2':
        NewDirector = [0, 0, 1]   
        
    print('Computing nematic director...')
    Director = (np.mean((md.compute_directors(traj, chain_indices)), 1))
    Director = Director/np.linalg.norm(Director)
    print('Rotating Frames...')
    
    for n in tqdm(range(0, traj.n_frames), unit=' Frames'):         # perform rotation of COM

        rot = rotation_matrix_from_vectors([Director[n, :]], NewDirector)
        COM[n, :, :] = np.dot(COM[n, :, :], rot)
        
        if args.selection_mode == 'hybrid':                         # apply same transformation to SEL, if required by -mode "hybrid"
             SEL[n, :, :] = np.dot(SEL[n, :, :], rot)

if args.orientation == 'user':
    SkipOrient = 0

    Orientation_Vector = np.array((args.orientation_vector).split(','),dtype=float) 

    if np.all(Orientation_Vector == [1,0,0]) == True:               # special case #1; if we want to orient along [1,0,0] we don't need to rotate at all
        print('User vector is [1,0,0]; no rotation needed, as this is the default')
        SkipOrient = 1

    if np.all(Orientation_Vector == [0,0,0]) == True:               #special case #2; if user forgets to supply -vec, just go for [1,0,0] and print a note
        print('You didn\'t supply a vector to use, so no rotation will be performed')
        Orientation_Vector = np.array([1,0,0])
        SkipOrient = 1

    if SkipOrient != 1:                                             # if the user supplied a vector that isn't [1,0,0] then do rotation now
        print('Orienting along ' + str(Orientation_Vector))
        print(str(type(Orientation_Vector)))
        Director = [1, 0, 0]                                        # Here, we want to rotate the simulation from X towards the users vector.       

        print('Rotating Frames...')
        for n in tqdm(range(0, traj.n_frames), unit=' Frames'):     # perform rotation of COM

            rot = rotation_matrix_from_vectors(Director, Orientation_Vector)
            COM[n, :, :] = np.dot(COM[n, :, :], rot)

            if args.selection_mode == 'hybrid':                     # apply same transformation to SEL, if required by -mode "hybrid"
                SEL[n, :, :] = np.dot(SEL[n, :, :], rot)    

print('Calculating inter-COM distances...')

Dist_L = np.zeros((np.size(COM, 0), np.size(COM, 1), np.size(COM, 1)), dtype=np.float16)  
Dist_R = np.zeros((np.size(COM, 0), np.size(COM, 1), np.size(COM, 1)), dtype=np.float16)

if args.selection_mode == 'default':                           # default mode: use inter-COM distance    
    for n in tqdm(range(0, np.size(COM, 0)), unit=' timestep'):
        Dist_L[n, :, :] = distance_matrix(COM[n, :, 0])
        Dist_R[n, :, :] = ((distance_matrix(COM[n, :, 1]) ** 2) + (distance_matrix(COM[n, :, 2]) ** 2)) ** 0.5

if args.selection_mode == 'hybrid':                             # hybrid mode: selection-COM distances ignoring the intermolecular
    for n in tqdm(range(0, np.size(COM, 0)), unit=' timestep'):
        Dist_L[n, :, :] = dual_matrix_distance(SEL[n,:,0],COM[n, :, 0])
        Dist_R[n, :, :] = ((dual_matrix_distance(SEL[n,:,1],COM[n, :, 1]) ** 2) + (dual_matrix_distance(SEL[n,:,2],COM[n, :, 2]) ** 2)) ** 0.5

H, Vol = cylindrical_PCF(args.Cutoff_L, args.Cutoff_R, args.res, Dist_L, Dist_R)

# Plotting
print('Making CDF Plots...')

CDF = (np.sum(H, 0) / (Vol)) / np.mean(np.sum(H, 0) / (Vol))

# plotting engine
plt.rcParams.update({'font.size': 20})

plt.figure(figsize=(8, 10))

Mini = 0.0
Maxi = np.round(np.max(CDF),decimals=1)

if args.mini != 0:
    Mini = args.mini

if args.maxi != 0:
    Maxi = args.maxi

if args.log_scale != 'yes':     # Don't log scale
    MainPlot = plt.imshow((CDF[:,args.res:]), cmap='magma', vmin=str(Mini), vmax=str(Maxi), interpolation='bilinear', aspect='equal', extent=[1, np.size(CDF, 1), 0, np.size(CDF, 0)])

if args.log_scale == 'yes':     # Or, if requested, use log scale
    MainPlot = plt.imshow((CDF[:,args.res:]),cmap='magma', interpolation='bilinear', extent=[1, np.size(CDF, 1), 0, np.size(CDF, 0)], norm=colors.LogNorm(vmin=str(Mini+1), vmax=str(Maxi)))

XTicks = np.array(args.res * (np.flip(np.arange(args.Cutoff_R, 1, - 5))))
XTicks.astype(int)
XTickLabels = np.array((np.flip(np.arange(args.Cutoff_R, 1, - 5))))

YTicks = np.array(np.arange(0, args.Cutoff_L * 2 * args.res, 10 * args.res))
YTicks.astype(int)
YTickLabels = np.array(np.arange(-args.Cutoff_L, args.Cutoff_L, 10))

plt.xlabel('$\it{r}$ / Å')
plt.xticks(ticks=XTicks, labels=XTickLabels)
plt.ylabel('$\it{h}$ / Å')
plt.yticks(ticks=YTicks, labels=YTickLabels)
plt.colorbar(MainPlot, extend='max', shrink=0.7, label='$\it{g_{(h, r)}}$')  # set colourbar options

plt.grid(color='w', linestyle='-', linewidth=0.5)   # set grid options

LogStr = ''
if args.log_scale == 'yes': 
    LogStr ='_log'
    
# need to encode selection string in output filename...

# orientation is tricky, multiple options... this works OK...
if args.orientation == 'user':
    Out_Str_Ori = ('_user_vec_' + str(args.orientation_vector))
if args.orientation != 'user':
    Out_Str_Ori = str(args.orientation)
    
Out_String = (args.output_name + Out_Str_Ori + '_res' + str(args.res) + '_L' + str(args.Cutoff_L) + '_R' + str(args.Cutoff_R) + LogStr + '_sel' + args.selection + '_' + args.selection_name + args.selection_element + '_' + args.selection_mode)

# save as raster and vector images, save figure as pickle for later editing.
plt.savefig(Out_String + '_CDF.pdf')
plt.savefig(Out_String + '_CDF.png')
pickle.dump(MainPlot, open(Out_String + '_CDF.fig.pickle', 'wb'))  # save as pickle, we can edit this later if we want.
print('Plotting done :)')

# Save compressed distance data for later reriveal (if requested)
if args.save == 'yes':
    print('\n User requested data to be saved')
    print('Saving CDF Data...')
    np.savez_compressed(Out_String, Dist_L=Dist_L, Dist_R=Dist_R)

if args.save != 'yes':
    print('\n CDF data is NOT being saved, to save rerun with -save "yes"')
