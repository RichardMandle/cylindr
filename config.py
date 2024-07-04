import argparse

# Config class - just argument parsing
class Config:
    def __init__(self):
        self.parser = argparse.ArgumentParser(description='Calculate Cylindrical Distribution Function (CDF) From an MD Trajectory')
        self._add_arguments()
        self.args = self.parser.parse_args()

    def _add_arguments(self):
        # input and output file names
        self.parser.add_argument('-top', '--topology', default='', type=str, help='Input topology filename')
        self.parser.add_argument('-traj', '--trajectory', default='', type=str, help='Input trajectory filename')
        self.parser.add_argument('-o', '--output_name', default='', type=str, help='name used for output CDF plots')
        self.parser.add_argument('-replot', '--replot', action='store_true', help='Replot from saved data; specify pickle file to read with -o flag.')

        # atom/molecule/residue selection arguments
        self.parser.add_argument('-mode', '--selection_mode', default='default', type=str, help='default is to look at distances between centres of mass (COM). "-mode hybrid" computes distances between the COM and neighbouring atoms, specified according to -sel and -selname/-selelement')
        self.parser.add_argument('-sel', '--selection', default='com', type=str, help='Selection type; defaults to centre of mass (COM). selecting "name" or "element" requires additional -selname or -selelement input.')
        self.parser.add_argument('-selname', '--selection_name', default='', type=str, help='Used when -sel is "name"; pass atom name(s) as string, e.g. C1, O4, N9 etc.')
        self.parser.add_argument('-selelement', '--selection_element', default='', type=str, help='Used when -sel is "element"; pass element name(s) as string, e.g. H, C, O, N etc.')

        # parameters for calculating the CDF
        self.parser.add_argument('-l', '--cutoff_length', default=40, type=int, help='cylindrical shell length cutoff in angstroms (default = 40)')
        self.parser.add_argument('-r', '--cutoff_radius', default=15, type=int, help='cylindrical shell radial cutoff in angstroms (default = 15)')
        self.parser.add_argument('-b', '--first_frame', default=0, type=int, help='frame to start at (default = 0)')
        self.parser.add_argument('-e', '--end_frame', default=-1, type=int, help='frame to end at (default = 1)')
        self.parser.add_argument('-res', '--res', default=4, type=int, help='spatial resolution of integration grid in points per Angstrom (default = 4)')
        self.parser.add_argument('-ori', '--orientation', default='nem', type=str, help='specify orientation of cylinder length: "nem" = length along nematic director (default); "perp1" and "perp2" = perpendicular to nematic director (i.e. if nem=x, perp1=y and perp2=z); "user" = user supplied vector, please use the "-vec" flag and provide a 3vector (e.g. "1,0,0" = X axis)')
        self.parser.add_argument('-vec', '--orientation_vector', default='0,0,0', type=str, help='vector describing the orientation of the cylinder length; i.e. 1,0,0 = X, 0,1,0 = Y, 0,0,1= Z. Defaults to being along the nematic director')

        # output options: save data, plot formatting
        self.parser.add_argument('-save', '--save', action='store_true', help='save intermolecular distances as a compressed .npz file for reading by lookup.py. Default is no, because the files are big, but specify "yes" if you want them')
        self.parser.add_argument('-log', '--log_scale', action='store_true', help='If called, plot CDF using a log scale')
        self.parser.add_argument('-min', '--mini', default='0', type=float, help='specify the minimum value to use in CDF plot (default = 0; use min of dataset)')
        self.parser.add_argument('-max', '--maxi', default='0', type=float, help='specify the maximum value to use in CDF plot (default = 0; use max of dataset)')
        self.parser.add_argument('-cmap', '--cmap', default='magma', type=str, help='cmap to use in plotting (default = magma)')
        self.parser.add_argument('-interp', '--interp', default='bilinear', type=str, help='interpolation to use in plot (default = bilinear')
        self.parser.add_argument('-font_size', '--font_size', default=20, type=int, help='Font size for the plot')
        self.parser.add_argument('-fig_size', '--figure_size', default= 6, type=float, help='Figure size for the plot (height; inches?)')