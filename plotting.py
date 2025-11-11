import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import pickle
import json


class Plotting:
    def __init__(self, config, cdf_analysis=None):
        self.config = config
        self.cdf_analysis = cdf_analysis
        self.cdf = None

    def generate_cdf(self):
        # generate the cdf but averaging over all timesteps and normalizing against volume
        
        H = np.mean(self.cdf_analysis.cdf, axis=0)  
        vol = self.cdf_analysis.volume
        self.cdf = (H / vol) / np.mean(H / vol)     
        
        return self.cdf

    def generate_filename(self, log_scale=False, diameter = False):
        log_str = '_log' if log_scale else ''   
        dia_str = '_diam_' if diameter else ''
        out_string = (self.config.args.output_name + dia_str + log_str)
        return out_string

    def save_plot_data(self, filename):
        # save data required for replotting flag
         np.savez_compressed(filename, cdf=self.cdf, volume=self.cdf_analysis.volume, config=self.config.args.__dict__)

    def load_plot_data(self, filename):
        # handle loading data for replotting if -replot is called
        data = np.load(filename, allow_pickle=True)
        self.cdf = data['cdf']

    def save_config(self, filename):
        # write the current config to {filename}.txt so we can see what was done.
        with open(f"{filename}.txt", 'w') as f:
            json.dump(self.config.args.__dict__, f, indent=4)

    def plot_cdf(self):
        args = self.config.args
        plt.rcParams.update({'font.size': args.font_size})

        def get_data():
            base = self.cdf[:, args.res:]
            if getattr(args, "diameter", False):
                # mirror to show -R ... +R (full diameter)
                data = np.concatenate([np.fliplr(base), base], axis=1)
                return data, True
            else:
                return base, False

        def plot_and_save(log_scale: bool):
            data, is_diam = get_data()

            v_min = args.mini
            finite = data[np.isfinite(data)]
            if args.maxi != 0:
                v_max = args.maxi
            else:
                v_max = np.round(np.percentile(finite, 99), 1)

            ny, nx = data.shape
            y_min, y_max = -args.cutoff_length, args.cutoff_length
            if is_diam:
                x_min, x_max = -args.cutoff_radius, args.cutoff_radius
            else:
                x_min, x_max = 0.0, args.cutoff_radius

            extent = [x_min, x_max, y_min, y_max]

            plt.figure(figsize=(args.figure_size, args.figure_size),dpi = 300)

            if log_scale:
                norm = colors.LogNorm(vmin=max(v_min, 1e-10), vmax=v_max)
                main_plot = plt.imshow(
                    data, cmap=args.cmap, interpolation=args.interp,
                    aspect='equal', extent=extent, origin='lower', norm=norm
                )
            else:
                main_plot = plt.imshow(
                    data, cmap=args.cmap, interpolation=args.interp,
                    aspect='equal', extent=extent, origin='lower',
                    vmin=v_min, vmax=v_max
                )
            if is_diam:
                xticks = np.arange(-args.cutoff_radius, args.cutoff_radius + 0.1, 10)
                xlabel = r'$x$ / $\mathrm{\AA}$'
            else:
                xticks = np.arange(0.0, args.cutoff_radius + 0.1, 5)
                xlabel = r'$r$ / $\mathrm{\AA}$'

            yticks = np.arange(y_min, y_max + 0.1, 10)

            plt.xticks(xticks)
            plt.yticks(yticks)
            plt.xlabel(xlabel)
            plt.ylabel(r'$h$ / $\mathrm{\AA}$')
            plt.colorbar(main_plot, extend='max', shrink=0.7, label=r'$g(h,r)$')
            plt.grid(color='w', linestyle='-', linewidth=0.5, alpha=0.5)

            fname = self.generate_filename(log_scale)
            if is_diam:
                fname = fname.replace("_CDF", "") + "_diameter"

            plt.savefig(f"{fname}_CDF.pdf", bbox_inches='tight')
            plt.savefig(f"{fname}_CDF.png", bbox_inches='tight')
            pickle.dump(main_plot, open(f"{fname}_CDF.fig.pickle", 'wb'))
            plt.close()

        plot_and_save(log_scale=False)
        plot_and_save(log_scale=True)