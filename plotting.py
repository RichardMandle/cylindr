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
        H = np.mean(self.cdf_analysis.cdf, axis=0)  # should be average not cumulative...
        vol = self.cdf_analysis.volume
        self.cdf = (H / vol) / np.mean(H / vol)     # normalize against the mean density
        return self.cdf

    def generate_filename(self, log_scale=False):
        # simply put "_log" into a filename if needed. Does it need its own function?
        log_str = '_log' if log_scale else ''       
        out_string = (self.config.args.output_name + log_str)
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
        
        plt.rcParams.update({'font.size': self.config.args.font_size})
        plt.figure(figsize=(self.config.args.figure_size, self.config.args.figure_size)) # funny shape, but we'll tight_layout() it later on.

        def plot_and_save(log_scale):
            v_min = self.config.args.mini
            v_max = self.config.args.maxi if self.config.args.maxi != 0 else np.round(np.percentile(self.cdf, 99), decimals=1)
            
            
            if log_scale:
                main_plot = plt.imshow(
                    self.cdf[:, self.config.args.res:], cmap=self.config.args.cmap,
                    interpolation=self.config.args.interp, aspect='equal',
                    extent=[1, np.size(self.cdf, 1), 0, np.size(self.cdf, 0)],
                    norm=colors.LogNorm(vmin=max(v_min, 1e-10), vmax=v_max)
                )

            else:
                main_plot = plt.imshow(
                    self.cdf[:, self.config.args.res:], cmap=self.config.args.cmap, vmin=v_min, vmax=v_max,
                    interpolation=self.config.args.interp, aspect='equal', extent=[1, np.size(self.cdf, 1), 0, np.size(self.cdf, 0)]
                )


            x_ticks = np.array(self.config.args.res * (np.flip(np.arange(self.config.args.cutoff_radius, 1, -5))))
            x_tick_labels = np.array((np.flip(np.arange(self.config.args.cutoff_radius, 1, -5))))
            y_ticks = np.array(np.arange(0, self.config.args.cutoff_length * 2 * self.config.args.res, 10 * self.config.args.res))
            y_tick_labels = np.array(np.arange(-self.config.args.cutoff_length, self.config.args.cutoff_length, 10))

            plt.xlabel('$\it{r}$ / Å')
            plt.xticks(ticks=x_ticks, labels=x_tick_labels)
            plt.ylabel('$\it{h}$ / Å')
            plt.yticks(ticks=y_ticks, labels=y_tick_labels)
            plt.colorbar(main_plot, extend='max', shrink=0.7, label='$\it{g_{(h, r)}}$')
            plt.grid(color='w', linestyle='-', linewidth=0.5)

            filename = self.generate_filename(log_scale)
            plt.savefig(f"{filename}_CDF.pdf", bbox_inches='tight') # told you
            plt.savefig(f"{filename}_CDF.png", bbox_inches='tight')
            pickle.dump(main_plot, open(f"{filename}_CDF.fig.pickle", 'wb'))  # save as pickle for later reloading
            plt.close()

        # Plot and save both linear and log scale plots
        plot_and_save(log_scale=False)
        plot_and_save(log_scale=True)