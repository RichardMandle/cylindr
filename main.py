from config import Config
from trajectory import TrajectoryProcessor
from analysis import CDFAnalysis
from plotting import Plotting
import platform

if __name__ == "__main__":
    print('                 _  _             _        ')
    print('                | |(_)           | |       ')
    print('     ___  _   _ | | _  _ __    __| | _ __  ')
    print('    / __|| | | || || || \'_ \  / _` || \'__| ')
    print('   | (__ | |_| || || || | | || (_| || |    ')
    print('    \___| \__, ||_||_||_| |_| \__,_||_|    ')
    print('           __/ |                           ')
    print('          |___/                        V1.0')
    print("running on", platform.system(), platform.release(), platform.version())
    print("Please cite the paper: 10.1371/journal.pone.0279679")
    config = Config()

    if not config.args.replot:  # Normal mode
        traj_processor = TrajectoryProcessor(config)
        analysis = CDFAnalysis(config, traj_processor)
        analysis.calculate_distances()
        analysis.cylindrical_pcf()
        plotting = Plotting(config, analysis)
        filename = plotting.generate_filename()
        plotting.generate_cdf()
        plt = plotting.plot_cdf()       
        plotting.save_plot_data(f"{filename}_CDF_data.npz") # write cdf + vol + config to npz for retrieva/replot
        print(f'CDF data saved to {filename}_CDF_data.npz')
        plotting.save_config(filename)                      # Save configuration settings
        print(f'Config data written to {filename}.txt')
        
    else:  # Replot mode
        plotting = Plotting(config)
        filename = config.args.output_name
        plotting.load_plot_data(filename)
        plt = plotting.plot_cdf()

    print('Plotting done :)')
    if not config.args.replot and config.args.save:
        print('Saving data...')
        analysis.save_additional_data(f"{filename}_pair_data.npz")