import os
import shutil
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from src.vapor import calc_vapor_mass_fraction
from src.new_and_old_eos import plot

class VMFtimeseries:

    def __init__(self, new_phase_path, old_phase_path, min_time, max_time, txt_fname, output, animation_name="vmf_timeseries.mp4"):
        if os.path.exists(output):
            shutil.rmtree(output)
        if txt_fname in os.listdir(os.getcwd()):
            os.remove(txt_fname)
        self.text = open(txt_fname, 'w')
        self.text.write("new time,old time,new vmf,old vmf\n")
        self.to_path = output
        self.new_phase_path = new_phase_path
        self.old_phase_path = old_phase_path
        self.animation_name = animation_name
        self.min_time = min_time
        self.max_time = max_time
        self.vmf_new_list = []
        self.vmf_old_list = []
        self.new_time_list = []
        self.old_time_list = []

        self.parameter = "entropy"
        self.min_normalize = 0
        self.max_normalize = 8000
        self.square_scale = 4e7

        plt.style.use("dark_background")
        self.cmap = cm.get_cmap('jet')
        self.normalizer = Normalize(self.min_normalize, self.max_normalize)

    def __build_vmf_new_old_plot(self, iteration, new_particles, old_particles, new_time, old_time):
        # Create 2x2 sub plots
        gs = gridspec.GridSpec(2, 2)
        fig = plt.figure(figsize=(10, 10))
        fig.patch.set_facecolor('xkcd:black')
        ax1 = fig.add_subplot(gs[0, 0])  # row 0, col 0
        ax2 = fig.add_subplot(gs[0, 1])  # row 0, col 1
        ax3 = fig.add_subplot(gs[1, :])  # row 1, span all columns
        axs = [ax1, ax2, ax3]
        ax1.set_title("New EoS")
        ax2.set_title("Old EoS")
        ax3.set_xlabel("Time (hrs)")
        ax3.set_ylabel("Silicate Vapor Mass Fraction (%)")
        ax1 = plot(fig=fig, axs=axs, index=0, time=new_time, particles=new_particles, cmap=self.cmap,
                   normalizer=self.normalizer,
                   parameter=self.parameter, square_scale=self.square_scale, flatten=False)
        ax2 = plot(fig=fig, axs=axs, index=1, time=old_time, particles=old_particles, cmap=self.cmap,
                   normalizer=self.normalizer,
                   parameter=self.parameter, square_scale=self.square_scale, flatten=False)
        ax3.plot(
            self.new_time_list,
            self.vmf_new_list,
            color='aqua',
            linewidth=1.0,
            label="New EoS"
        )
        ax3.plot(
            self.old_time_list,
            self.vmf_old_list,
            color='magenta',
            linewidth=1.0,
            label="Old EoS"
        )
        ax3.text(
            self.max_time - (self.max_time * 0.25),
            50,
            "Silicate New EoS VMF: {}%\nSilicate Old EoS VMF: {}%".format(round(self.vmf_new_list[-1], 2), round(self.vmf_old_list[-1], 2)),
            c="white",
            fontsize=10
        )
        ax3.grid(alpha=0.4)
        ax3.set_xlim(self.min_time, self.max_time)
        ax3.set_ylim(0, 60)
        sm = cm.ScalarMappable(norm=self.normalizer, cmap=self.cmap)
        sm.set_array([])
        # cbar = fig.colorbar(sm, ax=axs.flatten()[1])
        cbaxes = inset_axes(ax1, width="30%", height="3%", loc=2, borderpad=1.8)
        cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
        cbar.ax.tick_params(labelsize=6)
        # cbar.ax.xaxis.set_ticks_position('top')
        cbar.ax.set_title(self.parameter.replace("_", " ").title(), fontsize=6)
        legend = ax3.legend(fontsize=6, loc='upper left')
        plt.savefig(self.to_path + "/{}.png".format(iteration), format='png', dpi=200)

    def get_vmfs(self, iteration, new_particles, old_particles, new_time, old_time):
        vmf_new = calc_vapor_mass_fraction(particles=new_particles, phase_path=self.new_phase_path, only_disk=True) * 100.0
        vmf_old = calc_vapor_mass_fraction(particles=old_particles, phase_path=self.old_phase_path, only_disk=True) * 100.0
        self.text.write("{},{},{},{}\n".format(new_time, old_time, vmf_new, vmf_old))
        self.__build_vmf_new_old_plot(
            iteration=iteration,
            new_particles=new_particles,
            old_particles=old_particles,
            new_time=new_time,
            old_time=old_time
        )
