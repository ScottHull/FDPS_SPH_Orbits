import os
import shutil
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from src.identify import ParticleMap
from src.combine import CombineFile
from src.new_and_old_eos import get_parameter_from_particles


class DiskProperties:

    def __init__(self, new_eos_path, old_eos_path, properties=None, formatted=True):
        if properties is None:
            properties = ['entropy', 'pressure', 'density', 'temperature']
        self.properties = properties
        self.formatted = formatted
        self.new_eos_path = new_eos_path
        self.old_eos_path = old_eos_path
        self.end_state_particles_new_eos = {}
        self.end_state_particles_old_eos = {}
        self.earth_radius = 6371 * 1000
        self.labels = {
            1: "Target Silicate",
            2: "Target Iron",
            3: "Impactor Silicate",
            4: "Impactor Iron"
        }

    def get_end_state_disk_particles(self, end_iteration, formatted=False, number_processes=200):
        if formatted:
            pm_new = ParticleMap(path=self.new_eos_path + "/{}.csv".format(end_iteration), center=True,
                                 relative_velocity=False, formatted=self.formatted).collect_particles()
            pm_old = ParticleMap(path=self.old_eos_path + "/{}.csv".format(end_iteration), center=True,
                                 relative_velocity=False, formatted=self.formatted).collect_particles()
        else:
            cf = CombineFile(num_processes=number_processes, time=end_iteration, output_path=self.new_eos_path)
            formatted_time = cf.sim_time
            combined_file = cf.combine()
            f = os.getcwd() + "/merged_{}.dat".format(end_iteration)
            pm_new = ParticleMap(path=f, center=True, relative_velocity=False).collect_particles()
            os.remove(f)
            
            cf = CombineFile(num_processes=number_processes, time=end_iteration, output_path=self.old_eos_path)
            formatted_time = cf.sim_time
            combined_file = cf.combine()
            f = os.getcwd() + "/merged_{}.dat".format(end_iteration)
            pm_old = ParticleMap(path=f, center=True, relative_velocity=False).collect_particles()
            os.remove(f)
        for i in list(zip(["new", "old"], [pm_new, pm_old])):
            ids = [p.particle_id for p in i[1] if p.label == "DISK"]
            label = [p.label for p in i[1] if p.label == "DISK"]
            if i[0] == "new":
                self.end_state_particles_new_eos = dict(zip(ids, label))
            else:
                self.end_state_particles_old_eos = dict(zip(ids, label))
        return pm_new, pm_old

    def __scatter2D(self, ax, particles, y_property, s=2):
        for tag in set([p.tag for p in particles]):
            ax.scatter(
                [p.distance / self.earth_radius for p in particles if p.tag == tag],
                [get_parameter_from_particles(particle=p, parameter=y_property) for p in particles if p.tag == tag],
                label=self.labels[tag],
                s=s
            )

    def profile_disk_at_time(self, new_eos_particles, old_eos_particles, fig, axs, iteration, name="disk_profile_{}.png"):
        axs = axs.flatten()
        new_disk_particles = [p for p in new_eos_particles if self.end_state_particles_new_eos[p.particle_id] == "DISK"]
        old_disk_particles = [p for p in old_eos_particles if self.end_state_particles_old_eos[p.particle_id] == "DISK"]
        index_tracker = 0  # for incrementing a duel-column new/old EOS setup
        for index, prop in enumerate(self.properties):
            new_index, old_index = index_tracker, index_tracker + 1
            new_ax, old_ax = axs[new_index], axs[old_index]
            self.__scatter2D(ax=new_ax, particles=new_disk_particles, y_property=prop)
            self.__scatter2D(ax=old_ax, particles=old_disk_particles, y_property=prop)
            y_dat_new_and_old = [get_parameter_from_particles(particle=p, parameter=prop) for p in
                                 new_disk_particles] + [get_parameter_from_particles(particle=p, parameter=prop) for p
                                                        in old_disk_particles]
            min_y, max_y = min(y_dat_new_and_old), max(y_dat_new_and_old)
            new_ax.set_ylabel(min_y, max_y)
            for ax in [new_ax, old_ax]:
                ax.grid(alpha=0.4)
                if index - -1 == len(self.properties):
                    ax.set_xlabel("r'Radius ($R_{\bigoplus}$)',")
            index_tracker += 2
        fig.save(name.format(iteration), format='png')
        return fig
