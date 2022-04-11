import os
import pandas as pd

from src.vapor import get_all_particle_vapor_fractions_from_formatted, \
    calc_vapor_mass_fraction_with_circularization_from_formatted

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
angles = ["b073", "b075"]
cutoff_densities = [5, 500, 1000, 2000]

new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"


def get_all_sims(high=True):
    fformat = "{}_{}_{}"
    tformat = "{}{}{}"
    names = []
    titles = []
    high_res_name = None
    high_res_title = None
    for angle in angles:
        for runs in ["new", "old"]:
            n = "n"
            if runs == "old":
                n = "o"
            for cd in cutoff_densities:
                output_name = fformat.format(cd, angle, runs)
                title_name = tformat.format(cd, angle, n)
                titles.append(title_name)
                names.append(output_name)
                if cd == 5 and high and runs == "new":
                    high_res_name = fformat.format(cd, angle, runs) + "_high"
                    high_res_title = tformat.format(cd, angle, n) + "-high"
        if high_res_name is not None and high_res_title is not None:
            names.append(high_res_name)
            titles.append(high_res_title)
        return names, titles


def reformat():
    sims, titles = get_all_sims(high=False)
    for sim, title in zip(sims, titles):
        phase_path = new_phase_path
        if "old" in sim:
            phase_path = old_phase_path
        formatted_path = base_path + "{}/circularized_{}/".format(sim, sim)
        report_path = base_path + "{}/{}_reports/".format(sim, sim)
        for i in os.listdir(report_path):
            df_formatted = pd.read_csv(formatted_path + i)
            df_report = pd.read_csv(report_path + i)
            vmf_uncirc = get_all_particle_vapor_fractions_from_formatted(df=df_formatted, phase_path=phase_path)
            vmf_circ = calc_vapor_mass_fraction_with_circularization_from_formatted(df=df_formatted,
                                                                                    phase_path=phase_path)
            print(vmf_uncirc, vmf_circ)
            del df_report['VMF']
            df_report['DISK_VMF_W_CIRC'] = [vmf_circ]
            df_report['DISK_VMF_WITHOUT_CIRC'] = [vmf_uncirc]
            print(df_report)


reformat()
