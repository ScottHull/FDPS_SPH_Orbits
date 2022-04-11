import os
import pandas as pd


base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
angle = ["b073", "b075"]
cutoff_densities = [5, 500, 1000, 2000]


def get_all_sims(high=True):
    fformat = "{}_{}_{}"
    tformat = "{}{}{}"
    names = []
    titles = []
    high_res_name = None
    high_res_title = None
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
        formatted_path = base_path + "{}/circularized_{}/".format(sim, sim)
        report_path = base_path + "{}/{}_reports/".format(sim, sim)
        df_formatted = pd.read_csv(formatted_path)
        df_report = pd.read_csv(report_path)
        vmf_
