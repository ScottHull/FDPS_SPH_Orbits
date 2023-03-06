import ast
import pandas as pd
import matplotlib.pyplot as plt

runs = [
    "5b073S", "5b073S-high", "500b073S", "1000b073S", "2000b073S", 
    "5b073N", "500b073N", "1000b073N", "2000b073N",
    "5b075S", "500b075S", "1000b075S", "2000b075S", 
    "5b075N", "500b075N", "1000b075N", "2000b075N", "2000b075N-low"
]

base_path = "/home/theia/scotthull/Paper1_SPH/gi"
hugoniot_results_path = "figures-selected"
end_iteration = 1800



for run in runs:
    angle = "b073"
    if "b075" in run:
        angle = "b075"
    hugoniot_file = f"{hugoniot_results_path}/{angle}_maxvals_at_endstate/{run}_max_vals.csv"
    circ_end_file = f"{base_path}/{run}/circularized_{run}2/{end_iteration}.csv"
    # read in the hugoniot file as a dictionary
    hugoniot_dict = {}
    # reading the data from the file
    with open(hugoniot_file) as f:
        hugoniot_dict = ast.literal_eval(f.read())
    f.close()
    # read data from the circularized file
    df = pd.read_csv(circ_end_file)
    disk_df = df[df["id"].isin(list(hugoniot_dict.keys()))]

