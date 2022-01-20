#!/usr/bin/env python3
import src.full_suite_profiling as profile

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
base_path_setups = "/home/theia/scotthull/Paper1_SPH/setups/"

gi_b73_runs = {
    "5_b073_new": {
        "name": "5b073n",
        "path": base_path + "5_b073_new/formatted_5_b073_new",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b073", "5_new_gi_setup_b_073.txt"))
    },
    "5_b073_old": {
        "name": "5b073o",
        "path": base_path + "5_b073_old/formatted_5_b073_old",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b073", "5_old_gi_setup_b_073.txt"))
    },
    "500_b073_new": {
        "name": "500b073n",
        "path": base_path + "500_b073_new/formatted_500_b073_new",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b073", "500_new_gi_setup_b_073.txt"))
    },
    "500_b073_old": {
        "name": "500b073o",
        "path": base_path + "500_b073_old/formatted_500_b073_old",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b073", "500_old_gi_setup_b_073.txt"))
    },
    "1000_b073_new": {
        "name": "1000b073n",
        "path": base_path + "1000_b073_new/formatted_1000_b073_new",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b073", "1000_new_gi_setup_b_073.txt"))
    },
    "1000_b073_old": {
        "name": "1000b073o",
        "path": base_path + "1000_b073_old/formatted_1000_b073_old",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b073", "1000_old_gi_setup_b_073.txt"))
    },
    "2000_b073_new": {
        "name": "2000b073n",
        "path": base_path + "2000_b073_new/formatted_2000_b073_new",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b073", "2000_new_gi_setup_b_073.txt"))
    },
    "2000_b073_old": {
        "name": "2000b073o",
        "path": base_path + "2000_b073_old/formatted_2000_b073_old",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b073", "2000_old_gi_setup_b_073.txt"))
    },
}

gi_b75_runs = {
    "5_b075_new": {
        "name": "5b075n",
        "path": base_path + "5_b075_new/formatted_5_b075_new",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b075", "5_new_gi_setup_b_075.txt"))
    },
    "5_b075_old": {
        "name": "5b075o",
        "path": base_path + "5_b075_old/formatted_5_b075_old",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b075", "5_old_gi_setup_b_075.txt"))
    },
    "500_b075_new": {
        "name": "500b075n",
        "path": base_path + "500_b075_new/formatted_500_b075_new",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b075", "500_new_gi_setup_b_075.txt"))
    },
    "500_b075_old": {
        "name": "500b075o",
        "path": base_path + "500_b075_old/formatted_500_b075_old",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b075", "500_old_gi_setup_b_075.txt"))
    },
    "1000_b075_new": {
        "name": "1000b075n",
        "path": base_path + "1000_b075_new/formatted_1000_b075_new",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b075", "1000_new_gi_setup_b_075.txt"))
    },
    "1000_b075_old": {
        "name": "1000b075o",
        "path": base_path + "1000_b075_old/formatted_1000_b075_old",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b075", "1000_old_gi_setup_b_075.txt"))
    },
    "2000_b075_new": {
        "name": "2000b075n",
        "path": base_path + "2000_b075_new/formatted_2000_b075_new",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b075", "2000_new_gi_setup_b_075.txt"))
    },
    "2000_b075_old": {
        "name": "2000b075o",
        "path": base_path + "2000_b075_old/formatted_2000_b075_old",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b075", "2000_old_gi_setup_b_075.txt"))
    },
}

profile.build_impact_angle_geometries("gi_b73_runs", gi_b73_runs, start_iteration=0, end_iteration=30, specified_imp_angle=0.73)
profile.build_impact_velocity_charts("gi_b73_runs", gi_b73_runs, start_iteration=0, end_iteration=30)
profile.build_vmf_timeplots("gi_b73_runs", gi_b73_runs, start_iteration=0, end_iteration=3000, increment=100)
profile.map_disk_to_phase_profile("gi_b73_runs", gi_b73_runs, end_iteration=3000)
profile.map_disk_to_phase_profile_eos_charts("gi_b73_runs", gi_b73_runs, end_iteration=3000)
profile.get_end_profile_reports(gi_b73_runs, end_iteration=3000, number_processes=200)
profile.disk_temperature_vs_radius(gi_b73_runs, gi_b73_runs, iteration=3000)
profile.build_scenes("gi_b73_runs", gi_b73_runs, start_iteration=0, end_iteration=3000, increment=1, to_path="gi_b73_runs_scenes")

profile.build_impact_angle_geometries("gi_b75_runs", gi_b75_runs, start_iteration=0, end_iteration=30, specified_imp_angle=0.75)
profile.build_impact_velocity_charts("gi_b75_runs", gi_b75_runs, start_iteration=0, end_iteration=30)
profile.build_vmf_timeplots("gi_b75_runs", gi_b75_runs, start_iteration=0, end_iteration=3000, increment=100)
profile.map_disk_to_phase_profile("gi_b75_runs", gi_b75_runs, end_iteration=3000)
profile.map_disk_to_phase_profile_eos_charts("gi_b75_runs", gi_b75_runs, end_iteration=3000)
profile.get_end_profile_reports(gi_b75_runs, end_iteration=3000, number_processes=200)
profile.disk_temperature_vs_radius(gi_b75_runs, gi_b75_runs, iteration=3000)
profile.build_scenes("gi_b75_runs", gi_b75_runs, start_iteration=0, end_iteration=3000, increment=1, to_path="gi_b75_runs_scenes")
