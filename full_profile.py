import src.full_suite_profiling as profile

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"

gi_b73_runs = {
    "5_b073_new": {
        "name": "5_b073_new",
        "path": base_path + "formatted_5_b073_new"
    },
    "5_b073_old": {
        "name": "5_b073_old",
        "path": base_path + "formatted_5_b073_old"
    },
    "500_b073_new": {
        "name": "500_b073_new",
        "path": base_path + "formatted_500_b073_new"
    },
    "500_b073_old": {
        "name": "500_b073_old",
        "path": base_path + "formatted_500_b073_old"
    },
    "1000_b073_new": {
        "name": "1000_b073_new",
        "path": base_path + "formatted_1000_b073_new"
    },
    "1000_b073_old": {
        "name": "1000_b073_old",
        "path": base_path + "formatted_1000_b073_old"
    },
    "2000_b073_new": {
        "name": "2000_b073_new",
        "path": base_path + "formatted_2000_b073_new"
    },
    "2000_b073_old": {
        "name": "2000_b073_old",
        "path": base_path + "formatted_2000_b073_old"
    },
}

gi_b75_runs = {
    "5_b075_new": {
        "name": "5_b075_new",
        "path": base_path + "formatted_5_b075_new"
    },
    "5_b075_old": {
        "name": "5_b075_old",
        "path": base_path + "formatted_5_b075_old"
    },
    "500_b075_new": {
        "name": "500_b075_new",
        "path": base_path + "formatted_500_b075_new"
    },
    "500_b075_old": {
        "name": "500_b075_old",
        "path": base_path + "formatted_500_b075_old"
    },
    "1000_b075_new": {
        "name": "1000_b075_new",
        "path": base_path + "formatted_1000_b075_new"
    },
    "1000_b075_old": {
        "name": "1000_b075_old",
        "path": base_path + "formatted_1000_b075_old"
    },
    "2000_b075_new": {
        "name": "2000_b075_new",
        "path": base_path + "formatted_2000_b075_new"
    },
    "2000_b075_old": {
        "name": "2000_b075_old",
        "path": base_path + "formatted_2000_b075_old"
    },
}

profile.build_impact_angle_geometries(gi_b73_runs, start_iteration=0, end_iteration=50, specified_imp_angle=0.73)
