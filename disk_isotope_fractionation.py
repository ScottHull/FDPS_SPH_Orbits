
mve_isotopes = {
    r"$\delta^{41}$K": {
        "41K": {'mass': 40.961},
        "39K": {'mass': 38.963},
    },
    r"$\delta^{87}$Rb": {
        "87Rb": {'mass': 86.909},
        "85Rb": {'mass': 84.911},
    },
    r"$\delta^{66}$Zn": {
        "66Zn": {'mass': 65.926},
        "64Zn": {'mass': 63.929},
    },
    r"$\delta^{65}$Cu": {
        "65Cu": {'mass': 64.927},
        "63Cu": {'mass': 62.929},
    },
    r"$\delta^{71}$Ga": {
        "71Ga": {'mass': 70.924},
        "69Ga": {'mass': 68.925},
    },
    r"$\delta^{124}$Sn": {
        "124Sn": {'mass': 123.905},
        "116Sn": {'mass': 115.901},
    },
}

def pressure_scaling_law(temeprature):
    return -51.17 - (1048.8 / temeprature) + (14.60 * temeprature)


