import pandas as pd

hug_headers = ["density", "pressure", "temperature", "energy", "sound speed", "entropy",
                                          "shock speed", "part speed", "phase"]


def entropy_increase_table_hugonoit(hugonoit_path):
    hugonoit_df = pd.read_fwf(hugonoit_path, skiprows=1, names=hug_headers)
    pressures = hugonoit_df['pressure']
    entropies = hugonoit_df['entropy']
    dS_list = []
    P_list = []
    for index, pressure in enumerate(pressures):
        entropy = entropies[index]
        if index != 0:
            dS = entropy - entropies[index - 1]
            dP = pressure - pressures[index - 1]
            dS_list.append(dS / (10 ** 3))
            P_list.append(pressure / (10 ** 9))
    return dS_list, P_list

def entropy_increase_analytical_hugonoit(v_p_list, s, T, P_i, C_0, rho_i):
    """
    dS/dP = (s * v_p^2) / (T v_s)
    v_s = C_0 + s v_p
    P = P_i v_s v_p
    :param v_p_list:
    :param s:
    :param T_i:
    :param P_i:
    :param C_0:
    :param rho_i:
    :return:
    """
    dS_list = []
    P_list = []
    for index, v_p in enumerate(v_p_list[1:]):
        dv_p = v_p - v_p_list[index]  # change in particle velocity
        v_s = C_0 + (s * v_p)  # shock velocity
        P = P_i + (rho_i * v_s * v_p)  # pressure
        dS_dv_p = (s * (v_p ** 2)) / (T * v_s)
        dS = ((s * (v_p ** 2)) / (T * v_s)) * dv_p  # change in entropy
        dS_list.append(dS / (10 ** 3))
        P_list.append(P / (10 ** 9))
    return dS_list, P_list
