from sys import exit
import pandas as pd
from scipy.interpolate import interp2d, RectBivariateSpline, LinearNDInterpolator
import warnings
warnings.filterwarnings("ignore")


class Hugoniot_Granite:

    def __init__(self):
        pass

    def read_SESAME(self, file_path):
        """
        Reads in the SESAME data file and reformats the data in actually usable Pandas dataframes.
        :param file_path:
        :return:
        """
        df = pd.read_csv(file_path, delimiter="\t", header=None)
        self.nt, self.nr = int(df[0][0]), int(df[1][0])
        grid_length = self.nt * self.nr

        self.temperatures = df[0][1:self.nt + 1]
        self.densities = df[0][self.nt + 1:self.nt + self.nr + 1] * 1000.0  # to kg/m3
        self.gridded_temperatures = [i for i in self.temperatures for j in self.densities]
        self.gridded_densities = [j for i in self.temperatures for j in self.densities]

        starting_point = self.nt + self.nr + 1
        end_point = starting_point + grid_length
        p_cc = df.iloc[starting_point:end_point]
        p_cc['temperature'] = self.gridded_temperatures
        p_cc['density'] = self.gridded_densities
        starting_point += grid_length
        end_point += grid_length
        u_s = df.iloc[starting_point:end_point]
        u_s['temperature'] = self.gridded_temperatures
        u_s['density'] = self.gridded_densities

        self.pressures = p_cc.drop(p_cc.columns[[1, 2, 3]], axis=1)
        self.soundspeeds = p_cc.drop(p_cc.columns[[0, 2, 3]], axis=1)
        self.internal_energies = u_s.drop(p_cc.columns[[1, 2, 3]], axis=1)
        self.entropies = u_s.drop(p_cc.columns[[0, 2, 3]], axis=1)

        self.get_interpolation_functions(self.gridded_temperatures, self.gridded_densities)

    def __get_interp2d_func(self, x_array, y_array, z_array):
        return LinearNDInterpolator((x_array, y_array), z_array[list(z_array.keys())[0]])

    def interpolate2D(self, variable: str, x_value, y_value):
        """
        Interpolates the Z array given a x and y array, and the x and y locations at which z must be interpolated.
        :param x_array:
        :param y_array:
        :param z_array:
        :param x_value:
        :param y_value:
        :return:
        """
        f = self.interp_funcs[variable]
        return f(x_value, y_value)  # return the interpolated value

    def get_interpolation_functions(self, x_array, y_array):
        self.interp_funcs = {
            "pressure": self.__get_interp2d_func(x_array, y_array, self.pressures),
            "soundspeed": self.__get_interp2d_func(x_array, y_array, self.soundspeeds),
            "internal_energy": self.__get_interp2d_func(x_array, y_array, self.internal_energies),
            "entropy": self.__get_interp2d_func(x_array, y_array, self.entropies),
        }


    def initial_conditions(self, aneos_hugoniot_filepath):
        """
        Gives the initial conditions for the calculated Hugoniot given the ANEOS Hugoniot file.
        :param aneos_hugoniot_filepath:
        :return:
        """
        # get ANEOS Hugoniot data and store it in order to get initial conditions
        df = pd.read_fwf(aneos_hugoniot_filepath, header=None)
        self.rho_h, self.T_h, self.P_h, self.U_h, self.S_h, self.Us_h, self.Up_h = df[0], df[1], df[2] * 10 ** 9, \
                                                                               df[4], df[5], df[6] * 1000, df[7] * 1000
        rho1, T1 = self.rho_h[0], self.T_h[0]
        P1 = self.interpolate2D("pressure", T1, rho1)
        C1 = self.interpolate2D("soundspeed", T1, rho1)
        U1 = self.interpolate2D("internal_energy", T1, rho1)
        S1 = self.interpolate2D("entropy", T1, rho1)
        return rho1, P1, C1, U1, S1

    def calculate_internal_energy(self, U1, P1, P2, rho1, rho2):
        """
        A function to calculate internal energy, U.
        Derived from the conservation of energy Rankine-Hugoniot condition:
        h1 + 0.5u1^2 = h2 + 0.5u2^2
        where u is is the internal energy and h is the specific enthalpy (i.e. per unit mass).
        Using both conservation of mass and conservation of momentum, the Hugoniot equation for conservation of energy
        becomes (because u1 and u2 are eliminated in the process):
        h2 - h1 = 0.5 * (1/rho2 + 1/rho1) * (P2 - P1)
        where u2 = h2 and u1 = h1.
        Here, the densities reduce to a specific volume element V = 1/rho to form the equation:
        h2 - h1 = 0.5 * (P2 + P1) * (V2 - V1)
        :param U1:
        :param P1:
        :param P2:
        :param rho1:
        :param rho2:
        :return:
        """
        return U1 + (0.5 * (P1 + P2) * (rho2 - rho1) / (rho1 * rho2))

    def calculate_T2(self, P2, rho2):
        """
        Linearly interpolates T2 from the 2D-interpolated pressure.
        :param P2:
        :param rho2:
        :return:
        """
        T_prev = self.temperatures.tolist()[0]  # the initial temperature for the interpolation
        P_prev = self.interpolate2D("pressure", T_prev, rho2)  # pressure corresponding to the initial temperature
        for T in self.temperatures.tolist()[1:]:
            P = self.interpolate2D("pressure", T, rho2)  # the pressure at T and rho2
            # begin temperature interpolation
            if (P - P2) * (P_prev - P2) < 0:  # this is a midpoint condition between P_prev and P2
                # assume that temperature scales with pressure according to (P2_p_prev) / (P - P_prev)
                return T_prev + ((T - T_prev) / (P - P_prev) * (P2 - P_prev))  # this is T2
            T_prev = T
            P_prev = P

    def update_tracked_list(self, list_to_update, new_value, index):
        """
        Tries to update the list with the new value at the index location. If the index is out of range, it will add the new value to the end of the list.
        :param list_to_update:
        :param new_value:
        :return:
        """
        if index < len(list_to_update) - 1:
            list_to_update[index] = new_value
        else:
            list_to_update.append(new_value)
        return list_to_update

    def rankine_hugoniot_equations(self, rho1, P1, U1):
        T_s = []
        Rho_s = []
        Us_s = []
        Up_s = []
        P_s = []
        U_s = []
        S_s = []

        for i in range(len(self.Up_h)):  # 50 is the length of the given Hugoniot data
            print("At loop number: ", i)
            Up = self.Up_h[i]  # why this index location?  This is a loop in the MATLAB code for no reason.
            rho2 = rho1 + 10  # bump rho2 from rho1 by 10

            counter = 0  # this is an index location of the previously unsuccessful rho2 convergence loop
            tracked_rho2 = []  # tracks previous rho2 values
            tracked_U2 = []  # tracks previous U2 values

            while rho2 < 9000:  # sets a limit for rho2 convergence
                print("At counter: ", counter)
                Us = rho2 / (rho2 - rho1) * Up  # conservation of momentum
                # P2 - P1 = rho1 (Us - u1) * (Up - U1)
                # if U1 = 0, then P2 = P1 * rho1 * Up * Us
                # U1 = 0 => material initially at rest?
                P2 = P1 + rho1 * Up * Us  # the pressure at the second state
                if self.interpolate2D("pressure", self.temperatures.tolist()[-1], rho2) > P2:  # enforce an increase in pressure from the previous iteration
                    T2 = self.calculate_T2(P2, rho2)  # the temperature at the second state
                    U2 = self.interpolate2D("internal_energy", T2, rho1)  # why rho1?
                    U2_calc = self.calculate_internal_energy(U1, P1, P2, rho1, rho2)

                    # these hold old values for some reason?
                    try:
                        tracked_rho2[counter] = rho2
                    except:
                        tracked_rho2.append(rho2)
                    try:
                        tracked_U2[counter] = U2
                    except:
                        tracked_U2.append(U2)

                    if counter > 0:  # if not the initial loop (i.e. we have U2 and rho2 values to compare to)
                        if (U2 - U2_calc) * (tracked_U2[-2] - U2_calc) < 0:
                            rho2 = tracked_rho2[-2] + ((rho2 - tracked_rho2[-2]) / (U2 - tracked_U2[-2]) * (U2_calc - tracked_U2[-2]))
                            Us = rho2 / (rho2 - rho1) * Up
                            P2 = P1 + (rho1 * Up * Us)
                            T2 = self.calculate_T2(P2, rho2)
                            
                            # update lists with new values
                            Rho_s.append(rho2)
                            T_s.append(T2)
                            Us_s.append(Us)
                            Up_s.append(Up)
                            P_s.append(P2)
                            U_s.append(U2)
                            S_s.append(self.interpolate2D("entropy", T2, rho2))
                            break  # exit while loop

                    counter += 1  # we did not find a rho2 or P2, must bump rho2 and go again
                rho2 += 100  # bump rho2 by 100


        return T_s, Rho_s, Us_s, Up_s, P_s, U_s, S_s



class Hugoniot_ANEOS:

    def __init__(self):
        pass


    def read_ANEOS(self, file_path):
        df = pd.read_fwf(file_path, sep="\t", header=None, skiprows=2)
        self.gridded_densities = df[0]
        self.gridded_internal_energies = df[1]
        self.temperatures = df[2]
        self.pressures = df[3]
        self.soundspeeds = df[4]
        self.entropies = df[5]

        self.set_temperatures = list(set(self.temperatures))

        self.densities = list(set(self.gridded_densities))
        self.internal_energies = list(set(self.gridded_internal_energies))

        self.get_interpolation_functions(self.gridded_densities.tolist(), self.gridded_internal_energies.tolist())

    def __get_interp2d_func(self, x_array, y_array, z_array):
        return LinearNDInterpolator((x_array, y_array), z_array)

    def interpolate2D(self, variable: str, x_value, y_value):
        """
        Interpolates the Z array given a x and y array, and the x and y locations at which z must be interpolated.
        :param x_array:
        :param y_array:
        :param z_array:
        :param x_value:
        :param y_value:
        :return:
        """
        f = self.interp_funcs[variable]
        return f(x_value, y_value)  # return the interpolated value

    def get_interpolation_functions(self, x_array, y_array):
        self.interp_funcs = {
            "pressure": self.__get_interp2d_func(x_array, y_array, self.pressures),
            "soundspeed": self.__get_interp2d_func(x_array, y_array, self.soundspeeds),
            "internal_energy": self.__get_interp2d_func(x_array, y_array, self.gridded_internal_energies),
            "entropy": self.__get_interp2d_func(x_array, y_array, self.entropies),
            "temperature": self.__get_interp2d_func(x_array, y_array, self.temperatures)
        }

    def initial_conditions(self, aneos_hugoniot_filepath):
        """
        Gives the initial conditions for the calculated Hugoniot given the ANEOS Hugoniot file.
        :param aneos_hugoniot_filepath:
        :return:
        """
        # get ANEOS Hugoniot data and store it in order to get initial conditions
        df = pd.read_fwf(aneos_hugoniot_filepath, header=None, skiprows=1)
        self.rho_h, self.P_h, self.T_h, self.U_h, self.C_h, self.S_h, self.Us_h, self.Up_h = [df[i] for i in list(df.keys())[:-1]]
        self.P_h /= 10 ** 9
        self.U_h *= 1000
        self.Us_h *= 1000
        self.Up_h *= 1000
        rho1, U1 = self.rho_h[0], self.U_h[0]
        P1 = self.interpolate2D("pressure", rho1, U1)
        C1 = self.interpolate2D("soundspeed", rho1, U1)
        T1 = self.interpolate2D("temperature", rho1, U1)
        S1 = self.interpolate2D("entropy", rho1, U1)
        return rho1, P1, C1, U1, S1

    def calculate_internal_energy(self, U1, P1, P2, rho1, rho2):
        """
        A function to calculate internal energy, U.
        Derived from the conservation of energy Rankine-Hugoniot condition:
        h1 + 0.5u1^2 = h2 + 0.5u2^2
        where u is is the internal energy and h is the specific enthalpy (i.e. per unit mass).
        Using both conservation of mass and conservation of momentum, the Hugoniot equation for conservation of energy
        becomes (because u1 and u2 are eliminated in the process):
        h2 - h1 = 0.5 * (1/rho2 + 1/rho1) * (P2 - P1)
        where u2 = h2 and u1 = h1.
        Here, the densities reduce to a specific volume element V = 1/rho to form the equation:
        h2 - h1 = 0.5 * (P2 + P1) * (V2 - V1)
        :param U1:
        :param P1:
        :param P2:
        :param rho1:
        :param rho2:
        :return:
        """
        return U1 + (0.5 * (P1 + P2) * (rho2 - rho1) / (rho1 * rho2))

    def calculate_T2(self, P2, rho2):
        """
        Linearly interpolates T2 from the 2D-interpolated pressure.
        :param P2:
        :param rho2:
        :return:
        """
        T_prev = self.temperatures.tolist()[0]  # the initial temperature for the interpolation
        P_prev = self.interpolate2D("pressure", T_prev, rho2)  # pressure corresponding to the initial temperature
        for T in self.temperatures.tolist()[1:]:
            P = self.interpolate2D("pressure", T, rho2)  # the pressure at T and rho2
            # begin temperature interpolation
            if (P - P2) * (P_prev - P2) < 0:  # this is a midpoint condition between P_prev and P2
                # assume that temperature scales with pressure according to (P2_p_prev) / (P - P_prev)
                return T_prev + ((T - T_prev) / (P - P_prev) * (P2 - P_prev))  # this is T2
            T_prev = T
            P_prev = P

    def update_tracked_list(self, list_to_update, new_value, index):
        """
        Tries to update the list with the new value at the index location. If the index is out of range, it will add the new value to the end of the list.
        :param list_to_update:
        :param new_value:
        :return:
        """
        if index < len(list_to_update) - 1:
            list_to_update[index] = new_value
        else:
            list_to_update.append(new_value)
        return list_to_update

    def rankine_hugoniot_equations(self, rho1, P1, U1):
        T_s = []
        Rho_s = []
        Us_s = []
        Up_s = []
        P_s = []
        U_s = []
        S_s = []

        for i in range(len(self.Up_h)):  # 50 is the length of the given Hugoniot data
            print("At loop number: ", i)
            Up = self.Up_h[i]  # why this index location?  This is a loop in the MATLAB code for no reason.
            rho2 = rho1 + 10  # bump rho2 from rho1 by 10

            counter = 0  # this is an index location of the previously unsuccessful rho2 convergence loop
            tracked_rho2 = []  # tracks previous rho2 values
            tracked_U2 = []  # tracks previous U2 values

            while rho2 < 9000:  # sets a limit for rho2 convergence
                print("At counter: ", counter)
                Us = rho2 / (rho2 - rho1) * Up  # conservation of momentum
                # P2 - P1 = rho1 (Us - u1) * (Up - U1)
                # if U1 = 0, then P2 = P1 * rho1 * Up * Us
                # U1 = 0 => material initially at rest?
                P2 = P1 + rho1 * Up * Us  # the pressure at the second state
                if self.interpolate2D("pressure", self.temperatures.tolist()[-1],
                                      rho2) > P2:  # enforce an increase in pressure from the previous iteration
                    T2 = self.calculate_T2(P2, rho2)  # the temperature at the second state
                    U2 = self.interpolate2D("internal_energy", rho2)  # why rho1?
                    U2_calc = self.calculate_internal_energy(U1, P1, P2, rho1, rho2)

                    # these hold old values for some reason?
                    try:
                        tracked_rho2[counter] = rho2
                    except:
                        tracked_rho2.append(rho2)
                    try:
                        tracked_U2[counter] = U2
                    except:
                        tracked_U2.append(U2)

                    if counter > 0:  # if not the initial loop (i.e. we have U2 and rho2 values to compare to)
                        if (U2 - U2_calc) * (tracked_U2[-2] - U2_calc) < 0:
                            rho2 = tracked_rho2[-2] + (
                                        (rho2 - tracked_rho2[-2]) / (U2 - tracked_U2[-2]) * (U2_calc - tracked_U2[-2]))
                            Us = rho2 / (rho2 - rho1) * Up
                            P2 = P1 + (rho1 * Up * Us)
                            T2 = self.calculate_T2(P2, rho2)

                            # update lists with new values
                            Rho_s.append(rho2)
                            T_s.append(T2)
                            Us_s.append(Us)
                            Up_s.append(Up)
                            P_s.append(P2)
                            U_s.append(U2)
                            S_s.append(self.interpolate2D("entropy", T2, rho2))
                            break  # exit while loop

                    counter += 1  # we did not find a rho2 or P2, must bump rho2 and go again
                rho2 += 100  # bump rho2 by 100

        return T_s, Rho_s, Us_s, Up_s, P_s, U_s, S_s