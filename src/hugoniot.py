from sys import exit
import pandas as pd
from scipy.interpolate import interp2d, RectBivariateSpline, LinearNDInterpolator
import warnings
warnings.filterwarnings("ignore")


class Hugoniot:
    def __init__(self):
        pass

    def read_SESAME(self, file_path):
        """
        Reads in the SESAME data file and reformats the data in actually usable Pandas dataframes.
        :param file_path:
        :return:
        """
        self.df = pd.read_csv(file_path, delimiter="\t", header=None)
        self.nt, self.nr = int(self.df[0][0]), int(self.df[1][0])
        grid_length = self.nt * self.nr

        self.temperatures = self.df[0][1:self.nt + 1]
        self.densities = self.df[0][self.nt + 1:self.nt + self.nr + 1] * 1000.0  # to kg/m3
        self.gridded_temperatures = [i for i in self.temperatures for j in self.densities]
        self.gridded_densities = [j for i in self.temperatures for j in self.densities]

        starting_point = self.nt + self.nr + 1
        end_point = starting_point + grid_length
        p_cc = self.df.iloc[starting_point:end_point]
        p_cc['temperature'] = self.gridded_temperatures
        p_cc['density'] = self.gridded_densities
        starting_point += grid_length
        end_point += grid_length
        u_s = self.df.iloc[starting_point:end_point]
        u_s['temperature'] = self.gridded_temperatures
        u_s['density'] = self.gridded_densities

        self.pressures = p_cc.drop(p_cc.columns[[1, 2, 3]], axis=1)
        self.soundspeeds = p_cc.drop(p_cc.columns[[0, 2, 3]], axis=1)
        self.internal_energies = u_s.drop(p_cc.columns[[1, 2, 3]], axis=1)
        self.entropies = u_s.drop(p_cc.columns[[0, 2, 3]], axis=1)

        self.get_interpolation_functions()

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
        return f(x_value, y_value)[0]  # return the interpolated value

    def get_interpolation_functions(self):
        self.interp_funcs = {
            "pressure": self.__get_interp2d_func(self.gridded_temperatures, self.gridded_densities, self.pressures),
            "soundspeed": self.__get_interp2d_func(self.gridded_temperatures, self.gridded_densities, self.soundspeeds),
            "internal_energy": self.__get_interp2d_func(self.gridded_temperatures, self.gridded_densities, self.internal_energies),
            "entropy": self.__get_interp2d_func(self.gridded_temperatures, self.gridded_densities, self.entropies),
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
        :param P1:
        :param P2:
        :param rho1:
        :param rho2:
        :return:
        """
        return U1 + (0.5 * (P1 + P2) * (rho2 - rho1) / (rho1 * rho2))

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
            rho2 = rho1 + 10  # why +10?

            counter = 0
            rho_r = []
            U2_r = []
            while rho2 < 9000:  # why 9000?
                print("At counter: ", counter)
                Us = rho2 / (rho2 - rho1) * Up
                P2 = P1 + rho1 * Up * Us
                if self.interpolate2D("pressure", self.temperatures.tolist()[-1], rho2) > P2:
                    T_prev = self.temperatures.tolist()[0]
                    P_prev = self.interpolate2D("pressure", T_prev, rho2)
                    for T in self.temperatures.tolist()[1:]:
                        P = self.interpolate2D("pressure", T, rho2)
                        if (P - P2) * (P_prev - P2) < 0:
                            T2 = T_prev + ((T - T_prev) / (P - P_prev) * (P2 - P_prev))
                            break
                        T_prev = T
                        P_prev = P
                    U2 = self.interpolate2D("internal_energy", T2, rho1)  # why rho1?
                    U2_calc = self.calculate_internal_energy(U1, P1, P2, rho1, rho2)

                    # these hold old values for some reason?
                    try:
                        rho_r[counter] = rho2
                    except:
                        rho_r.append(rho2)
                    try:
                        U2_r[counter] = U2
                    except:
                        U2_r.append(U2)

                    if counter > 0:
                        if (U2 - U2_calc) * (U2_r[-2] - U2_calc) < 0:
                            rho2 = rho_r[-2] + ((rho2 - rho_r[-2]) / (U2 - U2_r[-2]) * (U2_calc - U2_r[-2]))
                            Us = rho2 / (rho2 - rho1) * Up
                            P2 = P1 + (rho1 * Up * Us)
                            T_prev = self.temperatures.tolist()[0]
                            P_prev = self.interpolate2D("pressure", T_prev, rho2)
                            for T in self.temperatures.tolist()[1:]:
                                P = self.interpolate2D("pressure", T, rho2)
                                if (P - P2) * (P_prev - P2) < 0:
                                    T2 = T_prev + ((T - T_prev) / (P - P_prev) * (P2 - P_prev))
                                    break
                                T_prev = T
                                P_prev = P
                            Rho_s.append(rho2)
                            T_s.append(T2)
                            Us_s.append(Us)
                            Up_s.append(Up)
                            P_s.append(P2)
                            U_s.append(U2)
                            S_s.append(self.interpolate2D("entropy", T2, rho2))
                            # self.update_tracked_list(Rho_s, rho2, i)
                            # self.update_tracked_list(T_s, T2, i)
                            # self.update_tracked_list(Us_s, Us, i)
                            # self.update_tracked_list(Up_s, Up, i)
                            # self.update_tracked_list(P_s, P2, i)
                            # self.update_tracked_list(U_s, U2, i)
                            # S = self.interpolate2D("entropy", T2, rho2)
                            # self.update_tracked_list(S_s, S, i)
                            break

                    counter += 1
                rho2 += 100


        return T_s, Rho_s, Us_s, Up_s, P_s, U_s, S_s



