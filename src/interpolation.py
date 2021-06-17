import pandas as pd
from operator import itemgetter


def __get_neighbors(val, array):
    # returns array index value of neighbors
    if val < array[0]:
        return (1, 0)
    elif val > array[-1]:
        return (len(array) - 1, len(array) - 2)
    for index, i in enumerate(array):
        if index < len(array) + 1:
            if i < val <= array[index + 1]:
                return index + 1, index


def interpolate1d(val, val_array, interp_array):
    sorted_interp_array = [x for _, x in sorted(zip(val_array, interp_array))]
    sorted_val_array = list(reversed(val_array))

    neighbor_indices = __get_neighbors(val=val, array=sorted_val_array)
    val_neighbors = (sorted_val_array[neighbor_indices[0]], sorted_val_array[neighbor_indices[1]])
    interp_neighbors = (sorted_interp_array[neighbor_indices[0]], sorted_interp_array[neighbor_indices[1]])

    return (((val_neighbors[1] - val) / (val_neighbors[1] - val_neighbors[0])) * interp_neighbors[0]) + \
           (((val_neighbors[0] - val) / (val_neighbors[1] - val_neighbors[0])) * interp_neighbors[1])


class GenericTrilinearInterpolation:

    def __init__(self, var1_array, var2_array, var3_array, var1, var2, grid_length=120):
        self.var1 = var1
        self.var2 = var2
        self.var1_array = var1_array
        self.var2_array = var2_array
        self.var3_array = var3_array
        self.grid_length = grid_length

        self.var1_neighbors = None
        self.var2_neighbors = None
        self.var3_neighbors = None

        self.p1 = None
        self.p2 = None
        self.s11 = None
        self.s12 = None
        self.s21 = None
        self.s22 = None
        self.u11 = None
        self.u12 = None
        self.u21 = None
        self.u22 = None
        self.u1 = None
        self.u2 = None

    def restrict_var1_indices_to_single_var1(self, var1_array, given_var1, bound='lower'):
        """
        Gets the bounds on 3 values for the purpose of interpolating var2 when given the ordered var1 array.
        :param var1_array:
        :param given_var1:
        :return:
        """
        b1 = None
        b2 = None

        if not bound == 'lower':
            for index, i in enumerate(var1_array):
                if i > given_var1:
                    b1 = index
                    b2 = index + self.grid_length
                    if b2 > len(var1_array) - 1:
                        b2 = len(var1_array) - 1
                    break
        else:
            for index, i in enumerate(reversed(var1_array)):
                if i < given_var1:
                    # b2 = len(var1_array) - (index - 1) - self.grid_length - 1
                    # b1 = len(var1_array) - (index) - (2 * self.grid_length)
                    b1 = (len(var1_array) - 1) - index - (self.grid_length - 1)
                    b2 = b1 + self.grid_length
                    if b1 < 0:
                        b1 = 0
                    if b2 == 0:
                        b2 = self.grid_length
                    break
        if b1 is None and b2 is None:
            b1 = 0
            b2 = self.grid_length
        return b1, b2

    def calc_distance(self, given, sample):
        """
        A simple function for calculating directional distances between a given value and a sample value.
        :param given:
        :param sample:
        :return:
        """
        distance = given - sample
        return distance

    def get_var2_neighbors(self, restriced_indices):
        """
        Get the nearest var2 neighbors to a given var2 value in a restricted array.
        :param d1_indices:
        :param d2_indices:
        :param self.var2_array:
        :param self.var2:
        :return:
        """

        var2_array_restricted = self.var2_array[restriced_indices[0]:restriced_indices[1]]

        min_distance = None
        min_distance_index = None

        for index, i in enumerate(var2_array_restricted):
            distance = self.calc_distance(given=self.var2, sample=i)
            if min_distance is None:
                min_distance = distance
                min_distance_index = index
            elif abs(distance) < abs(min_distance):
                min_distance = distance
                min_distance_index = index
        if min_distance < 0:
            if min_distance_index <= 0:
                return (min_distance_index, min_distance_index + 1)
            elif min_distance_index + 1 == self.grid_length:
                return (min_distance_index - 2, min_distance_index - 1)
            return (min_distance_index - 1, min_distance_index)
        elif min_distance > 0:
            if min_distance_index + 1 == self.grid_length:
                return (min_distance_index - 2, min_distance_index - 1)
            return (min_distance_index, min_distance_index + 1)
        else:
            if min_distance_index >= self.grid_length - 1:
                return (min_distance_index - 2, min_distance_index - 1)
            elif min_distance_index < 0:
                return (min_distance_index + 1, min_distance_index + 2)
            return (min_distance_index, min_distance_index + 1)

    def get_var3_neighbor_values(self, s11, s12, s21, s22, lower_var1_restricted_indices,
                                 upper_var1_restricted_indices):
        """
        Get the 4 nearest var3 values
        s11: the index position of the lower var2 neighbor at var1 d1
        s12: the index position of the upper var2 neighbor at var1 d1
        s21: the index position of the lower var2 neighbor at var1 d2
        s2: the index position of the upper var2 neighbor at var1 d2
        :param s11:
        :param s12:
        :param s21:
        :param s22:
        :param self.var3_array:
        :return:
        """

        e11 = self.var3_array[lower_var1_restricted_indices[0]:lower_var1_restricted_indices[1]][s11]
        e12 = self.var3_array[lower_var1_restricted_indices[0]:lower_var1_restricted_indices[1]][s12]
        e21 = self.var3_array[upper_var1_restricted_indices[0]:upper_var1_restricted_indices[1]][s21]
        e22 = self.var3_array[upper_var1_restricted_indices[0]:upper_var1_restricted_indices[1]][s22]

        return (e11, e12, e21, e22)

    def bilinear_interpolate(self, x1, x2, x, y1, y2, y, q11, q12, q21, q22):
        f = (q11 * (x2 - self.var1) * (y2 - self.var2) +
             q21 * (self.var1 - x1) * (y2 - self.var2) +
             q12 * (x2 - self.var1) * (self.var2 - y1) +
             q22 * (self.var1 - x1) * (self.var2 - y1)
             ) / ((x2 - x1) * (y2 - y1) + 0.0)
        return f

    def linear_interpolate(self, x1, x2, x, q1, q2):
        f = (((x2 - x) / (x2 - x1)) * q1) + (((x - x1) / (x2 - x1)) * q2)
        return f

    def restrict(self):

        # now, given that we'll have var3 values within a range of a single var1 in df, we must restrict the var1 array
        # the following 2 functions will return the 'upper' and 'lower' nearest neighbor index ranges to given_var1
        # we can use these to restrict the arrays to within these index ranges
        # d1 indices will give the index range for var1 which gives the 'lower' nearest neighbor
        # d2 indices will give the index range for var1 which gives the 'upper' nearest neighbor
        d1_indices = self.restrict_var1_indices_to_single_var1(var1_array=self.var1_array, given_var1=self.var1,
                                                               bound='lower')
        d2_indices = self.restrict_var1_indices_to_single_var1(var1_array=self.var1_array, given_var1=self.var1,
                                                               bound='upper')

        # now, restrict the var1 array based on d1_indices and d2_indices
        var1_1_array = self.var1_array[d1_indices[0]:d1_indices[1]]
        var1_2_array = self.var1_array[d2_indices[0]:d2_indices[1]]

        # we will restrict the var2 array also based on the index ranges given by d1_indices and d2_indices
        # the following 2 functions will return the nearest var2 neighbors to given_var2 with the restricted upper and lower array
        lower_var2_neighbors = self.get_var2_neighbors(restriced_indices=d1_indices)
        upper_var2_neighbors = self.get_var2_neighbors(restriced_indices=d2_indices)

        # because we need nearest var3 neighbors for interpolation, we get the var3 values at the same index location as the var2 values
        var3_neighbors = self.get_var3_neighbor_values(s11=lower_var2_neighbors[0], s12=lower_var2_neighbors[1],
                                                       s21=upper_var2_neighbors[0], s22=upper_var2_neighbors[1],
                                                       lower_var1_restricted_indices=d1_indices,
                                                       upper_var1_restricted_indices=d2_indices)

        # package the neighbor values up for use for interpolation
        var1_neighbor_values = (
            self.var1_array[d1_indices[0]:d1_indices[1]][0], self.var1_array[d2_indices[0]:d2_indices[1]][0])
        restricted_var2_array_lower = self.var2_array[d1_indices[0]:d1_indices[1]]
        restricted_var2_array_upper = self.var2_array[d2_indices[0]:d2_indices[1]]
        var2_neighbor_values = (restricted_var2_array_lower[lower_var2_neighbors[0]],
                                restricted_var2_array_lower[lower_var2_neighbors[1]],
                                restricted_var2_array_upper[upper_var2_neighbors[0]],
                                restricted_var2_array_upper[upper_var2_neighbors[1]])

        return (var1_neighbor_values, var2_neighbor_values, var3_neighbors)

    def interpolate(self):

        r = self.restrict()

        var1_neighbor_values = r[0]
        var2_neighbor_values = r[1]
        var3_neighbor_values = r[2]

        self.var1_neighbors = var1_neighbor_values
        self.var2_neighbors = var2_neighbor_values
        self.var3_neighbors = var3_neighbor_values

        self.p1 = var1_neighbor_values[0]
        self.p2 = var1_neighbor_values[1]
        self.s11 = var2_neighbor_values[0]
        self.s12 = var2_neighbor_values[1]
        self.s21 = var2_neighbor_values[2]
        self.s22 = var2_neighbor_values[3]
        self.u11 = var3_neighbor_values[0]
        self.u12 = var3_neighbor_values[1]
        self.u21 = var3_neighbor_values[2]
        self.u22 = var3_neighbor_values[3]

        if self.s11 == self.s12:
            self.u1 = self.u11
        else:
            self.u1 = self.linear_interpolate(x1=self.s11, x2=self.s12, x=self.var2, q1=self.u11, q2=self.u12)
        if self.s21 == self.s22:
            self.u2 = self.u21
        else:
            self.u2 = self.linear_interpolate(x1=self.s21, x2=self.s22, x=self.var2, q1=self.u21, q2=self.u22)
        if self.p1 == self.p2:
            u = self.u1
        else:
            u = self.linear_interpolate(x1=self.p1, x2=self.p2, x=self.var1, q1=self.u1, q2=self.u2)

        # print("***************")
        # print(self.var1)
        # print(self.var2)
        # print(self.u1, self.u2)
        # print("var1 neighbors: {}".format(var1_neighbor_values))
        # print("var2 neighbors: {}".format(var2_neighbor_values))
        # print("var3 neighbors: {}".format(var3_neighbor_values))
        # print(u)

        return u


class NearestNeighbor1D:

    def __init__(self):
        pass

    def __distance(self, given_val, test_val):
        return abs(given_val - test_val)

    def neighbor(self, given_val, array):
        m = min(enumerate([self.__distance(given_val=given_val, test_val=i) for i in array]), key=itemgetter(1))[0]
        return array[m]

    def neighbor_index(self, given_val, array):
        neighbor = self.neighbor(given_val=given_val, array=array)
        return array.index(neighbor)
