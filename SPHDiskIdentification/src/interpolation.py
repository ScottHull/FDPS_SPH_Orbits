from operator import itemgetter

class NearestNeighbor1D:
    """
    Nearest Neighbor Interpolation for 1D data.
    """

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
