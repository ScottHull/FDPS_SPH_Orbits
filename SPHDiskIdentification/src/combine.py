import os
import csv
from random import randint
import pandas as pd


class CombinedFile:
    """
    Class for combining SPH outputs from different processes.
    """

    def __init__(self, path: str, iteration: int, number_of_processes: int, to_fname: str):
        self.total_particles = None
        self.sim_time = None
        self.file_format = "results.{}_{}_{}.dat"
        self.path = path
        self.iteration = iteration
        self.number_of_processes = number_of_processes
        self.to_fname = to_fname
        self.current_process = 0

    def __get_filename(self):
        """
        Return the st
        :return:
        """
        return self.path + "/" + self.file_format.format(str(self.iteration).zfill(5),
                                                         str(self.number_of_processes).zfill(5),
                                                         str(self.current_process).zfill(5))

    def __read_sph_file(self):
        df = pd.read_csv(self.__get_filename(), sep='\t', skiprows=2, header=None)
        return df

    def format_time(self, time):
        """
        Takes the SPH time given in seconds and converts it to hours with 2 decimal places.
        """
        return round(time / 3600, 2)

    def write_combined_file(self):
        dfs = []
        total_N = 0
        time = 0
        for proc in range(0, self.number_of_processes, 1):
            self.current_process = proc
            with open(self.__get_filename(), 'r') as infile:
                reader = csv.reader(infile, delimiter="\t")
                time = float(next(reader)[0])
                total_N += int(next(reader)[0])
                infile.close()
            dfs.append(self.__read_sph_file())
        merged_df = pd.concat(dfs)
        merged_df.to_csv(self.to_fname.format(self.iteration), index=False, header=False, sep='\t')
        tmp_fname = "temp{}.dat".format(randint(0, int(1e12)))
        temp = open(tmp_fname, 'w')
        temp.write("{}\n{}\n".format(time, total_N))
        with open(self.to_fname.format(self.iteration)) as infile:
            for line in infile:
                temp.write(line)
        infile.close()
        temp.close()
        os.remove(self.to_fname.format(self.iteration))
        os.rename(tmp_fname, self.to_fname.format(self.iteration))
        self.sim_time = self.format_time(time)
        self.total_particles = total_N

    def combine_to_memory(self):
        print("Combining to memory...")
        dfs = []
        total_N = 0
        time = 0
        for proc in range(0, self.number_of_processes, 1):
            self.current_process = proc
            with open(self.__get_filename(), 'r') as infile:
                reader = csv.reader(infile, delimiter="\t")
                time = float(list(next(reader))[0])
                total_N += int(list(next(reader))[0])
                infile.close()
            dfs.append(self.__read_sph_file())
        self.sim_time = self.format_time(time)
        merged_df = pd.concat(dfs)
        print(f"Done combining to memory.  {total_N} particles found in {self.number_of_processes} processes.")
        return merged_df
