import pandas as pd
import csv
import os


class CombineFile:

    def __init__(self, num_processes, time, output_path):
        self.num_processes = num_processes
        self.time = time
        self.sim_time = None
        self.output_path = output_path
        self.file_format = "results.{}_{}_{}.dat"
        self.curr_process = 0

    def __get_filename(self):
        return self.output_path + "/" + self.file_format.format(str(self.time).zfill(5),
                                                                str(self.num_processes).zfill(5),
                                                                str(self.curr_process).zfill(5))

    def __read_sph_file(self):
        df = pd.read_csv(self.__get_filename(), sep='\t', skiprows=2, header=None)
        return df

    def combine(self):
        dfs = []
        total_N = 0
        time = 0
        for proc in range(0, self.num_processes, 1):
            self.curr_process = proc
            with open(self.__get_filename(), 'r') as infile:
                reader = csv.reader(infile, delimiter="\t")
                time = float(next(reader)[0])
                total_N += int(next(reader)[0])
                infile.close()
            dfs.append(self.__read_sph_file())
        merged_df = pd.concat(dfs)
        merged_df.to_csv("merged_{}.dat".format(self.time), index=False, header=False, sep='\t')
        temp = open("temp.dat", 'w')
        temp.write("{}\n{}\n".format(time, total_N))
        with open("merged_{}.dat".format(self.time)) as infile:
            for line in infile:
                temp.write(line)
        infile.close()
        temp.close()
        os.remove("merged_{}.dat".format(self.time))
        os.rename("temp.dat", "merged_{}.dat".format(self.time))

    def combine_df(self):
        dfs = []
        total_N = 0
        time = 0
        for proc in range(0, self.num_processes, 1):
            self.curr_process = proc
            with open(self.__get_filename(), 'r') as infile:
                reader = csv.reader(infile, delimiter="\t")
                time = float(next(reader)[0])
                total_N += int(next(reader)[0])
                infile.close()
            dfs.append(self.__read_sph_file())
        self.sim_time = time
        merged_df = pd.concat(dfs)
        return merged_df
