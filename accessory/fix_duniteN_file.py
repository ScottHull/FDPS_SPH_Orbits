import os
import pandas as pd
import csv

f = "/Users/scotthull/Documents - Scott’s MacBook Pro/PhD Research/Hugoniot/data/duniteN.table_broken.txt"
if os.path.exists("text.txt"):
    os.remove("text.txt")
outfile = open("test.txt", "w")
df = pd.read_fwf(f, sep='\t', header=None, skiprows=1)
df = df.dropna()
for index2, row in df.iterrows():
    for index, i in enumerate(row):
        if type(i) is str:
            if "-" in i and i[0] != "-" and "E-" not in i:
                i = i.replace("-", "E-")
                row[index] = float(i)
                print(i)
    df.loc[index2] = row
df.to_csv("test.txt", sep='\t', header=False, index=False)

# for index, row in df.iterrows():
#     formatted = []
#     for i in row:
#         # if "-" != i[0] and "-" in i and "E-" not in i and "e-" not in i:
#         #     i = i.replace("-", "E-")
#         formatted.append(i)
#     outfile.write("\t".join(i for i in formatted) + "\n")
# outfile.close()