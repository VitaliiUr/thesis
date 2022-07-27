import os

import pandas as pd
import numpy as np
from tqdm import tqdm

FORCES = ["LO", "NLO", "N2LO", "N3LO", "N4LO", "N4LO+"]
src = "Deuteron/data/"


def read_data(fname):
    df = read_av18(fname) if "av18" in fname else read_chiral(fname)
    if fname.startswith("0."):
        df["WAVE"] = "1NC"
    else:
        df["WAVE"] = "SIEGERT"
    return df


def read_chiral(fname):
    add = 0
    skpfoot = 0
    if "dat1" in fname:
        add = 1
        skpfoot = 6
    with open(src + fname, "r") as f:
        cols = f.readlines()[0+add]
        cols = [c.strip() for c in cols.split("#")[-1].split(",")]
    df = pd.read_csv(src + fname, skiprows=1+add, delimiter=" ",
                     skipinitialspace=True, skipfooter=skpfoot, engine='python')
    df.columns = cols
    df["angle"] = (df["THCM"]/np.pi*180).round(2)
    df["CUTOFF"] = 350 + 50*int(fname.split(".")[3])
    ostat = int(fname.split(".")[2])
    df["FORCE"] = FORCES[ostat]
    # if FORCES[ostat] == "N4LO+":
    #     print (fname, 350 + 50*int(fname.split(".")[3]), FORCES[ostat])
    return df


def read_av18(fname):
    add = 0
    skpfoot = 0
    if "dat1" in fname:
        add = 1
        skpfoot = 6
    with open(src + fname, "r") as f:
        cols = f.readlines()[0+add]
        cols = [c.strip() for c in cols.split("#")[-1].split(",")]
    df = pd.read_csv(src + fname, skiprows=1+add, delimiter=" ",
                     skipinitialspace=True, skipfooter=skpfoot, engine='python')
    df.columns = cols
    df["angle"] = (df["THCM"]/np.pi*180).round(2)
    df["CUTOFF"] = 0
    df["FORCE"] = "AV18"
    return df

if __name__ == "__main__":
    energies = set(map(lambda x: int(x.split(".")[-1].rstrip("mev")),
                       [f for f in os.listdir(src) if f.endswith("mev")]))
    print(energies)
    df_all = pd.DataFrame()
    for energy in tqdm(energies):
        data = [f for f in os.listdir(src) if f.endswith(
            f".{energy}mev")]
        df1 = pd.concat([read_data(fname) for fname in data if "dat1" in fname])
        df2 = pd.concat([read_data(fname) for fname in data if "dat2" in fname])
        df3 = pd.concat([read_data(fname) for fname in data if "dat3" in fname])
        df4 = pd.concat([read_data(fname) for fname in data if "dat4" in fname])
        df5 = pd.concat([read_data(fname) for fname in data if "dat5" in fname])
        
        df = pd.merge(left=df1, right=df2.drop(columns=["THCM"]), on=["angle", "CUTOFF", "FORCE", "WAVE"])
        df = pd.merge(left=df, right=df3.drop(columns=["THCM"]), on=["angle", "CUTOFF", "FORCE", "WAVE"])
        df = pd.merge(left=df, right=df4.drop(columns=["THCM"]), on=["angle", "CUTOFF", "FORCE", "WAVE"])
        df = pd.merge(left=df, right=df5.drop(columns=["THCM"]), on=["angle", "CUTOFF", "FORCE", "WAVE"])
        df["Energy"] = energy
        df_all = df_all.append(df, ignore_index=True)

    df_all.to_csv("deuteron_all_data.csv", index=False)

    src = "Deuteron/ExpData/"
    df_exp = pd.DataFrame()
    for f in os.listdir(src):
        if not os.path.isfile(src+f) or not f.endswith(".dat"):
            print(f"not a file {src+f}")
            continue
        print(f)
        df1 = pd.read_csv(src+f, delim_whitespace=True, comment="#",
                        skip_blank_lines=True,
                        header=None, engine='python',
                        names=["angle", "value", "error"])
        df1["fname"] = f
        for energy in [30, 100, 140]:
            if str(energy) in f:
                df1["energy"] = energy
                break
        df_exp = df_exp.append(df1)
    df_exp.to_csv("deuteron_all_exp.csv", index=False)
