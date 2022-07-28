import os

import pandas as pd
import numpy as np
from tqdm import tqdm

FORCES = ["LO", "NLO", "N2LO", "N3LO", "N4LO", "N4LO+"]
# src = "../../../Deuteron/PRODUCTION1_DEUTERON_2022/output/"
src = "./Deuteron/data/"


def read_data(fname):
    df = read_av18(fname) if "av18" in fname else read_chiral(fname)
    if fname.startswith("0."):
        df["WAVE"] = "1NC"
    else:
        df["WAVE"] = "SIEGERT"
    return df


def read_chiral(fname):
    with open(src + fname, "r") as f:
        cols = f.readlines()[int("dat1" in fname)]
        cols = [c.strip() for c in cols.split("#")[-1].split(",")]
    df = pd.read_csv(src + fname, delim_whitespace=True,
                     skipinitialspace=True, engine='python',
                     comment="#", skip_blank_lines=True,
                     names=cols)
    df["angle"] = (df["THCM"]/np.pi*180).round(2)
    df["CUTOFF"] = 350 + 50*int(fname.split(".")[3])
    ostat = int(fname.split(".")[2])
    df["FORCE"] = FORCES[ostat]
    return df


def read_av18(fname):
    with open(src + fname, "r") as f:
        cols = f.readlines()[int("dat1" in fname)]
        cols = [c.strip() for c in cols.split("#")[-1].split(",")]
    df = pd.read_csv(src + fname, delim_whitespace=True,
                     skipinitialspace=True, engine='python',
                     comment="#", skip_blank_lines=True,
                     names=cols)
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

        df = pd.merge(left=df1, right=df2.drop(columns=["THCM"]),
                      on=["angle", "CUTOFF", "FORCE", "WAVE"])
        df = pd.merge(left=df, right=df3.drop(columns=["THCM"]),
                      on=["angle", "CUTOFF", "FORCE", "WAVE"])
        df = pd.merge(left=df, right=df4.drop(columns=["THCM"]),
                      on=["angle", "CUTOFF", "FORCE", "WAVE"])
        df = pd.merge(left=df, right=df5.drop(columns=["THCM"]),
                      on=["angle", "CUTOFF", "FORCE", "WAVE"])
        df["Energy"] = energy
        df_all = df_all.append(df, ignore_index=True)

    df_all.to_csv("Deuteron/deuteron_all_data.csv", index=False)

    src = "Deuteron/ExpData/DiffCross/"
    df_exp = pd.DataFrame()
    for f in os.listdir(src):
        if not os.path.isfile(src+f):
            print(f"not a file {src+f}")
            continue
        if not f.endswith(".dat"):
            print(f"wrong file {src+f}")
            continue
        df1 = pd.read_csv(src+f, delim_whitespace=True, comment="#",
                          skip_blank_lines=True,
                          header=None, engine='python',
                          names=["angle", "CROSS", "error"])
        df1["fname"] = f
        for energy in [30, 100, 140]:
            if str(energy) in f:
                df1["energy"] = energy
                break
        df_exp = pd.concat((df_exp, df1), ignore_index=True)
    df_exp.to_csv("Deuteron/deuteron_exp_diffcross.csv", index=False)

    src = "./Deuteron/ExpData/TotalCross/"
    df_exp_tot = pd.DataFrame()
    for f in os.listdir(src):
        if not os.path.isfile(src+f):
            print(f"not a file {src+f}")
            continue
        if not f.endswith(".csv"):
            print(f"wrong file {src+f}")
            continue
        df1 = pd.read_csv(src+f, comment="#",
                          skip_blank_lines=True,engine='python')
        df1["fname"] = f
        df_exp_tot = pd.concat((df_exp_tot, df1), ignore_index=True)
    df_exp_tot.to_csv("Deuteron/deuteron_exp_totcross.csv", index=False)
