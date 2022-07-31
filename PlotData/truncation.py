import numpy as np
import pandas as pd


FORCES = ["LO", "NLO", "N2LO", "N3LO", "N4LO", "N4LO+"]
CUTOFF = [400, 450, 500, 550]

FMM1 = 197.327
MN = 939.5653
MP = 938.2720
MM = 2.0*MP*MN/(MP+MN)
EBDEUTALL = np.array([-2.1189454304609652, -2.1830488263808228,
                      -2.1998863525318981, -2.2232629900695660,
                      -2.2232629164870863, -2.2232628519800381])
PIONMASS = 139.570  # MeV

def get_truncation(df, energy, cutoff=450,
                   observable="CROSS2", wave="SIEGERT"):

    ELABNN = energy-np.abs(EBDEUTALL)
    ECMNN = ELABNN-energy**2/(4.0*MM)
    P0MEV = np.sqrt(ECMNN*MM)
    # denom = cutoff
    denom = 600

    # EPSILON = max([P0MEV[-1]/cutoff, PIONMASS/cutoff]) # ??????????
    EPSILON = max([P0MEV[-1]/denom, PIONMASS/denom]) # ??????????

    df_tmp = df[(df.CUTOFF == cutoff) &
                (df.Energy == energy) &
                (df.WAVE == wave)]
    df_tmp.sort_values("angle", inplace=True)

    df_piv = df_tmp.pivot(index="angle", columns="FORCE")[observable]
    df_diff = pd.DataFrame(columns=FORCES)
    orders = np.arange(len(FORCES))
    orders[1:] += 1

    for i, (order, force) in enumerate(zip(orders, FORCES)):
        maxdiff = df_piv[FORCES[i:]].max(axis=1) - \
            df_piv[FORCES[i:]].min(axis=1)
        base = df_piv["LO"]*EPSILON**(order+1+int(order == 0))
        diff = []
        maxdiff2 = []
        for j in range(2, order+1):
            diff.append((df_piv[FORCES[np.where(orders == j)[0][0]]] -
                        df_piv[FORCES[np.where(orders == j)[0][0]-1]])
                        * EPSILON**(order+1-j))
        for j in range(i):
            maxdiff2.append(df_diff[FORCES[j]]*EPSILON**(i-j))
        # print(order, force)
        df_diff[force] = np.max((base, maxdiff, *diff, *maxdiff2), axis=0)

    return df_diff
