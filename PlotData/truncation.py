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


def get_truncation(df, energy, Lambda=650):
    """Calculate truncation error

    Parameters
    ----------
    df : pd.DataFrame
        pandas data frame with index - x value of observables (angle) and each columns is 
        a prediction in subsequent chiral orders
    energy : float
        Energy of the incoming particle (photon)
    Returns
    -------
    pd.DataFrame
        Truncation errors for each order
    """
    FORCES = df.columns
    ELABNN = energy-np.abs(EBDEUTALL)
    ECMNN = ELABNN-energy**2/(4.0*MM)
    P0MEV = np.sqrt(ECMNN*MM)

    EPSILON = max([P0MEV[-1]/Lambda, PIONMASS/Lambda]) # ??????????
    print(EPSILON)
    df_piv = df.sort_index()

    # df_piv = df_tmp.pivot(index="angle", columns="FORCE")[observable]
    # \delta X - truncation
    df_diff = pd.DataFrame(columns=FORCES)
    orders = np.arange(len(FORCES))
    orders[1:] += 1

    for force_num, (order, force) in enumerate(zip(orders, FORCES)):
        # Q^{i+1} X^{(0)}
        base = df_piv["LO"]*EPSILON**(order+1+int(order == 0))
        # Q^{i+1-j} \Delta X^{(j)}
        diff = []
        for j in range(2, order+1):
            diff.append((df_piv[FORCES[np.where(orders == j)[0][0]]] -
                        df_piv[FORCES[np.where(orders == j)[0][0]-1]])
                        * EPSILON**(order+1-j))
        # maximal difference max(|X^{j>=i} - X^{k>=i}|)
        maxdiff = df_piv[FORCES[force_num:]].max(axis=1) - \
            df_piv[FORCES[force_num:]].min(axis=1)
        # maxdiff2 = []
        # for j in range(i):
        #     maxdiff2.append(df_diff[FORCES[j]]*EPSILON**(i-j))
        # \delta X^(i) > Q \delta X^(i-1)
        if order >= 2:
            maxdiff.append(df_diff[FORCES[force_num-1]]*EPSILON)
        df_diff[force] = np.max((base, maxdiff, *diff), axis=0)

    return df_diff
