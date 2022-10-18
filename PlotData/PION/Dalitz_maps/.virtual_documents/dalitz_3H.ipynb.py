import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import seaborn as sns


plt.style.use(['science', 'retro'])

mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'DejaVu Sans'
mpl.rcParams['mathtext.it'] = 'DejaVu Sans:italic'
mpl.rcParams['mathtext.bf'] = 'DejaVu Sans:bold'

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "DejaVu Sans"

plt.rc('text', usetex=True)
plt.rcParams['axes.linewidth'] = 1.2


# plt.style.reload_library()
# plt.style.use('retro')
# plt.style.use(['science', 'retro'])
# plt.rcParams['axes.linewidth'] = 1.2
# mpl.rc('text', usetex=True)


df_xy = pd.read_csv("./data_3H_N4LO+_cut=2/for_Dalitz_plot_3h_with_v4_xy2_w.txt",
                    sep=" ", skipinitialspace=True, header=None, index_col=0,
                   names=["ELAB1", "ELAB2", "ELAB3", "x", "y", "val"])
df_e = pd.read_csv("./data_3H_N4LO+_cut=2/for_Dalitz_plot_3h_v4_E1E2.txt",
                    sep=" ", skipinitialspace=True, header=None, index_col=0,
                   names=["ELAB1", "ELAB2", "ELAB3","val"])
df_xy_pw = pd.read_csv("./data_3H_N4LO+_cut=2/for_Dalitz_plot_3h_with_v4_xy2_w_PWIAS.txt",
                    sep=" ", skipinitialspace=True, header=None, index_col=0,
                   names=["ELAB1", "ELAB2", "ELAB3", "x", "y", "val"])
df_e_pw = pd.read_csv("./data_3H_N4LO+_cut=2/for_Dalitz_plot_3h_v4_E1E2_PWIAS.txt",
                    sep=" ", skipinitialspace=True, header=None, index_col=0,
                   names=["ELAB1", "ELAB2", "ELAB3","val"])
df_xy_1nc = pd.read_csv("./data_3H_N4LO+_cut=2/for_Dalitz_plot_3h_with_v4_xy2_w_1NC.txt",
                    sep=" ", skipinitialspace=True, header=None, index_col=0,
                   names=["ELAB1", "ELAB2", "ELAB3", "x", "y", "val"])
df_e_1nc = pd.read_csv("./data_3H_N4LO+_cut=2/for_Dalitz_plot_3h_v4_E1E2_1NC.txt",
                    sep=" ", skipinitialspace=True, header=None, index_col=0,
                   names=["ELAB1", "ELAB2", "ELAB3","val"])
df_xy_2nf = pd.read_csv("./data_3H_N4LO+_cut=2/for_Dalitz_plot_3h_xy2_w.txt",
                    sep=" ", skipinitialspace=True, header=None, index_col=0,
                   names=["ELAB1", "ELAB2", "ELAB3", "x", "y", "val"])
df_e_2nf = pd.read_csv("./data_3H_N4LO+_cut=2/for_Dalitz_plot_3h_E1E2.txt",
                    sep=" ", skipinitialspace=True, header=None, index_col=0,
                   names=["ELAB1", "ELAB2", "ELAB3","val"])


newcolors= np.array([[224,224,224], [153, 204, 255], [102, 255, 255],
                     [51, 153, 255], [0, 128, 255], [0, 102, 204],
                     [0, 51, 102], [0, 25, 51], [25, 0, 51], [0, 0, 0]])
bounds_x_pw = np.array([1.6E-9, 2.0E-09, 4.0E-9, 7.0E-09, 1.0E-08, 5.0E-08, 1.0E-07, 5.0E-07, 1.0E-06])*1e22
bounds_x_1nc = np.array([0.3E-9, 0.4E-9, 0.6E-9, 1.0E-9, 2.0E-09, 5.0E-09, 1.0E-08, 3.0E-08, 6.0E-08])*1e22
# bounds_x_1nc = np.array([0.3E-9, 0.6E-9, 0.9E-9, 2.0E-9, 5.0E-09, 8.0E-09, 1.0E-08, 3.0E-08, 6.0E-08])*1e22
# bounds_x_2nf = bounds_x_pw
# bounds_x_2nf = np.array([2.0E-08, 3.0E-08, 5.0E-08, 8.0E-08, 2.0E-07, 5.0E-07, 8.0E-07, 2.0E-06, 5.0E-06])*1e21
# bounds_x_2nf = np.array([1.3E-09, 2.0E-09, 3.0E-08, 8.0E-08, 2.0E-07, 5.0E-07, 8.0E-07, 2.0E-06, 5.0E-06])*1e21
bounds_x = np.array([1.3E-09, 2.0E-09, 5.0E-09, 1.0E-08, 5.0E-08, 1.0E-07, 2.0E-07, 4.0E-07, 1.0E-06])*1e22
# bounds_x = bounds_x_2nf
bounds_x_2nf = bounds_x
bounds_x_pw = bounds_x


# bounds_e = np.array([5.0E-08, 1.0E-07, 4.0E-07, 7.0E-07, 1.0E-06, 4.0E-06, 7.0E-06, 1.0E-05, 2.0E-05])*1e22
cmap = mpl.colors.ListedColormap(newcolors/255)
# cmap = mpl.colors.Colormap("CMRmap_r")
norm_x = [mpl.colors.BoundaryNorm(b, cmap.N, extend='both') for b in [bounds_x_pw, bounds_x_1nc, bounds_x_2nf, bounds_x]]


df_e_pw.sort_values("val", ascending=False).head(1)*1e22


df_e_1nc.sort_values("val", ascending=False).head(1)*1e22


df_e_2nf.sort_values("val", ascending=False).head(1)*1e22


df_e.sort_values("val", ascending=False).head(1)*1e22


fig = plt.figure(figsize=(7.5,7))
fig.suptitle(r"d$^2\Gamma_\text{nnn}$/dxdy[s$^\text{-1}$] (N4LO+ with $\Lambda$=450 [MeV])", size=14)
fig.supxlabel("x", y=0.05, size=12)
fig.supylabel("y", x=0.03, size=12)

for i, (df, title) in enumerate(zip([df_xy_pw, df_xy_1nc, df_xy_2nf, df_xy],
                                  ["PWIAS", "1NC only", "2NF only", "2NF + 3NF"])):
    ax = plt.subplot(221 + i)
    p = plt.scatter(df.x, df.y, c=df.val*1e22, cmap=cmap, norm=norm_x[i], marker="s", s=6.5, edgecolors='none')
    ax.set_aspect('equal', adjustable='box')
    # plt.title(title)
    plt.minorticks_on()
    plt.xticks(ticks=np.arange(-1, 1.1, 0.5))
    plt.yticks(ticks=np.arange(-1, 1.1, 0.5))
    plt.tick_params(which='minor', direction='in', length=4, width=0.5, top=True, right=True)
    plt.tick_params(which='major', direction='in', length=7, width=1, top=True, right=True)
    plt.grid(which='major', alpha=0.3)
    cbar = plt.colorbar(p, ax=ax, shrink=0.8, pad=0.03)
    cbar.ax.get_yaxis().get_offset_text().set_x(3.0)
    cbar.ax.minorticks_off()
    # if i < 2:
    #     ax.set_xticklabels([])
    # if i % 2 get_ipython().getoutput("= 0:")
    #     ax.set_yticklabels([])
        # cbar.set_label('$d^2\Gamma_{pnn}/dxdy[s^\text{-1}]$')
plt.tight_layout()
# plt.savefig("./figures/Dalitz_map_nnn_xy.pdf", dpi=600, facecolor="white")


bounds_e_pw = np.array([3.0E-08, 5.0E-08, 1.0E-07, 5.0E-07,
                        1.0E-06, 4.0E-06, 7.0E-06, 1.0E-05,
                        1.2E-05])*1e22
bounds_e_1nc = np.array([6.0E-09, 1.0E-8, 1.5E-8, 3.0E-08,
                         6.0E-08, 1.0E-07, 7.0E-07, 1.0E-06,
                        4.0E-06])*1e22
bounds_e_2nf = np.array([5.0E-08, 1.0E-07, 4.0E-07, 7.0E-07,
                         1.0E-06, 4.0E-06, 7.0E-06, 1.0E-05,
                        1.3E-05])*1e22
# bounds_e_2nf = bounds_e_pw
# bounds_e = np.array([5.0E-07, 7.0E-07, 1.0E-06, 4.0E-06, 7.0E-06, 1.0E-05, 5.0E-05, 1.0E-04, 3.0E-04])*1e22
bounds_e = bounds_e_pw
bounds_e_2nf = bounds_e_pw


norm_e = [mpl.colors.BoundaryNorm(b, cmap.N, extend='both') for b in [bounds_e_pw, bounds_e_1nc, bounds_e_2nf, bounds_e]]


fig = plt.figure(figsize=(7.5,7))
fig.suptitle(r"d$^2\Gamma_\text{nnn}$/dE$_\text{1}$dE$_\text{2}$[fm$^\text{2}$~s$^\text{-1}$] (N4LO+ with $\Lambda$=450 [MeV])", size=14)
fig.supxlabel(r"E$_\text{1}$ [MeV]", y=0.05, size=12)
fig.supylabel(r"E$_\text{2}$ [MeV]", x=0.03, size=12)

for i, (df, title) in enumerate(zip([df_e_pw, df_e_1nc, df_e_2nf, df_e],
                                  ["PWIAS", "1NC only", "2NF only", "2NF + 3NF"])):
    ax = plt.subplot(221 + i)
    p = plt.scatter(df.ELAB1, df.ELAB2, c=df.val*1e22, cmap=cmap, norm=norm_e[i], marker="s", s=6.5, edgecolors='none')
    ax.set_aspect('equal', adjustable='box')
    # plt.title(title)
    plt.minorticks_on()
    # plt.xticks(ticks=np.arange(-1, 1.1, 0.5))
    # plt.yticks(ticks=np.arange(-1, 1.1, 0.5))
    plt.tick_params(which='minor', direction='in', length=4, width=0.5, top=True, right=True)
    plt.tick_params(which='major', direction='in', length=7, width=1, top=True, right=True)
    plt.grid(which='major', alpha=0.3)
    cbar = plt.colorbar(p, ax=ax, shrink=0.8, pad=0.03)
    cbar.formatter.set_powerlimits([16,16])
    cbar.ax.get_yaxis().get_offset_text().set_x(3.0)
    # cbar.ax.minorticks_off()
    # if i < 2:
    #     ax.set_xticklabels([])
    # if i % 2 get_ipython().getoutput("= 0:")
    #     ax.set_yticklabels([])
        # cbar.set_label('$d^2\Gamma_{pnn}/dxdy[s^\text{-1}]$')
plt.tight_layout()
plt.savefig("./figures/Dalitz_map_nnn_E1E2.pdf", dpi=600, facecolor="white")


df_en = pd.read_csv("./data_3H_N4LO+_cut=2/for_Dalitz_plot_3h_with_v4_dGammadEn.txt",
                    sep=" ", skipinitialspace=True, header=None, index_col=None, skiprows=1,
                    names=["energy", "dGdq", "dGdE"])
df_en_pw = pd.read_csv("./data_3H_N4LO+_cut=2/for_Dalitz_plot_3h_with_v4_dGammadEn_PWIAS.txt",
                    sep=" ", skipinitialspace=True, header=None, index_col=None, skiprows=1,
                    names=["energy", "dGdq", "dGdE"])
df_en_1nc = pd.read_csv("./data_3H_N4LO+_cut=2/for_Dalitz_plot_3h_with_v4_dGammadEn_1NC.txt",
                    sep=" ", skipinitialspace=True, header=None, index_col=None, skiprows=1,
                    names=["energy", "dGdq", "dGdE"])
df_en_2nf = pd.read_csv("./data_3H_N4LO+_cut=2/for_Dalitz_plot_3h_dGammadEn.txt",
                    sep=" ", skipinitialspace=True, header=None, index_col=None, skiprows=1,
                    names=["energy", "dGdq", "dGdE"])



plt.figure(figsize=(4,3))
plt.plot(df_en.energy, df_en.dGdE*1e22, label="with 3NF", c="k")
plt.plot(df_en_2nf.energy, df_en_2nf.dGdE*1e22, label="2NF only", c="k", ls="--")
plt.plot(df_en_1nc.energy, df_en_1nc.dGdE*1e22, label="1NC", c="k", ls=":")
plt.plot(df_en_pw.energy, df_en_pw.dGdE*1e22, label="PWIAS", c="k", ls="-.")
# plt.legend()

plt.yscale("log")
plt.tick_params(which='minor', direction='in', length=3, width=0.5, top=True, right=True)
plt.tick_params(which='major', direction='in', length=5, width=1, top=True, right=True)
plt.ylabel(r"d$\Gamma_\text{nnn}$/dE$_\text{n}$[fm/s]")
plt.xlabel(r"E$_\text{n}$[1/fm]")

plt.savefig("./figures/3H_dGdEn.pdf", dpi=600, facecolor="white")
plt.show()


colnames=["r", "dGdr", "DUMMY"]
df_r = pd.read_csv("./data_3H_N4LO+_cut=2/dgammadr_3h_with_v4.txt",
                    sep=" ", skipinitialspace=True, header=None, index_col=None, skiprows=1,
                    names=colnames)
df_r_pw = pd.read_csv("./data_3H_N4LO+_cut=2/dgammadr_3h_with_v4_PWIAS.txt",
                    sep=" ", skipinitialspace=True, header=None, index_col=None, skiprows=1,
                    names=colnames)
df_r_1nc = pd.read_csv("./data_3H_N4LO+_cut=2/dgammadr_3h_with_v4_1NC.txt",
                    sep=" ", skipinitialspace=True, header=None, index_col=None, skiprows=1,
                    names=colnames)
df_r_2nf = pd.read_csv("./data_3H_N4LO+_cut=2/dgammadr_3h.txt",
                    sep=" ", skipinitialspace=True, header=None, index_col=None, skiprows=1,
                    names=colnames)


plt.figure(figsize=(4,3))
plt.plot(df_r.r, df_r.dGdr*1e22, label="with 3NF", c="k")
plt.plot(df_r_2nf.r, df_r_2nf.dGdr*1e22, label="2NF only", c="k", ls="--")
plt.plot(df_r_1nc.r, df_r_1nc.dGdr*1e22, label="1NC", c="k", ls=":")
plt.plot(df_r_pw.r, df_r_pw.dGdr*1e22, label="PWIAS", c="k", ls="-.")
# plt.legend()

plt.yscale("log")
plt.tick_params(which='minor', direction='in', length=3, width=0.5, top=True, right=True)
plt.tick_params(which='major', direction='in', length=5, width=1, top=True, right=True)
plt.ylabel(r"d$\Gamma_\text{nnn}$/dr[s$^\text{-1}$]")
plt.xlabel(r"r")

plt.savefig("./figures/3H_dGdr.pdf", dpi=600, facecolor="white")
plt.show()


colnames=["phi", "dGdphi", "DUMMY"]
df_phi = pd.read_csv("./data_3H_N4LO+_cut=2/dgammadphi_3h_with_v4.txt",
                    sep=" ", skipinitialspace=True, header=None, index_col=None, skiprows=1,
                    names=colnames)
df_phi_pw = pd.read_csv("./data_3H_N4LO+_cut=2/dgammadphi_3h_with_v4_PWIAS.txt",
                    sep=" ", skipinitialspace=True, header=None, index_col=None, skiprows=1,
                    names=colnames)
df_phi_1nc = pd.read_csv("./data_3H_N4LO+_cut=2/dgammadphi_3h_with_v4_1NC.txt",
                    sep=" ", skipinitialspace=True, header=None, index_col=None, skiprows=1,
                    names=colnames)
df_phi_2nf = pd.read_csv("./data_3H_N4LO+_cut=2/dgammadphi_3h.txt",
                    sep=" ", skipinitialspace=True, header=None, index_col=None, skiprows=1,
                    names=colnames)


plt.figure(figsize=(4,3))
plt.plot(df_phi.phi*180/np.pi, df_phi.dGdphi*1e22, label="with 3NF", c="k")
plt.plot(df_phi_2nf.phi*180/np.pi, df_phi_2nf.dGdphi*1e22, label="2NF only", c="k", ls="--")
plt.plot(df_phi_1nc.phi*180/np.pi, df_phi_1nc.dGdphi*1e22, label="1NC", c="k", ls=":")
plt.plot(df_phi_pw.phi*180/np.pi, df_phi_pw.dGdphi*1e22, label="PWIAS", c="k", ls="-.")
# plt.legend()

plt.yscale("log")
plt.tick_params(which='minor', direction='in', length=3, width=0.5, top=True, right=True)
plt.tick_params(which='major', direction='in', length=5, width=1, top=True, right=True)
plt.ylabel(r"d$\Gamma_\text{nnn}\,^\text{ring}$/d$\Phi$[(s~rad)$^\text{-1}$]")
plt.xlabel(r"$\Phi$ [deg]")
plt.xlim([0, 120])

plt.savefig("./figures/3H_dGdphi.pdf", dpi=600, facecolor="white")
plt.show()



