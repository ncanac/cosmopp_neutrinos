from getdist import plots
from getdist import MCSamples
from getdist import loadMCSamples

dirs = [
        "../completed_runs/standard_Mnu_nEff_mn_build/",
        "../completed_runs/standard_Mnu_nEff_bao_mn_build/",
        "../completed_runs/standard_Mnu_nEff_lrg_mn_build/",
        "../completed_runs/standard_Mnu_nEff_wigglez_mn_build/"
        ]

roots = [
        "standard_neff_summnu_planck",
        "standard_neff_summnu_planck_bao_",
        "standard_neff_summnu_planck_lrg_",
        "standard_neff_summnu_planck_wigglez_"
        ]

names = [
        "ombh2",
        "omch2",
        "h",
        "tau",
        "ns",
        "as",
        "n_eff",
        "sum_mnu",
        "A_planck",
        ]

labels = [
        #"\sigma_8",
        "\Omega_b h^2",
        "\Omega_c h^2",
        "h",
        "\Tau",
        "n_s",
        "a_s",
        "n_{eff}",
        "\Sigma M_{nu}"
        "A_{planck}"
        ]

samples = []
for root in [d+r for d, r in zip(dirs, roots)]:
    samples.append(loadMCSamples(root))
    #samples.append(MCSamples(root=root, names=names, labels=labels))

for samp in samples:
    p = samp.getParams()
    samp.addDerived((p.ombh2+p.omch2)/p.h**2, name="omm", label="\Omega_M")

g = plots.getSinglePlotter()

g.plot_2d(samples, "omm", "ombh2", filled=True)#, lims=[0.20, 0.40, 0.55, 0.80])
g.settings.legend_fontsize = 12
g.add_legend(["Planck15", "Planck15+BAO", "Planck15+LRG", "Planck15+WiggleZ"], legend_loc="lower left")

g.export("sigma8_v_om_nEff+Mnu.pdf")
