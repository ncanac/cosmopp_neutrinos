from getdist import plots
from getdist import loadMCSamples

dirs = [
        "../completed_runs/standard_Mnu_nEff_mn/",
        "../completed_runs/standard_Mnu_nEff_bao_mn/",
        "../completed_runs/standard_Mnu_nEff_lrg_mn/",
        "../completed_runs/standard_Mnu_nEff_wigglez_mn/"#,
        #"../completed_runs/standard_Mnu_nEff_szclusters_mn/"
        ]

roots = [
        "standard_neff_summnu_planck",
        "standard_neff_summnu_planck_bao",
        "standard_neff_summnu_planck_lrg",
        "standard_neff_summnu_planck_wigglez"#,
        #"standard_neff_summnu_planck_wigglez_sz"
        ]

samples = []
for root in [d+r for d, r in zip(dirs, roots)]:
    samples.append(loadMCSamples(root))
    #samples.append(MCSamples(root=root, names=names, labels=labels))

for samp in samples:
    p = samp.getParams()
    samp.addDerived((p.ombh2+p.omch2)/p.h**2, name="omm", label="\Omega_M")

g = plots.getSinglePlotter()

g.plot_2d(samples, "omm", "sigma8", filled=True)#, lims=[0.20, 0.40, 0.55, 0.80])
g.settings.legend_fontsize = 12
g.add_legend(["Planck15", "Planck15+BAO", "Planck15+LRG", "Planck15+WiggleZ"], legend_loc="lower left")

g.export("sigma8_v_om_nEff+Mnu.pdf")
