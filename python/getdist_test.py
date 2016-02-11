from getdist import plots

dirs = [
        "../completed_runs/standard_Mnu_nEff_mn_build/",
        "../completed_runs/standard_Mnu_nEff_bao_mn_build/",
        "../completed_runs/standard_Mnu_nEff_lrg_mn_build/",
        "../completed_runs/standard_Mnu_nEff_wigglez_mn_build"
       ]

roots = [
        "standard_neff_summnu_planck",
        "standard_neff_summnu_planck_bao_",
        "standard_neff_summnu_planck_lrg_",
        "standard_neff_summnu_planck_wigglez_"
        ]

g = plots.getSinglePlotter(chain_dir=dirs)

g.plot_2d(roots, "h", "ombh2", filled=True, lims=[0.55, 0.80, 0.020, 0.024])
g.settings.legend_fontsize = 12
g.add_legend(["Planck2015", "Planck2015+BAO", "Planck2015+LRG", "Planck2015+WiggleZ"], legend_loc="lower right")

g.export("test.pdf")
