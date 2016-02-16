from getdist import plots
from getdist import loadMCSamples

dirs = [
        "../completed_runs/standard_pps/standard_Mnu_nEff_mn/",
        "../completed_runs/standard_pps/standard_Mnu_nEff_bao_mn/",
        "../completed_runs/standard_pps/standard_Mnu_nEff_lrg_mn/",
        "../completed_runs/standard_pps/standard_Mnu_nEff_wigglez_mn/",
        "../completed_runs/standard_pps/standard_Mnu_nEff_szclusters_mn/"
        ]

roots = [
        "standard_neff_summnu_planck",
        "standard_neff_summnu_planck_bao",
        "standard_neff_summnu_planck_lrg",
        "standard_neff_summnu_planck_wigglez",
        "standard_neff_summnu_planck_szclusters"
        ]

colors=['blue', 'green', 'red', 'magenta', 'gray']
labels = ['Planck15', 'Planck15+BAO', 'Planck15+LRG', 'Planck15+WZ', 'Planck15+Clusters']

# Reorder list for plotting
order = [0, 4, 3, 1, 2]
dirs = [dirs[i] for i in order]
roots = [roots[i] for i in order]
colors = [colors[i] for i in order]
labels = [labels[i] for i in order]

samples = []
for root in [d+r for d, r in zip(dirs, roots)]:
    samples.append(loadMCSamples(root))

#for samp in samples:
#    p = samp.getParams()
#    samp.addDerived((p.ombh2+p.omch2)/p.h**2, name="omm", label="\Omega_M")

# Single plot
#g = plots.getSinglePlotter()
#g.plot_2d(samples, "omm", "sigma8", filled=True)#, lims=[0.20, 0.40, 0.55, 0.80])
#g.settings.legend_fontsize = 12
#g.add_legend(["Planck15", "Planck15+BAO", "Planck15+LRG", "Planck15+WiggleZ", "Planck+Clusters"], legend_loc="lower left")

# Subplots
g = plots.getSubplotPlotter()
g.settings.figure_legend_frame = False
g.settings.legend_frac_subplot_margin = 0.2
g.settings.legend_fontsize = 12
g.settings.fig_width_inch = 5
#g.settings.axes_fontsize = 20
g.rectangle_plot(['sum_mnu', 'sum_mnu'], ['sigma8', 'sigma8'], roots=samples, filled=True,
                colors=colors, param_limits={'sum_mnu': (0, 1.25), 'sigma8': (0.63, 0.97)},
                plot_texts=[['$nKnots=1$', '$nKnots=2$'], ['$nKnots=3$', '$nKnots=4$']],
                legend_labels = labels, legend_ncol=2);

g.export("sigma8_v_Mnu.pdf")
