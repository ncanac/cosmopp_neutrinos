#include <memory>
#include <cstdlib>
#include <string>
#include <sstream>
#include <cmath>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <combined_like.hpp>
#include <power_spectrum.hpp>
#include <neutrino_cosmology_params.hpp>
#include <mn_scanner.hpp>

// argv[1] = c for cubic, anything else for linear
// argv[2] = number of knots
// argv[3+] =   neff to vary n effective
//              sum_mnu to vary sum of neutrino masses
//              planck to use planck data
//              wmap to use wmap data
//              bao to use bao data
//              lrg to use lrg data
int main(int argc, char *argv[])
{
    try {
        StandardException exc;

        bool varyNEff = false;
        bool varySumMNu = false;
        bool usePlanck = false;
        bool useWMAP = false;
        bool useBAO = false;
        bool useLRG = false;

        for(int i = 1; i < argc; ++i)
        {
            if(std::string(argv[i]) == "neff")
                varyNEff = true;
            if(std::string(argv[i]) == "sum_mnu")
                varySumMNu = true;
            if(std::string(argv[i]) == "planck")
                usePlanck = true;
            if(std::string(argv[i]) == "wmap")
                useWMAP = true;
            if(std::string(argv[i]) == "bao")
                useBAO = true;
            if(std::string(argv[i]) == "lrg")
                useLRG = true;
        }

        // omBH2, omCH2, h, tau, ns, as
        int nPar = 6;

        if(varyNEff)
        {
            ++nPar;
            output_screen("Varying N_eff!" << std::endl);
        }
        else
        {
            output_screen("N_eff not being varied. To vary it give \"neff\" as an extra argument." << std::endl);
        }
        
        if(varySumMNu)
        {
            ++nPar;
            output_screen("Varying sum_mnu!" << std::endl);
        }
        else
        {
            output_screen("sum_mnu not being varied. To vary it give \"sum_mnu\" as an extra argument." << std::endl);
        }

        if(usePlanck)
        {
            // for A_planck
            ++nPar;
            output_screen("Using Planck." << std::endl);
        }
        if(useWMAP)
        {
            output_screen("Using WMAP." << std::endl);
        }
        if(useBAO)
        {
            output_screen("Using BAO." << std::endl);
        }
        if(useLRG)
        {
            output_screen("Using LRG." << std::endl);
        }

        // Starting values for cosmological parameters
        const double h = 0.6731;
        const double omBH2 = 0.02222;
        const double omCH2 = 0.1197;
        const double tau = 0.078;
        const double ns = 0.9655;
        const double as = 3.089; // ln(10^10 as), as ~ 2.196e-9
        const double pivot = 0.05;

        const double r = 1e-10;
        const double nt = 0;

        const double nEff = 3.046; 
        const int nMassive = 1;
        const double sumMNu = 0.5;

        StandardPSDegenNuParams params(omBH2, omCH2, h, tau, ns, std::exp(as)/1e10, pivot, nEff, nMassive, sumMNu, varyNEff, varySumMNu);

//#ifdef COSMO_PLANCK_15
//        PlanckLikelihood planckLike(true, true, true, false, true, false, false, false, 100);
//#else
//        ERROR NOT IMPLEMENTED;
//#endif
//        planckLike.setModelCosmoParams(&params);
        CombinedLikelihood like(usePlanck, useWMAP, useBAO, useLRG);
        like.setModelCosmoParams(&params);

        std::stringstream root;
        root << "standard_";
        if(varyNEff)
            root << "_neff";
        if(varySumMNu)
            root << "_summnu";
        if(usePlanck)
            root << "_planck";
        if(useWMAP)
            root << "_wmap";
        if(useBAO)
            root << "_bao";
        if(useLRG)
            root << "_lrg";
        root << "_";
        MnScanner scanner(nPar, like, 300, root.str());

        int paramIndex = 0;

        output_screen("Setting parameters" << std::endl);
        scanner.setParam(paramIndex++, "ombh2", 0.02, 0.025);
        scanner.setParam(paramIndex++, "omch2", 0.1, 0.2);
        scanner.setParam(paramIndex++, "h", 0.55, 0.85);
        scanner.setParam(paramIndex++, "tau", 0.02, 0.20);
        scanner.setParam(paramIndex++, "ns", 0.9, 1.1);
        scanner.setParam(paramIndex++, "as", 2.0, 4.0);
        if(varyNEff)
            scanner.setParam(paramIndex++, "n_eff", 2.0, 5.0);
        if(varySumMNu)
            scanner.setParam(paramIndex++, "sum_mnu", 0.001, 3.0);

        scanner.setParamGauss(paramIndex++, "A_planck", 1, 0.0025);

        check(paramIndex == nPar, "");

        output_screen("Running scan." << std::endl);
        scanner.run();
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
