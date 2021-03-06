#include <memory>
#include <cstdlib>
#include <string>
#include <sstream>
#include <cmath>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <combined_like.hpp>
#include <power_spectrum.hpp>
#include <cosmological_params.hpp>
#include <mn_scanner.hpp>
#include <timer.hpp>

// argv[1+] =   planck to use planck data
//              wmap to use wmap data
//              bao to use bao data
//              lrg to use lrg data
//              wigglez to use wigglez data
int main(int argc, char *argv[])
{
    try {
        StandardException exc;
        if(argc < 2)
        {
            std::string exceptionString = "Need to specify at least one data set to use (e.g., planck, wmap, bao, or lrg).";
            exc.set(exceptionString);
            throw exc;
        }

        bool usePlanck = false;
        bool useWMAP = false;
        bool useBAO = false;
        bool useLRG = false;
        bool useWiggleZ = false;

        for(int i = 1; i < argc; ++i)
        {
            if(std::string(argv[i]) == "planck")
                usePlanck = true;
            if(std::string(argv[i]) == "wmap")
                useWMAP = true;
            if(std::string(argv[i]) == "bao")
                useBAO = true;
            if(std::string(argv[i]) == "lrg")
                useLRG = true;
            if(std::string(argv[i]) == "wigglez")
                useWiggleZ = true;
        }

        // omBH2, omCH2, h, tau, ns, as
        int nPar = 6;

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
        if(useWiggleZ)
        {
            output_screen("Using WiggleZ." << std::endl);
        }

        // Starting values for cosmological parameters
        output_screen("Running with Lambda CDM cosmological parameters." << std::endl);
        const double h = 0.6731;
        const double omBH2 = 0.02222;
        const double omCH2 = 0.1197;
        const double tau = 0.078;
        const double ns = 0.9655;
        const double as = 3.089; // ln(10^10 as), as ~ 2.196e-9
        const double pivot = 0.05;
        // From Planck 2015, Table 3, Column 4
        //const double h = 0.6727;
        //const double omBH2 = 0.02225;
        //const double omCH2 = 0.1198;
        //const double tau = 0.079;
        //const double ns = 0.9645;
        //const double as = 3.094; // ln(10^10*as)
        //const double pivot = 0.05;
        LambdaCDMParams params(omBH2, omCH2, h, tau, ns, std::exp(as)/1e10, pivot);

//#ifdef COSMO_PLANCK_15
//        PlanckLikelihood planckLike(true, true, true, false, true, false, false, false, 100);
//#else
//        ERROR NOT IMPLEMENTED;
//#endif
//        planckLike.setModelCosmoParams(&params);

        Cosmo cosmo;
        cosmo.preInitialize(5000, false, false, false, 0, 100, 1e-6, 1.0);

        std::string datapath = "/Volumes/Data1/ncanac/cosmopp_neutrinos";
        CombinedLikelihood like(datapath, cosmo, usePlanck, useWMAP, useBAO, useLRG, useWiggleZ);
        like.setModelCosmoParams(&params);

        std::stringstream root;
        root << "LCDM";
        if(usePlanck)
            root << "_planck";
        if(useWMAP)
            root << "_wmap";
        if(useBAO)
            root << "_bao";
        if(useLRG)
            root << "_lrg";
        if(useWiggleZ)
            root << "_wigglez";
        //PolyChord pc(nPar, planckLike, 300, root.str(), 3 * (4 + (varyNEff ? 1 : 0) + (varySumMNu ? 1 : 0)));
        MnScanner scanner(nPar, like, 300, root.str());

        int paramIndex = 0;

        output_screen("Setting parameters" << std::endl);
        scanner.setParam(paramIndex++, "ombh2", 0.02, 0.025);
        scanner.setParam(paramIndex++, "omch2", 0.1, 0.2);
        scanner.setParam(paramIndex++, "h", 0.55, 0.85);
        scanner.setParam(paramIndex++, "tau", 0.02, 0.20);
        scanner.setParam(paramIndex++, "ns", 0.9, 1.1);
        scanner.setParam(paramIndex++, "as", 2.0, 4.0);

        if(usePlanck)
            scanner.setParamGauss(paramIndex++, "A_planck", 1, 0.0025);

        check(paramIndex == nPar, "");

        Timer timer("Multinest scan");
        timer.start();
        scanner.run();
        const unsigned long time = timer.end();
        output_screen("Multinest scan took " << time / 1000000 << " seconds." << std::endl);
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
