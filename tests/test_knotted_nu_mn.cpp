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
#include <timer.hpp>

// argv[1] = c for cubic, anything else for linear
// argv[2] = number of knots
// argv[3+] =   neff to vary n effective
//              sum_mnu to vary sum of neutrino masses
//              planck to use planck data
//              wmap to use wmap data
//              bao to use bao data
//              lrg to use lrg data
//              wigglez use to wigglez data
int main(int argc, char *argv[])
{
    try {
        StandardException exc;
        if(argc < 3)
        {
            std::string exceptionString = "Need to specify the spline type (linear or cubic) and the number of (internal) knots.";
            exc.set(exceptionString);
            throw exc;
        }

        bool isLinear = true;
        if(argv[1][0] == 'c')
            isLinear = false;

        const int nKnots = std::atoi(argv[2]);

        if(nKnots < 0 || nKnots > 100)
        {
            std::stringstream exceptionStr;
            exceptionStr << "The number of knots " << nKnots << " is invalid. Needs to be positive and not larger than 100.";
            exc.set(exceptionStr.str());
            throw exc;
        }

        bool varyNEff = false;
        bool varySumMNu = false;
        bool usePlanck = false;
        bool useWMAP = false;
        bool useBAO = false;
        bool useLRG = false;
        bool useWiggleZ = false;

        for(int i = 3; i < argc; ++i)
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
            if(std::string(argv[i]) == "wigglez")
                useWiggleZ = true;
        }

        int nPar = 4 + 2 * (nKnots + 2) - 2;

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
        if(useWiggleZ)
        {
            output_screen("Using WiggleZ." << std::endl);
        }

        const double h = 0.702;
        const double omBH2 = 0.02262;
        const double omCH2 = 0.1138;
        const double tau = 0.088;
        const double as = 2.1955e-9;
        const double ns = 0.9655;
        const double pivot = 0.05;

        const double nEff = 3.046; 
        const int nMassive = (varySumMNu ? 1 : 0);
        const double sumMNu = 0.0;

        const double kMin = 0.8e-6;
        const double kMax = 12.0;
        const double aMin = -2;
        const double aMax = 4;

        std::vector<double> kVals(nKnots + 2);
        std::vector<double> amplitudes(nKnots + 2);

        kVals[0] = kMin;
        kVals.back() = kMax;

        const double deltaLogK = (std::log(kMax) - std::log(kMin)) / (nKnots + 1);

        for(int i = 1; i < kVals.size() - 1; ++i)
            kVals[i] = std::exp(std::log(kMin) + i * deltaLogK);

        for(int i = 0; i < amplitudes.size(); ++i)
            amplitudes[i] = as * pow(kVals[i]/pivot, ns - 1.0);

        SplineWithDegenerateNeutrinosParams params(isLinear, omBH2, omCH2, h, tau, kVals, amplitudes, nEff, nMassive, sumMNu, varyNEff, varySumMNu);

//#ifdef COSMO_PLANCK_15
//        PlanckLikelihood planckLike(true, true, true, false, true, false, false, false, 100);
//#else
//        ERROR NOT IMPLEMENTED;
//#endif
//        planckLike.setModelCosmoParams(&params);

        //Cosmo cosmo;
        //cosmo.preInitialize(5000, false, true, false, 0, 100, kMin/0.8, kMax/1.2);

        std::string datapath = "/Volumes/Data1/ncanac/cosmopp_neutrinos";
        CombinedLikelihood like(datapath, usePlanck, useWMAP, useBAO, useLRG, useWiggleZ);
        like.setModelCosmoParams(&params);

        std::stringstream root;
        root << "knotted_";
        if(isLinear)
            root << "linear_";
        else
            root << "cubic_";
        root << nKnots;
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
        if(useWiggleZ)
            root << "_wigglez";
        root << "_";
        MnScanner scanner(nPar, like, 300, root.str());

        int paramIndex = 0;

        output_screen("Setting parameters" << std::endl);
        scanner.setParam(paramIndex++, "ombh2", 0.02, 0.025);
        scanner.setParam(paramIndex++, "omch2", 0.1, 0.14);
        scanner.setParam(paramIndex++, "h", 0.55, 0.80);
        scanner.setParam(paramIndex++, "tau", 0.04, 0.12);
        if(varyNEff)
            scanner.setParam(paramIndex++, "n_eff", 2.0, 5.0);
        if(varySumMNu)
            scanner.setParam(paramIndex++, "sum_mnu", 0.01, 2.0);

        for(int i = 1; i < kVals.size() - 1; ++i)
        {
            std::stringstream paramName;
            paramName << "k_" << i;
            scanner.setParam(paramIndex++, paramName.str(), std::log(kMin), std::log(kMax));
        }

        for(int i = 0; i < amplitudes.size(); ++i)
        {
            std::stringstream paramName;
            paramName << "a_" << i;
            scanner.setParam(paramIndex++, paramName.str(), aMin, aMax);
        }

        if(usePlanck)
            scanner.setParamGauss(paramIndex++, "A_planck", 1, 0.0025);

        check(paramIndex == nPar, "");

        output_screen("Running scan." << std::endl);
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
