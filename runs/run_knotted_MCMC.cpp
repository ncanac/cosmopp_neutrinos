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
#include <mcmc.hpp>

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

        const double kMin = 0.8e-6;
        const double kMax = 1.2;
        const double aMin = -2;
        const double aMax = 4;

        std::vector<double> kVals(nKnots + 2);
        std::vector<double> amplitudes(nKnots + 2);

        kVals[0] = kMin;
        kVals.back() = kMax;

        const double deltaLogK = (std::log(kMax) - std::log(kMin)) / (nKnots + 1);

        for(int i = 1; i < kVals.size() - 1; ++i)
            kVals[i] = std::exp(std::log(kMin) + i * deltaLogK);
    
        const double as = 2.1955e-9;
        const double ns = 0.9655;
        const double pivot = 0.5;
        for(int i = 0; i < amplitudes.size(); ++i)
            amplitudes[i] = as * pow(kVals[i]/pivot, ns - 1.0);

        int nMassive = 1;
        double sumMNu = 0.5;
        SplineWithDegenerateNeutrinosParams params(isLinear, 0.022, 0.12, 0.7, 0.1, kVals, amplitudes, 3.046, nMassive, sumMNu, varyNEff, varySumMNu);

//#ifdef COSMO_PLANCK_15
//        PlanckLikelihood planckLike(true, true, true, false, true, false, false, false, 100);
//#else
//        ERROR NOT IMPLEMENTED;
//#endif
//        planckLike.setModelCosmoParams(&params);
        CombinedLikelihood like(true, usePlanck, useWMAP, useBAO, useLRG);
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
        //PolyChord pc(nPar, planckLike, 300, root.str(), 3 * (4 + (varyNEff ? 1 : 0) + (varySumMNu ? 1 : 0)));
        Math::MetropolisHastings scanner(nPar, like, root.str());

        int paramIndex = 0;

        output_screen("Setting parameters" << std::endl);
        scanner.setParam(paramIndex++, "ombh2", 0.02, 0.025, 0.022, 0.001, 0.0002);
        scanner.setParam(paramIndex++, "omch2", 0.1, 0.14, 0.12, 0.003, 0.0006);
        scanner.setParam(paramIndex++, "h", 0.6, 0.75, 0.67, 0.02, 0.003);
        scanner.setParam(paramIndex++, "tau", 0.02, 0.12, 0.07, 0.02, 0.003);
        if(varyNEff)
            scanner.setParam(paramIndex++, "n_eff", 2.0, 5.0, 3.0, 0.6, 0.006);
        if(varySumMNu)
            scanner.setParam(paramIndex++, "sum_mnu", 0.001, 3.0, 0.5, 0.2, 0.02);

        for(int i = 1; i < kVals.size() - 1; ++i)
        {
            std::stringstream paramName;
            paramName << "k_" << i;
            scanner.setParam(paramIndex++, paramName.str(), std::log(kMin), std::log(kMax), std::log(kVals[i]), 1, 0.1);
        }

        for(int i = 0; i < amplitudes.size(); ++i)
        {
            std::stringstream paramName;
            paramName << "a_" << i;
            scanner.setParam(paramIndex++, paramName.str(), aMin, aMax, std::log(amplitudes[i]*1e10), 0.1, 0.01);
        }
        scanner.setParamGauss(paramIndex++, "A_planck", 1, 0.0025);

        check(paramIndex == nPar, "");

        const unsigned long burnin = 500;
        output_screen("Running scan." << std::endl);
        scanner.run(10000, 10, burnin, Math::MetropolisHastings::GELMAN_RUBIN, 0.01, true);
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
