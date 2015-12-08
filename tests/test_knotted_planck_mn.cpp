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

// argv[1] = c* for cubic, anything else for linear
// argv[2] = number of knots
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

        int nPar = 5 + 2 * (nKnots + 2) - 2;

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
        {
            amplitudes[i] = as * pow(kVals[i]/pivot, ns - 1.0);
        }

        SplineParams params(isLinear, 0.022, 0.12, 0.7, 0.1, kVals, amplitudes);

        PlanckLikelihood like(true, true, true, false, true, false, false, false, 100);
        like.setModelCosmoParams(&params);

        MnScanner scanner(nPar, like, 300, std::string("knotted_planck_mn_"));

        int paramIndex = 0;

        output_screen("Setting parameters" << std::endl);
        scanner.setParam(paramIndex++, "ombh2", 0.02, 0.025);
        scanner.setParam(paramIndex++, "omch2", 0.1, 0.14);
        scanner.setParam(paramIndex++, "h", 0.55, 0.80);
        scanner.setParam(paramIndex++, "tau", 0.04, 0.12);

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
