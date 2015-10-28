#include <iostream>
#include <fstream>

#include <macros.hpp>
#include <neutrino_cosmology_params.hpp>
#include <cosmo.hpp>

int main(int argc, char *argv[])
{
    // Set cosmological parameters
    // WMAP5 recommended LCDM values:
    const double h = 0.702;
    const double omBH2 = 0.02262;
    const double omCH2 = 0.1138;
    const double tau = 0.088;
    const double ns = 0.962;
    const double as = 2.2154e-9;
    const double pivot = 0.05;

    const double r = 1e-10;
    const double nt = 0;

    const double nEff = 3.046; 
    const int nMassive = 1;
    const double sumMNu = 0.5;

    const int nKnots = 2;
    const double kMin = 0.5e-6;
    const double kMax = 1.5;
    std::vector<double> kVals(nKnots + 2);
    std::vector<double> amplitudes(nKnots + 2, 2e-9);

    kVals[0] = kMin;
    kVals.back() = kMax;

    const double deltaLogK = (std::log(kMax) - std::log(kMin)) / (nKnots + 1);

    for(int i = 1; i < kVals.size() - 1; ++i)
        kVals[i] = std::exp(std::log(kMin) + i * deltaLogK);

    for(int i = 0; i < amplitudes.size(); ++i)
    {
        amplitudes[i] = as * pow(kVals[i]/pivot, ns - 1.0);
        output_screen(kVals[i] << " " << amplitudes[i] << std::endl);
    }

    // Create an instance of CosmologicalParams
    //LambdaCDMParams params(omBH2, omCH2, h, tau, ns, as, pivot);
    //LCDMWithDegenerateNeutrinosParams params(omBH2, omCH2, h, tau, ns, as, pivot, nEff, nMassive, sumMNu);
    SplineWithDegenerateNeutrinosParams params(true, omBH2, omCH2, h, tau, kVals, amplitudes, nEff, nMassive, sumMNu, false, false);

    int lMax = 3000;

    // Create a Cosmo instance.
    Cosmo cosmo;

    output_screen("Pre-initializing CLASS..." << std::endl);
    // pre-initialize cosmo
    cosmo.preInitialize(lMax, false, true, false, 0, 100, 0.8e-6, 1.2);
    output_screen("OK" << std::endl);

    output_screen("Initializing CLASS..." << std::endl);
    // initialize cosmo
    cosmo.initialize(params, true, false, false, true, 0.5);
    output_screen("OK" << std::endl);

    Math::TableFunction<double, double> Ps;
    Math::TableFunction<double, double> P_lrg;

    cosmo.getMatterPs(0.0, &Ps);
    cosmo.getLRGHaloPs(&P_lrg);

    std::ofstream outPs("Ps.txt");
    for(auto const point : Ps)
        outPs << point.first << " " << point.second << std::endl;
    outPs.close();

    std::ofstream outPlrg("P_lrg.txt");
    for(auto const point : P_lrg)
        outPlrg << point.first << " " << point.second << std::endl;
    outPlrg.close();

    return 0;
}

