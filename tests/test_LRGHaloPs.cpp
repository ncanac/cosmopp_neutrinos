#include <iostream>
#include <fstream>

#include <macros.hpp>
#include <neutrino_cosmology_params.hpp>
#include <cosmo.hpp>

int main(int argc, char *argv[])
{
    const double zNEAR = 0.235, zMID = 0.342, zFAR = 0.421;
    // Set cosmological parameters
    // WMAP5 recommended LCDM values:
    double h = 0.702;
    double omBH2 = 0.02262;
    double omCH2 = 0.1138;
    double tau = 0.088;
    const double ns = 0.962;
    const double as = 2.2154e-9;
    const double pivot = 0.05;

    const double r = 1e-10;
    const double nt = 0;

    double nEff = 3.046; 
    int nMassive = 1;
    double sumMNu = 0.1;

    const int nKnots = 2;
    const double kMin = 0.8e-6;
    const double kMax = 1.2;
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
        //output_screen(kVals[i] << " " << amplitudes[i] << std::endl);
    }

    // Create an instance of CosmologicalParams
    //LambdaCDMParams params(omBH2, omCH2, h, tau, ns, as, pivot);
    //LCDMWithDegenerateNeutrinosParams params(omBH2, omCH2, h, tau, ns, as, pivot, nEff, nMassive, sumMNu);
    SplineWithDegenerateNeutrinosParams params(true, omBH2, omCH2, h, tau, kVals, amplitudes, nEff, nMassive, sumMNu, true, true);

    int lMax = 5000;

    // Create a Cosmo instance.
    Cosmo cosmo;

    output_screen("Pre-initializing CLASS..." << std::endl);
    // pre-initialize cosmo
    cosmo.preInitialize(lMax, false, true, false, 0, 100, 1.0e-6, 1.0);
    output_screen("OK" << std::endl);

    output_screen("Initializing CLASS..." << std::endl);
    // initialize cosmo
    cosmo.initialize(params, true, false, false, true, 0.5);
    output_screen("OK" << std::endl);

    //Math::TableFunction<double, double> Ps;
    Math::TableFunction<double, double> P_lrg;

    //cosmo.getMatterPs(zNEAR, &Ps);
    //std::ofstream outPs("Ps.txt");
    //for(auto const &point : Ps)
    //    outPs << point.first << " " << point.second << std::endl;
    //outPs.close();

    cosmo.getLRGHaloPs("/Volumes/Data1/ncanac/cosmopp_neutrinos/data/LRGDR7/", &P_lrg);
    std::ofstream outPlrg("P_lrg.txt");
    for(auto const &point : P_lrg)
        outPlrg << point.first << " " << point.second << std::endl;
    outPlrg.close();

    // Choose different values of parameters and output P_lrg

    h = 0.756994;
    omBH2 = 0.0215047;
    omCH2 = 0.119089;
    tau = 0.0420033;
    kVals = {0.8e-6, std::exp(-10.8105), std:exp(-5.56935), 1.2};
    amplitudes = {std::exp(3.23493)/1e10, std::exp(3.03256)/1e10, std::exp(3.04898)/1e10, std::exp(3.09344)/1e10};
    nMassive = 1;
    nEff = 3.06841;
    sumMNu = 0.96237;

    std::vector<double> v {omBH2, omCH2, h, tau, nEff, sumMNu, kVals[1], kVals[2], amplitudes[0], amplitudes[1], amplitudes[2], amplitudes[3]};

    params.setAllParameters(v);
    cosmo.initialize(params, true, false, false, true, 0.5);

    cosmo.getLRGHaloPs("/Volumes/Data1/ncanac/cosmopp_neutrinos/data/LRGDR7/", &P_lrg);
    outPlrg.open("P_lrg2.txt");
    for(auto const &point : P_lrg)
        outPlrg << point.first << " " << point.second << std::endl;
    outPlrg.close();

    return 0;
}

