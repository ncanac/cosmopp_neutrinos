#include <macros.hpp>
#include <wigglez_like.hpp>
#include <neutrino_cosmology_params.hpp>

int main(int argc, char *argv[])
{
    // Choose values of the cosmological parameters
    const double h = 0.702;
    const double omBH2 = 0.02262;
    const double omCH2 = 0.1138;
    const double tau = 0.088;
    const double ns = 0.962;
    const double as = 2.2154e-9;
    const double pivot = 0.05;
    
    output_screen("Creating instance of cosmological parameters..." << std::endl);
    // Create cosmological params
    // LambdaCDMParams params(omBH2, omCH2, h, tau, ns, as, pivot);
    bool isLinear = true;
    const double kMin = 0.8e-6;
    const double kMax = 1.2;
    std::vector<double> kVals {std::exp(std::log(kMin)), std::exp(std::log(kMax))};
    std::vector<double> amplitudes {as * pow(kVals[0]/pivot, ns - 1.0), as * pow(kVals[1]/pivot, ns - 1.0)};
    const double nEff = 3.046;
    const double nMassive = 0;
    const double sumMNu = 0.0;
    const double varyNEff = false;
    const double varySumMNu = false;
    SplineWithDegenerateNeutrinosParams params(isLinear, omBH2, omCH2, h, tau, kVals, amplitudes, nEff, nMassive, sumMNu, varyNEff, varySumMNu);

    output_screen("Initializing WiggleZ likelihood..." << std::endl);
    // Create likelihood
    bool primordialInit = true;
    WiggleZLikelihood likea(primordialInit, "/Volumes/Data1/ncanac/cosmopp_neutrinos", 'a');
    WiggleZLikelihood likeb(primordialInit, "/Volumes/Data1/ncanac/cosmopp_neutrinos", 'b');
    WiggleZLikelihood likec(primordialInit, "/Volumes/Data1/ncanac/cosmopp_neutrinos", 'c');
    WiggleZLikelihood liked(primordialInit, "/Volumes/Data1/ncanac/cosmopp_neutrinos", 'd');
    
    output_screen("Setting the cosmological parameters for likelihood calculation..." << std::endl);
    // Set the cosmological parameters
    likea.setCosmoParams(params);
    likeb.setCosmoParams(params);
    likec.setCosmoParams(params);
    liked.setCosmoParams(params);

    output_screen("Calculating the likelihood..." << std::endl);
    // Calculate likelihood
    double lnlike = 0.0;
    double lnliketot = 0.0;
    lnlike = likea.likelihood();
    lnliketot += lnlike;
    output_screen("WiggleZ a likelihood = " << lnlike << std::endl);
    lnlike = likeb.likelihood();
    lnliketot += lnlike;
    output_screen("WiggleZ b likelihood = " << lnlike << std::endl);
    lnlike = likec.likelihood();
    lnliketot += lnlike;
    output_screen("WiggleZ c likelihood = " << lnlike << std::endl);
    lnlike = liked.likelihood();
    lnliketot += lnlike;
    output_screen("WiggleZ d likelihood = " << lnlike << std::endl);

    // Output the likelihood
    output_screen("Likelihood = " << lnliketot << std::endl);

    return 0;
}
