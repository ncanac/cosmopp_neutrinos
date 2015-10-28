#include <fstream>
#include <vector>
#include <sstream>

#include <macros.hpp>
#include <neutrino_cosmology_params.hpp>
#include <combined_like.hpp>

int main(int argc, char *argv[])
{
    // Choose values of the cosmological parameters
    const double omBH2 = 0.022032;
    const double omCH2 = 0.12038;
    const double h = 0.6704;
    const double tau = 0.0925;
    const double ns = 0.9619;
    const double as = 2.2154e-9;
    const double pivot = 0.05;

    const double kMin = 0.8e-6;
    const double kMax = 1.2;

    int nKnots = 2;
    std::vector<double> kVals(nKnots + 2);
    std::vector<double> amplitudes(nKnots + 2);

    kVals[0] = kMin;
    kVals.back() = kMax;

    const double deltaLogK = (std::log(kMax) - std::log(kMin)) / (nKnots + 1);
    for(int i = 1; i < kVals.size() - 1; ++i)
        kVals[i] = std::exp(std::log(kMin) + i * deltaLogK);
    
    for(int i = 0; i < amplitudes.size(); ++i)
        amplitudes[i] = as * pow(std::exp(kVals[i]/pivot), ns - 1);
    
    // Create cosmological params
    LambdaCDMParams params(omBH2, omCH2, h, tau, ns, as, pivot);
    bool isLinear = true, varyNEff = true, varySumMNu = true;
    int nMassive = 1;
    double nEff = 3.046, sumMNu = 0.5;
    // SplineWithDegenerateNeutrinosParams params(isLinear, omBH2, omCH2, h, tau, kVals, amplitudes, nEff, nMassive, sumMNu, varyNEff, varySumMNu);


    // Create likelihoods
    CombinedLikelihood like(true, false, false, false);
    
    // Set the cosmological parameters
    like.setCosmoParams(params);

    // Calculate likelihoods
    double lnlike = like.likelihood();

    // Output the likelihoods
    output_screen("Combined likelihood = " << lnlike << std::endl);

    return 0;
}
