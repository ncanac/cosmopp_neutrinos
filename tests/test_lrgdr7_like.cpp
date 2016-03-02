#include <macros.hpp>
#include <combined_like.hpp>
#include <neutrino_cosmology_params.hpp>

int main(int argc, char *argv[])
{
    //bool isLinear = true;
    //if(argv[1][0] == 'c')
    //    isLinear = false;

    //const int nKnots = std::atoi(argv[2]);

    //bool varyNEff = false;
    //bool varySumMNu = false;

    //for(int i = 3; i < argc; ++i)
    //{
    //    if(std::string(argv[i]) == "neff")
    //        varyNEff = true;
    //    if(std::string(argv[i]) == "sum_mnu")
    //        varySumMNu = true;
    //}

    //int nPar = 4 + 2 * (nKnots + 2) - 2;
    //
    //if(isLinear)
    //    output_screen("Using a linear spline primordial power spectrum!" << std::endl);

    //if(varyNEff)
    //{
    //    ++nPar;
    //    output_screen("Varying N_eff!" << std::endl);
    //}
    //
    //if(varySumMNu)
    //{
    //    ++nPar;
    //    output_screen("Varying sum_mnu!" << std::endl);
    //}

    // Choose values of the cosmological parameters
    const double h = 0.702;
    const double omBH2 = 0.02262;
    const double omCH2 = 0.1138;
    const double tau = 0.088;
    const double ns = 0.9655;
    const double as = 2.1955e-9;
    const double pivot = 0.05;
    const int nMassive = 0;
    const double nEff = 3.046;
    const double sumMNu = 0.0;

    //const double kMin = 0.8e-6;
    //const double kMax = 1.2;
    //const double aMin = 2.9;
    //const double aMax = 3.6;
    //std::vector<double> kVals(nKnots + 2);
    //std::vector<double> amplitudes(nKnots + 2);

    //kVals[0] = kMin;
    //kVals.back() = kMax;

    //const double deltaLogK = (std::log(kMax) - std::log(kMin)) / (nKnots + 1);

    //for(int i = 1; i < kVals.size() - 1; ++i)
    //    kVals[i] = std::exp(std::log(kMin) + i * deltaLogK);
    //
    //for(int i = 0; i < amplitudes.size(); ++i)
    //    amplitudes[i] = as * pow(kVals[i]/pivot, ns - 1.0);

    //const double h = ;
    //const double omBH2 = ;
    //const double omCH2 = ;
    //const double tau = ;
    //std::vector<double> kVals {0.8e-6, std::exp(), std:exp(), 1.2};
    //std::vector<double> amplitudes {std::exp()/1e10, std::exp()/1e10, std::exp()/1e10, std::exp()/1e10};
    //const int nMassive = 1;
    //const double nEff = ;
    //const double sumMNu = ;

    // Does not work:

    //const double h = 0.694561;
    //const double omBH2 = 0.0246605;
    //const double omCH2 = 0.115974;
    //const double tau = 0.0650107;
    //std::vector<double> kVals {0.8e-6, std::exp(-12.0654), std:exp(-5.22511), 1.2};
    //std::vector<double> amplitudes {std::exp(3.30985)/1e10, std::exp(3.0059)/1e10, std::exp(3.19711)/1e10, std::exp(2.91141)/1e10};
    //const int nMassive = 1;
    //const double nEff = 3.99409;
    //const double sumMNu = 0.796424;

    //const double h = 0.695958;
    //const double omBH2 = 0.0207389;
    //const double omCH2 = 0.121261;
    //const double tau = 0.0866845;
    //std::vector<double> kVals {0.8e-6, std::exp(-5.38061), std:exp(-1.61162), 1.2};
    //std::vector<double> amplitudes {std::exp(3.99405)/1e10, std::exp(3.1774)/1e10, std::exp(3.13913)/1e10, std::exp(3.64804)/1e10};
    //const int nMassive = 1;
    //const double nEff = 2.40787;
    //const double sumMNu = 0.329839;

    output_screen("Creating instance of cosmological parameters..." << std::endl);
    // Create cosmological params
    LambdaCDMParams params(omBH2, omCH2, h, tau, ns, as, pivot);
    //SplineWithDegenerateNeutrinosParams params(isLinear, omBH2, omCH2, h, tau, kVals, amplitudes, nEff, nMassive, sumMNu, varyNEff, varySumMNu);

    output_screen("Initializing LRG DR7 likelihood..." << std::endl);
    //Cosmo cosmo;
    //cosmo.preInitialize(3500, false, true, false, 0, 100, 1e-6, 1);

    std::string datapath = "/Volumes/Data1/ncanac/cosmopp_neutrinos";

    // Create likelihood
    CombinedLikelihood lrgLike(datapath, false, false, false, true, false);
    //LRGDR7Likelihood lrgLike(datapath, cosmo);
    
    output_screen("Setting the cosmological parameters for likelihood calculation..." << std::endl);
    // Set the cosmological parameters
    lrgLike.setCosmoParams(params);

    output_screen("Calculating the likelihood..." << std::endl);
    // Calculate likelihood
    double lnlike = lrgLike.likelihood();

    // Output the likelihood
    output_screen("Likelihood = " << lnlike/2.0 << std::endl);

    return 0;
}
