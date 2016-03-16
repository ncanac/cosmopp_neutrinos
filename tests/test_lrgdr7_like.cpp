#include <macros.hpp>
#include <lrgdr7_like.hpp>
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
    //const double h = 0.702;
    //const double omBH2 = 0.02262;
    //const double omCH2 = 0.1138;
    //const double tau = 0.088;
    //const double ns = 0.98;//9655;
    //const double as = 2.1955e-9;
    //const double pivot = 0.05;
    //const int nMassive = 0;
    //const double nEff = 3.046;
    //const double sumMNu = 0.0;

    // BR09 sample parameters
    //const double h = 0.690167;
    //const double omBH2 = 0.0225740;
    ////const double omBH2 = 0.0280;
    //const double omCH2 = 0.116197;
    //const double tau = 0.090;
    //const double ns = 0.959959;
    //const double as = 2.15547e-9;
    //const double pivot = 0.05;
    //const int nMassive = 0;
    //const double nEff = 3.04;
    //const double sumMNu = 0.0;

    // Planck 2015 parameters
    const double h = 0.6727;
    const double omBH2 = 0.02225;
    //const double omCH2 = 0.1198;
    const double omCH2 = 0.11443;
    const double tau = 0.079;
    const double ns = 0.9645;
    const double as = 2.2065e-9;
    const double pivot = 0.05;
    const int nMassive = 1;
    const double nEff = 3.046;
    const double sumMNu = 0.5;

    // BR09 fiducial model parameters
    //const double h = 0.701;
    //const double omBH2 = 0.0227027;
    //const double omCH2 = 0.114496;
    //const double tau = 0.09;
    //const double ns = 0.96;
    //const double as = 21.76059e-10;
    ////const double as = 20e-10;
    //const double pivot = 0.05;
    //const int nMassive = 0;
    //const double nEff = 3.04;
    //const double sumMNu = 0.0;

    // Model with OmegaM = 0.25 used for computing a_scl
    //const double h = 0.690167;
    //const double omBH2 = 0.01937127431;
    //const double omCH2 = 0.09971134765;
    //const double tau = 0.09;
    //const double ns = 0.959959;
    //const double as = 2.15547e-9;
    //const double pivot = 0.05;
    //const int nMassive = 0;
    //const double nEff = 3.04;
    //const double sumMNu = 0.0;

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
    std::vector<double> kVals {1.0e-6, 10};
    std::vector<double> amplitudes {3.23979378e-9, 1.82817332e-9};
    //const int nMassive = 1;
    //const double nEff = ;
    //const double sumMNu = ;

    // Does not work:

    //output_screen("Creating instance of cosmological parameters..." << std::endl);
    // Create cosmological params
    //LambdaCDMParams params(omBH2, omCH2, h, tau, ns, as, pivot);
    SplineWithDegenerateNeutrinosParams params(true, omBH2, omCH2, h, tau, kVals, amplitudes, nEff, nMassive, sumMNu, false, false);

    //output_screen("Initializing LRG DR7 likelihood..." << std::endl);
    //Cosmo cosmo;
    //cosmo.preInitialize(3500, false, true, false, 0, 100, 1e-6, 1);

    Cosmo cosmo;
    cosmo.preInitialize(3500, false, true, false);
    std::string datapath = "/Volumes/Data1/ncanac/cosmopp_neutrinos";

    // Create likelihood
    //CombinedLikelihood lrgLike(datapath, false, false, false, true, false);
    LRGDR7Likelihood lrgLike(datapath, cosmo);
    
    //output_screen("Setting the cosmological parameters for likelihood calculation..." << std::endl);
    // Set the cosmological parameters
    lrgLike.setCosmoParams(params);

    //output_screen("Calculating the likelihood..." << std::endl);
    // Calculate likelihood
    double lnlike = lrgLike.likelihood();

    // Output the likelihood
    output_screen("Likelihood = " << lnlike/2.0 << std::endl);

    return 0;
}
