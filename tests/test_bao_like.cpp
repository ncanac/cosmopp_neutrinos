#include <macros.hpp>
#include <bao_like.hpp>
#include <cosmo.hpp>

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
    LambdaCDMParams params(omBH2, omCH2, h, tau, ns, as, pivot);

    output_screen("Creating a cosmology object..." << std::endl);
    // Create cosmo object
    Cosmo cosmo;
    cosmo.preInitialize(3000, false, false, false, 0);

    output_screen("Initializing BAO likelihood..." << std::endl);
    // Create likelihood
    BAOLikelihood baoLike(cosmo);
    
    output_screen("Setting the cosmological parameters for likelihood calculation..." << std::endl);
    // Set the cosmological parameters
    baoLike.setCosmoParams(params);

    output_screen("Calculating the likelihood..." << std::endl);
    // Calculate likelihood
    double lnlike = baoLike.likelihood();

    // Output the likelihood
    output_screen("Likelihood = " << lnlike << std::endl);

    return 0;
}
