#include <macros.hpp>
#include <lrgdr7_like.hpp>

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
    
    output_screen("Creating instance of cosmological parameters..." << std::endl);
    // Create cosmological params
    LambdaCDMParams params(omBH2, omCH2, h, tau, ns, as, pivot);

    output_screen("Initializing LRG DR7 likelihood..." << std::endl);
    // Create likelihood
    LRGDR7Likelihood lrgLike;
    
    output_screen("Setting the cosmological parameters for likelihood calculation..." << std::endl);
    // Set the cosmological parameters
    lrgLike.setCosmoParams(params);

    output_screen("Calculating the likelihood..." << std::endl);
    // Calculate likelihood
    double lnlike = lrgLike.likelihood();

    // Output the likelihood
    output_screen("Likelihood = " << lnlike << std::endl);

    return 0;
}
