#include <fstream>
#include <vector>
#include <sstream>

#include <macros.hpp>
//#include <planck_like.hpp>
//#include <bao_like.hpp>
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
    
    // Create cosmological params
    LambdaCDMParams params(omBH2, omCH2, h, tau, ns, as, pivot);

    // Create likelihoods
    PlanckLikelihood planckLike(true, true, true, false, true, false, false, false, 5);
    BAOLikelihood BAOLike;
    CombinedLikelihood combLike(true, true);
    
    // Set the cosmological parameters
    planckLike.setCosmoParams(params);
    BAOLike.setCosmoParams(params);
    combLike.setCosmoParams(params);

    // Calculate likelihoods
    double plancklnlike = planckLike.likelihood();
    double baolnlike = BAOLike.likelihood();
    double comblnlike = combLike.likelihood();

    // Output the likelihoods
    output_screen("Planck likelihood = " << plancklnlike << std::endl);
    output_screen("BAO likelihood = " << baolnlike << std::endl);
    output_screen("Combined likelihood = " << comblnlike << std::endl);

    // Check that the individual likelihoods equals the combined likelihood
    check(plancklnlike + baolnlike == comblnlike, "Check failed: likelihoods not equal.");

    return 0;
}

