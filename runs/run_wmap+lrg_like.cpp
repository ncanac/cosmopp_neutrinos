#include <string>

#include <macros.hpp>
#include <mcmc.hpp>
#include <combined_like.hpp>
#include <timer.hpp>
#include <neutrino_cosmology_params.hpp>
#include <numerics.hpp>

int main(int argc, char *argv[])
{
    using namespace Math;

    std::string root = "wmap+bao_mh";

    CombinedLikelihood like(false, true, false, true);
    const double pivot = 0.05;

    ///////////////////////////////////////////////
    // LambdaCDM
    MetropolisHastings mh(6, like, root);
    
    //mh.setParam(0, "ombh2", 0.01, 0.04, 0.022, 0.005, 0.0005);
    //mh.setParam(1, "omch2", 0.05, 0.3, 0.12, 0.02, 0.002);
    //mh.setParam(2, "h", 0.4, 1.0, 0.7, 0.05, 0.01);
    //mh.setParam(3, "tau", 0.01, 0.2, 0.08, 0.03, 0.003);
    //mh.setParam(4, "ns", 0.9, 1.1, 1.0, 0.05, 0.01);
    //mh.setParam(5, "As", 2.7, 4.0, 3.0, 0.1, 0.01);

    mh.setParam(0, "ombh2", 0.01, 0.04, 0.022, 0.005, 0.0005);
    mh.setParam(1, "omch2", 0.05, 0.3, 0.12, 0.02, 0.002);
    mh.setParam(2, "h", 0.4, 1.0, 0.7, 0.05, 0.02);
    mh.setParam(3, "tau", 0.01, 0.2, 0.08, 0.01, 0.006);
    mh.setParam(4, "ns", 0.9, 1.1, 1.0, 0.05, 0.01);
    mh.setParam(5, "As", 2.7, 4.0, 3.0, 0.1, 0.01);

    LambdaCDMParams params(0.022, 0.12, 0.7, 0.1, 1.0, std::exp(3.0) / 1e10, pivot);
    ///////////////////////////////////////////////

    like.setModelCosmoParams(&params);

    Timer timer("MCMC WMAP+LRG");

    const unsigned long burnin = 500;
    timer.start();
    const int nChains = mh.run(25000, 10, burnin, MetropolisHastings::GELMAN_RUBIN, 0.01, true);
    const unsigned long time = timer.end();
    output_screen("MCMC WMAP+LRG took " << time / 1000000 << " seconds." << std::endl);
    return 0;
}
