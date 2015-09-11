#include <string>

#include <macros.hpp>
#include <mcmc.hpp>
#include <lrgdr7_like.hpp>
#include <numerics.hpp>
#include <timer.hpp>

int main(int argc, char *argv[])
{
    using namespace Math;

    // Create the likelihood
    LRGDR7Likelihood LRGLike;
    std::string root = "lrg_mh";

    // Create the Metropolis-Hastings sampler
    MetropolisHastings mh(6, LRGLike, root);

    // Assign parameter names and ranges

    mh.setParam(0, "ombh2", 0.02, 0.025, 0.022, 0.001, 0.0002);
    mh.setParam(1, "omch2", 0.1, 0.14, 0.12, 0.003, 0.0006);
    mh.setParam(2, "h", 0.60, 0.75, 0.67, 0.02, 0.004);
    mh.setParam(3, "tau", 0.02, 0.12, 0.07, 0.02, 0.005);
    mh.setParam(4, "ns", 0.9, 1.04, 0.97, 0.03, 0.006);
    mh.setParam(5, "As", 2.9, 3.3, 3.1, 0.02, 0.004);

    // Create the cosmological params
    const double pivot = 0.05;
    LambdaCDMParams params(0.022, 0.12, 0.67, 0.07, 0.97, std::exp(3.0) / 1e10, pivot);

    LRGLike.setModelCosmoParams(&params);

    Timer timer("MCMC LRG");

    const unsigned long burnin = 500;
    timer.start();
    const int nChains = mh.run(25000, 100, burnin, MetropolisHastings::GELMAN_RUBIN, 0.01, true);
    const unsigned long time = timer.end();
    output_screen("MCMC Planck took " << time / 1000000 << " seconds." << std::endl);
    return 0;
}
