#include <vector>
#include <string>
#include <sstream>

#include <mcmc.hpp>
#include <combined_like.hpp>
#include <markov_chain.hpp>
#include <numerics.hpp>
#include <timer.hpp>
#include <neutrino_cosmology_params.hpp>

int main(int argc, char *argv[])
{
    using namespace Math;

    std::string root = "planck+bao_mh";

    CombinedLikelihood combinedLike(true, true);
    const double pivot = 0.05;

    ///////////////////////////////////////////////
    // LambdaCDM
    //MetropolisHastings mh(7, combinedLike, root);

    //mh.setParam(0, "ombh2", 0.01, 0.04, 0.022, 0.002, 0.0002);
    //mh.setParam(1, "omch2", 0.05, 0.3, 0.12, 0.004, 0.0004);
    //mh.setParam(2, "h", 0.4, 1.0, 0.7, 0.02, 0.002);
    //mh.setParam(3, "tau", 0.05, 0.2, 0.1, 0.01, 0.001);
    //mh.setParam(4, "ns", 0.9, 1.1, 1.0, 0.02, 0.002);
    //mh.setParam(5, "As", 2.7, 4.0, 3.0, 0.02, 0.002);

    //mh.setParamGauss(6, "A_planck", 1.0, 0.0025, 1.0, 0.01, 0.001);

    //LambdaCDMParams params(0.022, 0.12, 0.7, 0.1, 1.0, std::exp(3.0) / 1e10, pivot);
    ///////////////////////////////////////////////

    ///////////////////////////////////////////////
    // Uncomment this block for LCDM+neutrinos
    MetropolisHastings mh(9, combinedLike, root);

    mh.setParam(0, "ombh2", 0.01, 0.04, 0.022, 0.002, 0.0002);
    mh.setParam(1, "omch2", 0.05, 0.3, 0.12, 0.004, 0.0004);
    mh.setParam(2, "h", 0.4, 1.0, 0.7, 0.02, 0.002);
    mh.setParam(3, "tau", 0.05, 0.2, 0.1, 0.01, 0.001);
    mh.setParam(4, "ns", 0.9, 1.1, 1.0, 0.02, 0.002);
    mh.setParam(5, "As", 2.7, 4.0, 3.0, 0.02, 0.002);
    mh.setParam(6, "nEff", 1.5, 5.0, 3.0, 0.1, 0.01);
    mh.setParam(7, "sumMNu", 0.01, 2.0, 0.5, 0.2, 0.02);

    mh.setParamGauss(8, "A_planck", 1.0, 0.0025, 1.0, 0.01, 0.001);

    DegenerateNeutrinosParams params(0.022, 0.12, 0.7, 0.1, 1.0, std::exp(3.0) / 1e10, pivot, 3.046, 1, 0.5);
    ///////////////////////////////////////////////

    combinedLike.setModelCosmoParams(&params);

    Timer timer("MCMC PLANCK");

    const unsigned long burnin = 500;
    timer.start();
    const int nChains = mh.run(25000, 10, burnin, MetropolisHastings::GELMAN_RUBIN, 0.01, true);
    const unsigned long time = timer.end();
    output_screen("MCMC Planck took " << time / 1000000 << " seconds." << std::endl);
    return 0;
}
