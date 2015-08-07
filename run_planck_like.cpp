#include <vector>
#include <string>
#include <sstream>

#include <mcmc.hpp>
#include <planck_like.hpp>
#include <markov_chain.hpp>
#include <numerics.hpp>
#include <timer.hpp>
#include <neutrino_cosmology_params.hpp>

int main(int argc, char *argv[])
{
    using namespace Math;

    std::string root = "planck_mh";

    PlanckLikelihood planckLike(true, true, true, false, true, false, false, false, 5);
    MetropolisHastings mh(9, planckLike, root);

    mh.setParam(0, "ombh2", 0.005, 0.1, 0.022, 0.001, 0.0001);
    mh.setParam(1, "omch2", 0.001, 0.99, 0.12, 0.005, 0.0005);
    mh.setParam(2, "h", 0.2, 1.0, 0.7, 0.02, 0.002);
    mh.setParam(3, "tau", 0.01, 0.8, 0.1, 0.01, 0.001);
    mh.setParam(4, "ns", 0.9, 1.1, 1.0, 0.02, 0.002);
    mh.setParam(5, "As", 2.7, 4.0, 3.0, 0.03, 0.003);
    mh.setParam(6, "nEff", 1.5, 5.0, 3.0, 0.1, 0.01);
    mh.setParam(7, "sumMNu", 0.01, 2.0, 0.5, 0.2, 0.02);

    mh.setParamGauss(8, "A_planck", 1.0, 0.0025, 1.0, 0.01, 0.001);

    const double pivot = 0.05;
    //LambdaCDMParams params(0.022, 0.12, 0.7, 0.1, 1.0, std::exp(3.0) / 1e10, pivot);
    DegenerateNeutrinosParams params(0.022, 0.12, 0.7, 0.1, 1.0, std::exp(3.0) / 1e10, pivot, 3.046, 1, 0.5);
    planckLike.setModelCosmoParams(&params);

    Timer timer("MCMC PLANCK");

    const unsigned long burnin = 500;
    timer.start();
    const int nChains = mh.run(25000, 10, burnin, MetropolisHastings::GELMAN_RUBIN, 0.01, true);
    const unsigned long time = timer.end();
    output_screen("MCMC Planck took " << time / 1000000 << " seconds." << std::endl);
    return 0;
}
