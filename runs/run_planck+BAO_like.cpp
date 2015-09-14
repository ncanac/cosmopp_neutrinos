#include <string>

#include <macros.hpp>
#include <mcmc.hpp>
#include <combined_like.hpp>
#include <timer.hpp>
#include <neutrino_cosmology_params.hpp>

int main(int argc, char *argv[])
{
    using namespace Math;

    std::string root = "planck+bao_mh";

    CombinedLikelihood combinedLike(true, false, true, false);
    const double pivot = 0.05;
    const double nMassive = 1;

    ///////////////////////////////////////////////////////
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
    ///////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////
    // Uncomment this block for LCDM+neutrinos
    //MetropolisHastings mh(9, combinedLike, root);

    //mh.setParam(0, "ombh2", 0.02, 0.025, 0.022, 0.001, 0.0002);
    //mh.setParam(1, "omch2", 0.1, 0.14, 0.12, 0.004, 0.0004);
    //mh.setParam(2, "h", 0.64, 0.74, 0.673, 0.02, 0.003);
    //mh.setParam(3, "tau", 0.02, 0.16, 0.07, 0.02, 0.003);
    //mh.setParam(4, "ns", 0.9, 1.04, 0.96, 0.02, 0.004);
    //mh.setParam(5, "As", 2.9, 3.3, 3.1, 0.03, 0.003);
    //mh.setParam(6, "nEff", 1.0, 4.0, 2.04, 0.2, 0.02);
    //mh.setParam(7, "sumMNu", 0.001, 2.0, 0.1, 0.2, 0.03);

    //mh.setParamGauss(8, "A_planck", 1.0, 0.0025, 1.0, 0.02, 0.002);

    //mh.setParam(0, "ombh2", 0.02, 0.025, 0.022, 0.001, 0.0002);
    //mh.setParam(1, "omch2", 0.1, 0.14, 0.12, 0.003, 0.0006);
    //mh.setParam(2, "h", 0.60, 0.75, 0.67, 0.02, 0.004);
    //mh.setParam(3, "tau", 0.02, 0.12, 0.07, 0.02, 0.005);
    //mh.setParam(4, "ns", 0.9, 1.04, 0.97, 0.03, 0.006);
    //mh.setParam(5, "As", 2.9, 3.3, 3.1, 0.02, 0.004);
    //mh.setParam(6, "nEff", 2.0, 4.0, 3.04, 0.1, 0.02);
    //mh.setParam(7, "sumMNu", 0.001, 2.0, 0.1, 0.1, 0.03);

    //mh.setParamGauss(8, "A_planck", 1.0, 0.0025, 1.0, 0.02, 0.002);

    //DegenerateNeutrinosParams params(0.022, 0.12, 0.67, 0.07, 0.97, std::exp(3.1) / 1e10, pivot, 3.04, nMassive, 0.1);
    ///////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////
    // LCDM with linear spline PPS
    // For a standard power law PPS of the form As(k/k*)^(ns-1):
    //As ~ 2.2154e-9
    //ns ~ 0.9619
    //fix kmin = 1e-6, kmax = 1
    //Choose ln k between ln kmin and ln kmax
    MetropolisHastings mh(9, combinedLike, root);

    double kmin = std::log(1e-6), kmax = std::log(1.0); // -13.8 < k < 0
    double kstartWidth = (kmax - kmin) / 20.0;
    double ksampWidth = (kmax - kmin) / 200.0;
    // Rough bounds for PPS given standard PPS evaluated at kmin and kmax
    double Amin = 3.0, Amax = 3.7;
    double AstartWidth = (Amax - Amin) / 20.0;
    double AsampWidth = (Amax - Amin)/ 200.0;
    mh.setParam(0, "ombh2", 0.02, 0.025, 0.022, 0.001, 0.0002);
    mh.setParam(1, "omch2", 0.1, 0.14, 0.12, 0.003, 0.0006);
    mh.setParam(2, "h", 0.60, 0.75, 0.67, 0.02, 0.003); // 0.004 -> 0.003
    mh.setParam(3, "tau", 0.02, 0.12, 0.07, 0.02, 0.003); // 0.005 -> 0.003
    mh.setParam(4, "A0", Amin, Amax, 3.60, AstartWidth, AsampWidth*10);
    mh.setParam(5, "k1", kmin, kmax, std::log(1e-3), kstartWidth, ksampWidth*6);
    mh.setParam(6, "A1", Amin, Amax, 3.33, AstartWidth, AsampWidth*4);
    mh.setParam(7, "An", Amin, Amax, 3.07, AstartWidth, AsampWidth*3);

    mh.setParamGauss(8, "A_planck", 1.0, 0.0025, 1.0, 0.02, 0.002);

    std::vector<double> kVals {1e-6, 1e-3, 1};
    std::vector<double> amplitudes {std::exp(3.60) / 1e10,
                                    std::exp(3.33) / 1e10,
                                    std::exp(3.07) / 1e10};
    LCDMwithLinearSplineParams params(0.022, 0.12, 0.67, 0.07, kVals, amplitudes);
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
