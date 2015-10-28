#include <string>

#include <macros.hpp>
#include <mn_scanner.hpp>
#include <combined_like.hpp>
#include <timer.hpp>
#include <neutrino_cosmology_params.hpp>
#include <numerics.hpp>

int main(int argc, char *argv[])
{
    using namespace Math;

    std::string root = "planck+lrg_mn";

    CombinedLikelihood like(true, false, true, false);
    const double pivot = 0.05;

    ///////////////////////////////////////////////
    // LambdaCDM
    MnScanner scanner(7, like, 300, root);
    
    //scanner.setParam(0, "ombh2", 0.01, 0.04, 0.022, 0.005, 0.0005);
    //scanner.setParam(1, "omch2", 0.05, 0.3, 0.12, 0.02, 0.002);
    //scanner.setParam(2, "h", 0.4, 1.0, 0.7, 0.05, 0.01);
    //scanner.setParam(3, "tau", 0.01, 0.2, 0.08, 0.03, 0.003);
    //scanner.setParam(4, "ns", 0.9, 1.1, 1.0, 0.05, 0.01);
    //scanner.setParam(5, "As", 2.7, 4.0, 3.0, 0.1, 0.01);

    scanner.setParam(0, "ombh2", 0.02, 0.025);
    scanner.setParam(1, "omch2", 0.1, 0.2);
    scanner.setParam(2, "h", 0.55, 0.85);
    scanner.setParam(3, "tau", 0.02, 0.20);
    scanner.setParam(4, "ns", 0.9, 1.1);
    scanner.setParam(5, "As", 2.7, 4.0);

    scanner.setParamGauss(6, "A_planck", 1, 0.0025);

    LambdaCDMParams params(0.022, 0.12, 0.7, 0.1, 1.0, std::exp(3.0) / 1e10, pivot);
    ///////////////////////////////////////////////

    like.setModelCosmoParams(&params);

    Timer timer("Multinest Planck+LRG");

    timer.start();
    scanner.run();
    const unsigned long time = timer.end();
    output_screen("Multinest Planck+LRG took " << time / 1000000 << " seconds." << std::endl);
    return 0;
}
