#include <fstream>
#include <string>
#include <sstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <likelihood_function.hpp>
#include <mn_scanner.hpp>
#include <planck_like.hpp>

int main(int argc, char *argv[])
{
    try {
        StandardException exc;

        // Create the likelihood
        PlanckLikelihood like(true, true, true, false, true, false, false);

        // Create the Multinest scanner with 7 parameters and 300 live points
        MnScanner scanner(7, like, 300, std::string("planck_multinest_"));

        // Set the parameter names and ranges
        scanner.setParam(0, std::string("ombh2"), 0.02, 0.025);
        scanner.setParam(1, std::string("omch2"), 0.1, 0.2);
        scanner.setParam(2, std::string("h"), 0.55, 0.85);
        scanner.setParam(3, std::string("tau"), 0.01, 0.30);
        scanner.setParam(4, std::string("ns"), 0.9, 1.1);
        scanner.setParam(5, std::string("as"), 2.7, 3.5);

        scanner.setParamGauss(6, "A_planck", 1.0, 0.0025);

        // Set the model for planck likelihood by specifying an example parameter set.
        const double pivot = 0.05;
        LambdaCDMParams par(0.022, 0.12, 0.7, 0.1, 1.0, std::exp(3.0) / 1e10, pivot);
        like.setModelCosmoParams(&par);

        // Run the scanner. The Results will be output in corresponding files at the end
        scanner.run();
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}

