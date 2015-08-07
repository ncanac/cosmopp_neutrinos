#include <string>
#include <fstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <mcmc.hpp>

// Simple two dimensional Gaussian likelihood function
class ExampleMHLikelihood : public Math::LikelihoodFunction
{
public:
    ExampleMHLikelihood() {}
    ~ExampleMHLikelihood() {}

    virtual double calculate(double* params, int nParams)
    {
        check(nParams == 2, "");
        const double x = params[0], y = params[1];

        const double x1 = (x + y) / 2, y1 = (x - y) / 2;
        const double x0 = 0, y0 = 0;
        const double sigmaX = 1, sigmaY = 2;
        const double deltaX = x1 - x0;
        const double deltaY = y1 - y0;

        return deltaX * deltaX / (sigmaX * sigmaX) + deltaY * deltaY / (sigmaY * sigmaY);
    }
};

int main(int argc, char *argv[])
{
    try {
        StandardException exc;

        using namespace Math;

        // Create the likelihood
        ExampleMHLikelihood like;
        std::string root = "mh";

        // Create the Metropolis-Hastings sampler
        MetropolisHastings mh(2, like, root);

        // Assign parameter names and ranges
        const double xMin = -10, xMax = 10, yMin = -10, yMax = 10;
        mh.setParam(0, "x", xMin, xMax, 0, 5, 0.5, 0.05);
        mh.setParam(1, "y", yMin, yMax, 0, 5, 0.5, 0.05);

        // Set the Metropolis-Hastings sampler to vary both parameters together as one block
        std::vector<int> blocks(1, 2);
        mh.specifyParameterBlocks(blocks);

        // Choose the burnin and run
        const unsigned long burnin = 1000;
        const int nChains = mh.run(1000000, 0, burnin, MetropolisHastings::GELMAN_RUBIN, 0.00001);
	output_screen("nChains: " << nChains << std::endl);
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
