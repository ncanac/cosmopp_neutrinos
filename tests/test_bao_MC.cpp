#include <vector>
#include <string>
#include <sstream>

#include <markov_chain.hpp>
#include <numerics.hpp>
#include <timer.hpp>

int main(int argc, char *argv[])
{
    using namespace Math;
    std::string root = "test_files/bao_mh";
    const unsigned long burnin = 500;
    const unsigned int thin = 2;
    MarkovChain chain(4, root.c_str(), burnin, thin);

    const int nPoints = 1000;

    std::ofstream outParamLimits("test_files/bao_mh_n4_tol0.01_limits.txt");
    for(int i = 0; i < 6; ++i)
    {
        std::stringstream fileName;
        fileName << "test_files/mcmc_bao_" << i << ".txt";
        Posterior1D* p = chain.posterior(i, Posterior1D::GAUSSIAN_SMOOTHING);

        p->writeIntoFile(fileName.str().c_str(), nPoints);

        const double median = p->median();
        double lower, upper;
        p->get1SigmaTwoSided(lower, upper);
        const double sigma = (upper - lower) / 2.0;

        outParamLimits << i << " = " << median << "+-" << sigma << std::endl;
    }
    outParamLimits.close();

    return 0;
}
