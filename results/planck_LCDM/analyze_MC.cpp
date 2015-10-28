#include <string>
#include <sstream>
#include <fstream>

#include <markov_chain.hpp>

int main()
{
    int nChains = 10;
    std::string root = "/Volumes/Data1/ncanac/cosmopp_neutrinos/results/planck_LCDM/planck_mh";
    int burnin = 500;
    int thin = 2;
    MarkovChain chain(nChains, root.c_str(), burnin, thin);

    const int nPoints = 1000;
    const int nPar = 7;
    
    std::ofstream outParamLimits(root + "_param_limits.txt");
    for(int i = 0; i < nPar; ++i)
    {
        std::stringstream fileName;
        fileName << root << "_" << "param" << std::to_string(i) << ".txt";
        Posterior1D* p = chain.posterior(i, Posterior1D::GAUSSIAN_SMOOTHING);

        p->writeIntoFile(fileName.str().c_str(), nPoints);

        const double median = p->median();
        double lower, upper;
        p->get1SigmaTwoSided(lower, upper);
        const double sigma = (upper - lower) / 2.0;

        outParamLimits << std::to_string(i) << " = " << median << " +/- " << sigma << std::endl;
    }
    outParamLimits.close();
    return 0;
}
