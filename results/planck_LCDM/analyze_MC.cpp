#include <string>
#include <sstream>

#include <markov_chain.hpp>

int main()
{
    int nChains = 10;
    std::string root = "/home/nicolas/Dropbox/cosmopp_neutrinos/results/planck_LCDM/planck_mh_";
    int burnin = 500;
    int thin = 2;
    MarkovChain chain(nChains, root.c_str(), burnin, thin);

    const int nPoints = 1000;
    const int nPar = 7;
    
    std::ofstream outParamLimits(root + "param_limits.txt");
    for(int i = 0; i < nPar; ++i)
    {
        std::stringstream fileName;
        fileName << root << std::str(i) << ".txt";
        Posterior1D* p = chain.posterior(i, Posterior1D::GAUSSIAN_SMOOTHING);

        p->writeIntoFile(fileName.str().c_str(), nPoints);

        const double median = p->median();
        double lower, upper;
        p->get1SigmaTwoSided(lower, upper);
        const double sigma = (upper - lower) / 2.0;

        outParamLimits << std::str(i) << " = " << median << "+/-" << sigma << std::endl;
    }
    outParamLimits.close();
    return 0;
}
