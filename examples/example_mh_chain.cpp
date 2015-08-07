#include <string>
#include <fstream>

#include <numerics.hpp>
#include <macros.hpp>
#include <markov_chain.hpp>

int main(int argc, char *argv[])
{
    std::string root = "mh";
    int nChains = 4;
    int burnin = 1000;
    // Read the resulting chain(s) with thinning
    const unsigned int thin = 1;
    MarkovChain chain(nChains, root.c_str(), burnin, thin);
    
    // Get the one dimensional marginalized posterior distributions, Gaussian smoothed with a scale of 0.3
    Posterior1D* px = chain.posterior(0, Posterior1D::GAUSSIAN_SMOOTHING);
    Posterior1D* py = chain.posterior(1, Posterior1D::GAUSSIAN_SMOOTHING);
    
    // Get the two dimensional posterior distribution, gaussian smoothed with a scale of 0.25
    Posterior2D* pxy = chain.posterior(0, 1);
    
    // Write the distributions into text files
    
    output_screen("Writing the distributions into text files..." << std::endl);
    px->writeIntoFile("mh_px.txt");
    py->writeIntoFile("mh_py.txt");
    pxy->writeIntoFile("mh_pxy.txt");
    output_screen("OK" << std::endl);
    
    // Write the contour levels for the 2D distribution into a text file. This can be used later to make contour plots
    std::ofstream out("mh_contour_levels.txt");

    out << pxy->get1SigmaLevel() << std::endl;
    //out << pxy->get2SigmaLevel() << std::endl;
    out.close();
    
    // Delete the posterior distributions
    delete px;
    delete py;
    delete pxy;
    return 0;
} 
