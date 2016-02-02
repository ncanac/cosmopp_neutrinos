#include <fstream>
#include <vector>
#include <sstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <cosmological_params.hpp>
#include <cosmo.hpp>

int main(int argc, char *argv[])
{
    int nPar = 6;
    // Starting values for cosmological parameters
    // From Planck 2015 (http://arxiv.org/abs/1502.01589), Table 3, Column 4
    const double h = 0.6727;
    const double omBH2 = 0.02225;
    const double omCH2 = 0.1198;
    const double tau = 0.079;
    const double ns = 0.9645;
    const double as = 3.094; // ln(10^10*as)
    const double pivot = 0.05;
    LambdaCDMParams params(omBH2, omCH2, h, tau, ns, std::exp(as)/1e10, pivot);

    Cosmo cosmo;
    cosmo.preInitialize(5000, false, false, false, 0, 100, 1e-6, 1.0);

    std::string file_name = "/Volumes/Data1/ncanac/cosmopp_neutrinos/completed_runs/standard_LCDM_mn_build/LCDM_planckposterior.txt";
    std::ifstream datafile(file_name);
    std::string line;
    while(getline(datafile, line))
    {
        std::istringstream iss(line);
        std::vector<double> vec;
        double dummy;
        while(iss >> dummy)
            vec.push_back(dummy);
        std::vector<double> v(&vec[0], &vec[nPar]);
        params.setAllParameters(v);
        cosmo.initialize(params, true, false, false, true);
        for(int i = 0; i < vec.size(); ++i)
            output_screen(vec[i] << " ");
        output_screen(cosmo.sigma8() << std::endl);
    }
    datafile.close();

    return 0;
}
