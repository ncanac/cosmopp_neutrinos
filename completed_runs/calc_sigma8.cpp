#include <fstream>
#include <vector>
#include <sstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <neutrino_cosmology_params.hpp>
#include <cosmo.hpp>

// argv[1] = path + file name for root.txt file
// argv[2] = c* for cubic, anything else for linear
// argv[3] = number of knots
// argv[4+] =   neff to vary n effective
//              sum_mnu to vary sum of neutrino masses

int main(int argc, char *argv[])
{
    std::string file_name = argv[1];
    
    bool isLinear = true;
    if(argv[2][0] == 'c')
        isLinear = false;

    const int nKnots = std::atoi(argv[3]);

    bool varyNEff = false;
    bool varySumMNu = false;

    for(int i = 4; i < argc; ++i)
    {
        if(std::string(argv[i]) == "neff")
            varyNEff = true;
        if(std::string(argv[i]) == "sum_mnu")
            varySumMNu = true;
    }

    int nPar = 4 + 2 * (nKnots + 2) - 2;

    if(varyNEff)
        ++nPar;

    if(varySumMNu)
        ++nPar;

    // Starting values for cosmological parameters
    // From Planck 2015 (http://arxiv.org/abs/1502.01589), Table 3, Column 4
    const double h = 0.6727;
    const double omBH2 = 0.02225;
    const double omCH2 = 0.1198;
    const double tau = 0.079;
    const double ns = 0.9645;
    const double as = 3.094; // ln(10^10*as)
    const double pivot = 0.05;

    const double kMin = 1.0e-6;
    const double kMax = 10;
    const double aMin = -2;
    const double aMax = 4;

    std::vector<double> kVals(nKnots + 2);
    std::vector<double> amplitudes(nKnots + 2);

    kVals[0] = kMin;
    kVals.back() = kMax;

    const double deltaLogK = (std::log(kMax) - std::log(kMin)) / (nKnots + 1);

    for(int i = 1; i < kVals.size() - 1; ++i)
        kVals[i] = std::exp(std::log(kMin) + i * deltaLogK);

    for(int i = 0; i < amplitudes.size(); ++i)
        amplitudes[i] = (std::exp(as)/1e10) * pow(kVals[i]/pivot, ns - 1.0);

    int nMassive = (varySumMNu ? 3 : 0);
    double nEff = 3.046; 
    double sumMNu = 0.0;

    SplineWithDegenerateNeutrinosParams params(isLinear, omBH2, omCH2, h, tau, kVals, amplitudes, nEff, nMassive, sumMNu, varyNEff, varySumMNu);

    Cosmo cosmo;
    cosmo.preInitialize(3500, false, true, false, 0, 100, 1e-6, 1.0);

    std::ifstream datafile(file_name);
    std::ofstream outfile(file_name.substr(0, file_name.size()-4) + "sigma8.txt");
    std::string line;
    while(getline(datafile, line))
    {
        std::istringstream iss(line);
        std::vector<double> vec;
        double dummy;
        while(iss >> dummy)
            vec.push_back(dummy);
        std::vector<double> v(&vec[2], &vec[vec.size()-1]);
        params.setAllParameters(v);
        cosmo.initialize(params, true, true, true, true, 1.0);
        for(int i = 0; i < vec.size(); ++i)
            outfile << vec[i] << ' ';
        outfile << cosmo.sigma8() << std::endl;
    }
    outfile.close();
    datafile.close();

    return 0;
}
