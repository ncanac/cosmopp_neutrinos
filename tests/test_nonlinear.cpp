#include <iostream>
#include <fstream>

#include <macros.hpp>
#include <cosmological_params.hpp>
#include <cosmo.hpp>

int main(int argc, char *argv[])
{
    // Choose values of the cosmological parameters
    const double h = 0.6704;
    const double omBH2 = 0.022032;
    const double omCH2 = 0.12038;
    const double tau = 0.0925;
    const double ns = 0.9619;
    const double as = 2.2154e-9;
    const double pivot = 0.05;

    const double r = 1e-10;
    const double nt = 0;

    const double nEff = 3.046; 
    const int nMassive = 1;
    const double sumMNu = 0.5;

    // Create an instance of CosmologicalParams. Can use other examples below (just uncomment the appropriate line and comment this line)
    //LCDMWithTensorParams params(omBH2, omCH2, h, tau, ns, as, pivot, r, nt, pivot); 
    //LCDMWithDegenerateNeutrinosParams params(omBH2, omCH2, h, tau, ns, as, pivot, nEff, nMassive, sumMNu);
    LambdaCDMParams params(omBH2, omCH2, h, tau, ns, as, pivot);

    // Choose lMax and create vectors where the resulting cl values will be written
    int lMax = 3000;

    // Create a CMB instance.
    Cosmo cosmo;

    output_screen("Pre-initializing CLASS..." << std::endl);
    // pre-initialize cosmo
    cosmo.preInitialize(lMax, false, false, false, lMax);
    output_screen("OK" << std::endl);

    output_screen("Initializing CLASS..." << std::endl);
    // initialize cosmo
    cosmo.initialize(params, true, false, false, true);
    output_screen("OK" << std::endl);

    Math::TableFunction<double, double> P_lin;
    Math::TableFunction<double, double> P_nw;

    cosmo.getMatterPs(0.0, &P_lin);
    cosmo.getLRGPs(0.0, &P_nw);

    std::ofstream outPlin("P_lin.txt");
    for(auto const point : P_lin)
        outPlin << point.first << " " << point.second << std::endl;
    outPlin.close();

    std::ofstream outPnw("P_nw.txt");
    for(auto const point : P_nw)
        outPnw << point.first << " " << point.second << std::endl;
    outPnw.close();

    return 0;
}

