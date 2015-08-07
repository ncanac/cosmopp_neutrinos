#include <vector>
#include <string>
#include <sstream>

#include <mcmc.hpp>
#include <bao_like.hpp>
#include <markov_chain.hpp>
#include <numerics.hpp>
#include <timer.hpp>

int main(int argc, char *argv[])
{
    using namespace Math;
    

    // Create the likelihood
    BAOLikelihood BAOLike;
    std::string root = "bao_mh";

    // Create the Metropolis-Hastings sampler
    MetropolisHastings mh(6, BAOLike, root);

    // Assign parameter names and ranges
    // mh.setParam(0, "ombh2", 0.021, 0.023, 0.022032);
    // mh.setParam(1, "omch2", 0.11, 0.13, 0.12038);
    // mh.setParam(2, "h", 0.66, 0.68, 0.6704);
    // mh.setParam(3, "tau", 0.09, 0.1, 0.0925);
    // mh.setParam(4, "ns", 0.96, 0.97, 0.9619);
    // mh.setParam(5, "As", 3.096, 3.1, 3.098);

    // mh.setParam(0, "ombh2", 0.005, 0.1, 0.022, 0.0003, 0.00005);
    // mh.setParam(1, "omch2", 0.001, 0.99, 0.12, 0.003, 0.0005);
    // mh.setParam(2, "h", 0.2, 1.0, 0.7, 0.02, 0.002);
    // mh.setParam(3, "tau", 0.01, 0.8, 0.1, 0.01, 0.002);
    // mh.setParam(4, "ns", 0.9, 1.1, 1.0, 0.01, 0.002);
    // mh.setParam(5, "As", 2.7, 4.0, 3.0, 0.1, 0.002);

    double ombh2min = 0.01, ombh2max = 0.04, ombh2start = 0.022;
    double ombh2startwidth = (ombh2max - ombh2min)/300;
    double ombh2sampwidth = (ombh2max - ombh2min)/2000;
    mh.setParam(0, "ombh2", ombh2min, ombh2max, ombh2start,
                ombh2startwidth, ombh2sampwidth);
    double omch2min = 0.01, omch2max = 0.04, omch2start = 0.022;
    double omch2startwidth = (omch2max - omch2min)/300;
    double omch2sampwidth = (omch2max - omch2min)/2000;
    mh.setParam(1, "omch2", omch2min, omch2max, omch2start,
                omch2startwidth, omch2sampwidth);
    double hmin = 0.4, hmax = 1.0, hstart = 0.7;
    double hstartwidth = (hmax - hmin)/50;
    double hsampwidth = (hmax - hmin)/500;
    mh.setParam(2, "h", hmin, hmax, hstart, hstartwidth,
                hsampwidth);
    double taumin = 0.05, taumax = 0.2, taustart = 0.1;
    double taustartwidth = (taumax - taumin)/100;
    double tausampwidth = (taumax - taumin)/500;
    mh.setParam(3, "tau", taumin, taumax, taustart,
                taustartwidth, tausampwidth);
    double nsmin = 0.9, nsmax = 1.1, nsstart = 1.0;
    double nsstartwidth = (nsmax - nsmin)/50;
    double nssampwidth = (nsmax - nsmin)/200;
    mh.setParam(4, "ns", nsmin, nsmax, nsstart, nsstartwidth,
                nssampwidth);
    double Asmin = 2.7, Asmax = 4.0, Asstart = 3.0;
    double Asstartwidth = (Asmax - Asmin)/20;
    double Assampwidth = (Asmax - Asmin)/500;
    mh.setParam(5, "As", Asmin, Asmax, Asstart, Asstartwidth,
                Assampwidth);

    // Create the cosmological params
    const double pivot = 0.05;
    // LambdaCDMParams params(0.022032, 0.12038, 0.6704, 0.0925, 0.9619, std::exp(3.098) / 1e10, pivot);
    LambdaCDMParams params(0.022, 0.12, 0.7, 0.1, 1.0, std::exp(3.0) / 1e10, pivot);

    BAOLike.setModelCosmoParams(&params);

    const unsigned long burnin = 1000;
    const int nChains = mh.run(100000, 100, burnin, MetropolisHastings::GELMAN_RUBIN, 1e-4, true);

    output_screen("nChains: " << nChains << std::endl);
}
