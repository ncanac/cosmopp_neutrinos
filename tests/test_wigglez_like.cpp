#include <combined_like.hpp>

#include <cstdlib>

#include <macros.hpp>
#include <cosmological_params.hpp>

int main(int argc, char *argv[])
{
    try {
        double ombh2 = 0.02225;
        double omch2 = 0.1198;
        double h = 0.6727;
        double tau = 0.079;
        double ns = 0.9645;
        double as = 2.196e-9;
        const double pivot = 0.05;

        if(argc >= 7)
        {
            ombh2 = std::atof(argv[1]);
            omch2 = std::atof(argv[2]);
            h = std::atof(argv[3]);
            tau = std::atof(argv[4]);
            ns = std::atof(argv[6]);
            as = std::atof(argv[5]);
        }

        Cosmo cosmo;
        cosmo.preInitialize(3500, false, false, false);
        std::string datapath = "/Volumes/Data1/ncanac/cosmopp_neutrinos";
        WiggleZLikelihood likeA(datapath, cosmo, 'a');
        //WiggleZLikelihood likeB(datapath, cosmo, 'b');
        //WiggleZLikelihood likeC(datapath, cosmo, 'c');
        //WiggleZLikelihood likeD(datapath, cosmo, 'd');
        
        LambdaCDMParams params(ombh2, omch2, h, tau, ns, as, pivot);
        likeA.setCosmoParams(params);
        //likeB.setCosmoParams(params);
        //likeC.setCosmoParams(params);
        //likeD.setCosmoParams(params);

        const double like1 = likeA.likelihood();
        //const double like2 = likeB.likelihood();
        //const double like3 = likeC.likelihood();
        //const double like4 = likeD.likelihood();

        output_screen("Likelihood A = " << like1 << std::endl);
        //output_screen("Likelihood B = " << like2 << std::endl);
        //output_screen("Likelihood C = " << like3 << std::endl);
        //output_screen("Likelihood D = " << like4 << std::endl);
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
