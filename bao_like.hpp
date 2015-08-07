#ifdef COSMO_MPI
#include <mpi.h>
#endif

#include <string>
#include <fstream>

#include <macros.hpp>
#include <cosmo.hpp>
#include <likelihood_function.hpp>

class BAOLikelihood : public Math::LikelihoodFunction
{
public:
    BAOLikelihood()
    {
        cosmo_ = new Cosmo;
        lMax = 2;
        cosmo_->preInitialize(lMax, false, false, false, lMax);
    }

    ~BAOLikelihood()
    {
        delete cosmo_;
    }

    double calculate(double* params, int nPar)
    {
        // Check to see that modelParams_ is set. This is set by calling setModelCosmoParams().
        check(modelParams_, "model params must be set before calling this function");
        // Check to see that vModel_ is not empty
        check(!vModel_.empty(), "");
        const int nModel = vModel_.size();
        
        // Check that number of parameters passed in is same as number of parameters in nModel
        check(nPar == nModel, "");
    
        // Set all the parameters in vModel_ to the values in params
        for(int i = 0; i < nModel; ++i)
            vModel_[i] = params[i];
    
        // Sets the parameters in modelParams_ to the values in vModel_
        // NEED TO DO: Write setAllParameters for my cosmological parameters class
        modelParams_->setAllParameters(vModel_);
        // Set the cosmological parameters to modelParams_.
        setCosmoParams(*modelParams_);
    
        return likelihood();
    }

    void setCosmoParams(const CosmologicalParams& params)
    {
        params_ = &params;
        cosmo_->initialize(params, true, false, false, false);
    }

    double bao_boss_likelihood()
    {
        // Data at https://github.com/baudren/montepython_public/blob/2.1/data/bao_2014.txt
        return 0;
    }

    double bao_boss_aniso_likelihood()
    {
        // Refer to https://github.com/baudren/montepython_public/tree/2.1/montepython/likelihoods/bao_boss_aniso
        return 0;
    }

    double bao_boss_aniso_gauss_approx_likelihood()
    {
        // Refer to https://github.com/baudren/montepython_public/tree/2.1/montepython/likelihoods/bao_boss_aniso_gauss_approx
        return 0;
    }

    double likelihood()
    {
        double c = 2.99792458e8;
        // From data set in http://arxiv.org/abs/0907.1660
        double z1 = 0.2, z2 = 0.35;
        double rstodvz1 = 0.190533, rstodvz2 = 0.109715;
        double invcov[2][2];
        invcov[0][0] = 30124.1;
        invcov[0][1] = -17226.0;
        invcov[1][0] = invcov[0][1];
        invcov[1][1] = 86976.6;

        double OmK, OmLambda, OmM, h, H0, w, rsdrag, rsdrag_fudge;
        double dv1theory, dv2theory, hz1, hz2, DAz1, DAz2;
        double rstodvz1theorydelta, rstodvz2theorydelta, LnLike;

        h = params_->getH();
        H0 = 100.0*h;
        OmK = params_->getOmK();
        OmLambda = params_->getOmLambda();
        OmM = params_->getOmM();
        w = -1.0;
        rsdrag = cosmo_->getrsdrag();
        rsdrag_fudge = cosmo_->getrsdrag_fudge();
        DAz1 = cosmo_->getAngularDistance(z1);
        DAz2 = cosmo_->getAngularDistance(z2);

        // output_screen("H0: " << H0 << std::endl);
        // output_screen("OmK: " << OmK << std::endl);
        // output_screen("OmLambda: " << OmLambda << std::endl);
        // output_screen("OmM: " << OmM << std::endl);
        // output_screen("rsdrag: " << rsdrag << std::endl);
        // output_screen("rsdrag_fudge: " << rsdrag_fudge << std::endl);
        // output_screen("DAz1: " << DAz1 << std::endl);
        // output_screen("DAz2: " << DAz2 << std::endl);

        hz1 = sqrt( OmM*pow(1.0+z1,3.0) + OmK*pow(1.0+z1,2.0) + OmLambda*pow(1.0+z1,3.0*(1.0+w)) );
        hz2 = sqrt( OmM*pow(1.0+z2,3.0) + OmK*pow(1.0+z2,2.0) + OmLambda*pow(1.0+z2,3.0*(1.0+w)) );
        dv1theory = pow(((1.0+z1)*DAz1),2.0)*c*z1/H0/hz1/1000.0;
        dv2theory = pow(((1.0+z2)*DAz2),2.0)*c*z2/H0/hz2/1000.0;
        dv1theory = pow(dv1theory,(1.0/3.0));
        dv2theory = pow(dv2theory,(1.0/3.0));

        // output_screen("hz1: " << hz1 << std::endl);
        // output_screen("hz2: " << hz2 << std::endl);
        // output_screen("dv1theory: " << dv1theory << std::endl);
        // output_screen("dv2theory: " << dv2theory << std::endl);
    
        rstodvz1theorydelta = rsdrag_fudge/dv1theory - rstodvz1;
        rstodvz2theorydelta = rsdrag_fudge/dv2theory - rstodvz2;

        // output_screen("rstodvz1theorydelta: " << rstodvz1theorydelta << std::endl);
        // output_screen("rstodvz2theorydetla: " << rstodvz2theorydelta << std::endl);
    
        LnLike = 0.5*((rstodvz1theorydelta) * invcov[0][0] * (rstodvz1theorydelta)
                + 2.0 * (rstodvz1theorydelta) * invcov[0][1] * (rstodvz2theorydelta)
                + (rstodvz2theorydelta) * invcov[1][1] * (rstodvz2theorydelta));
    
        return LnLike;
    }

    void setModelCosmoParams(CosmologicalParams *params)
    {
        // Sets modelParams_ to initial cosmological parameters.
        modelParams_ = params;
        // Sets vModel_ to a vector containing the parameters in modelParams_.
        modelParams_->getAllParameters(vModel_);
    }

private:
    int lMax;
    Cosmo* cosmo_;

    const CosmologicalParams* params_; // Standard cosmological parameters
    
    CosmologicalParams* modelParams_; // Cosmological parameters in test models
    std::vector<double> vModel_; // modelParams_ as a vector
};
