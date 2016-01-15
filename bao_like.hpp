#pragma once
#ifdef COSMO_MPI
#include <mpi.h>
#endif

#include <string>
#include <fstream>

#include <macros.hpp>
#include <cosmo.hpp>
#include <cosmo_likelihood.hpp>

class BAOLikelihood : public Math::CosmoLikelihood
{
public:
    BAOLikelihood(Cosmo& cosmo)
    {
        cosmo_ = &cosmo;
    }

    ~BAOLikelihood() {}

    void setCosmoParams(const CosmologicalParams& params)
    {
        params_ = &params;
        cosmo_->initialize(params, true, false, false, false);
    }

    double bao_boss_likelihood()
    {
        // Data at https://github.com/baudren/montepython_public/blob/2.1/data/bao_2014.txt
        double chi2 = 0;
        double da, dr, dv, rsdrag, theo;
        double z, value, error;

        // For each data set, compute angular distance da, radial distance dr,
        // volume distance dv, sound horizon at baryon drag rsdrag, theoretical
        // prediction, and chi2 contribution.

        // 6DF
        z = 0.106;
        value = 0.327; // rs/D_V
        error = 0.015;

        da = cosmo_->getAngularDistance(z);
        dr = z / cosmo_->getHubble(z);
        dv = pow(da*da*(1 + z)*(1 + z)*dr,1.0/3.0);
        rsdrag = cosmo_->getrsdrag();
        theo = rsdrag/dv;
        chi2 = chi2 + pow((theo - value)/error, 2.0);

        // BOSS LOWZ DR10&11 Anderson et al. 1312.4877
        z = 0.32;
        value = 8.47; // D_V/rs
        error = 0.17;

        da = cosmo_->getAngularDistance(z);
        dr = z / cosmo_->getHubble(z);
        dv = pow(da*da*(1 + z)*(1 + z)*dr,1.0/3.0);
        rsdrag = cosmo_->getrsdrag();
        theo = dv/rsdrag;
        chi2 = chi2 + pow((theo - value)/error, 2.0);

        // BOSS CMASS DR10&11 Anderson et al. 1312.4877
        z = 0.57;
        value = 13.77; // D_V/rs
        error = 0.13;

        da = cosmo_->getAngularDistance(z);
        dr = z / cosmo_->getHubble(z);
        dv = pow(da*da*(1 + z)*(1 + z)*dr,1.0/3.0);
        rsdrag = cosmo_->getrsdrag();
        theo = dv/rsdrag;
        chi2 = chi2 + pow((theo - value)/error, 2.0);

        // MGS SDSS DR7 MGS, Ross et al. 1409.3242v1
        z = 0.15;
        value = 4.47; // D_V/rs
        error = 0.16;

        da = cosmo_->getAngularDistance(z);
        dr = z / cosmo_->getHubble(z);
        dv = pow(da*da*(1 + z)*(1 + z)*dr,1.0/3.0);
        rsdrag = cosmo_->getrsdrag();
        theo = dv/rsdrag;
        chi2 = chi2 + pow((theo - value)/error, 2.0);

        return chi2;
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

    double SDSS_likelihood()
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

        double OmK, OmLambda, OmM, h, H0, w, rsdrag;
        double dv1theory, dv2theory, hz1, hz2, DAz1, DAz2;
        double rstodvz1theorydelta, rstodvz2theorydelta, LnLike;

        h = params_->getH();
        H0 = 100.0*h;
        OmK = params_->getOmK();
        OmLambda = params_->getOmLambda();
        OmM = params_->getOmM();
        w = -1.0;
        rsdrag = cosmo_->getrsdrag();
        DAz1 = cosmo_->getAngularDistance(z1);
        DAz2 = cosmo_->getAngularDistance(z2);

        hz1 = sqrt( OmM*pow(1.0+z1,3.0) + OmK*pow(1.0+z1,2.0) + OmLambda*pow(1.0+z1,3.0*(1.0+w)) );
        hz2 = sqrt( OmM*pow(1.0+z2,3.0) + OmK*pow(1.0+z2,2.0) + OmLambda*pow(1.0+z2,3.0*(1.0+w)) );
        dv1theory = pow(((1.0+z1)*DAz1),2.0)*c*z1/H0/hz1/1000.0;
        dv2theory = pow(((1.0+z2)*DAz2),2.0)*c*z2/H0/hz2/1000.0;
        dv1theory = pow(dv1theory,(1.0/3.0));
        dv2theory = pow(dv2theory,(1.0/3.0));

        rstodvz1theorydelta = rsdrag/dv1theory - rstodvz1;
        rstodvz2theorydelta = rsdrag/dv2theory - rstodvz2;

        LnLike = 0.5*((rstodvz1theorydelta) * invcov[0][0] * (rstodvz1theorydelta)
                + 2.0 * (rstodvz1theorydelta) * invcov[0][1] * (rstodvz2theorydelta)
                + (rstodvz2theorydelta) * invcov[1][1] * (rstodvz2theorydelta));
    
        return LnLike;
    }

    double likelihood()
    {
        return bao_boss_likelihood();
    }

    void setModelCosmoParams(CosmologicalParams *params)
    {
        // Sets modelParams_ to initial cosmological parameters.
        modelParams_ = params;
        // Sets vModel_ to a vector containing the parameters in modelParams_.
        modelParams_->getAllParameters(vModel_);
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

private:
    Cosmo* cosmo_;

    const CosmologicalParams* params_; // Standard cosmological parameters
    
    CosmologicalParams* modelParams_; // Cosmological parameters in test models
    std::vector<double> vModel_; // modelParams_ as a vector
};
