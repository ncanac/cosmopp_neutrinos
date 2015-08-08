#ifdef COSMO_MPI
#include <mpi.h>
#endif

#include <string>
#include <fstream>

#include <macros.hpp>
#include <cosmo.hpp>
#include <likelihood_function.hpp>

class CombinedLikelihood : public Math::LikelihoodFunction
{
public:
    CombinedLikelihood()
    {
        cosmo_ = new Cosmo;
        lMax = 2;
        cosmo_->preInitialize(lMax, false, false, false, lMax);
    }

    ~CombinedLikelihood()
    {
        delete cosmo_;
    }

    void setCosmoParams(const CosmologicalParams& params)
    {
        params_ = &params;
        cosmo_->initialize(params, true, false, false, false);
    }

    double likelihood()
    {
        double lnLike = 0;
        if(wantPlanck)
            lnLike = lnLike + planck_.likelihood();
        if(wantBAO)
            lnLike = lnLike + bao_.likelihood();
        return lnLike;
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
    int lMax;
    Cosmo* cosmo_;

    const CosmologicalParams* params_; // Standard cosmological parameters
    
    CosmologicalParams* modelParams_; // Cosmological parameters in test models
    std::vector<double> vModel_; // modelParams_ as a vector

    // Likelihood objects
    PlanckLikelihood planck_;
    BAOLikelihood bao_;
};
