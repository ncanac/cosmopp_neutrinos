#ifdef COSMO_MPI
#include <mpi.h>
#endif

#include <string>
#include <fstream>

#include <macros.hpp>
#include <likelihood_function.hpp>
#include <planck_like.hpp>
#include <bao_like.hpp>

class CombinedLikelihood : public Math::LikelihoodFunction
{
public:
    CombinedLikelihood(bool usePlanck, bool useBAO) : usePlanck_(usePlanck), useBAO_(useBAO)
    {
        if(usePlanck_)
            planckLike_ = new PlanckLikelihood(true, true, true, false, true, false, false, false, 5);
        if(useBAO_)
            BAOLike_ = new BAOLikelihood;
    }

    ~CombinedLikelihood()
    {
        if(usePlanck_)
            delete planckLike_;
        if(useBAO_)
            delete BAOLike_;     
    }

    void setCosmoParams(const CosmologicalParams& params)
    {
        params_ = &params;
        if(usePlanck_)
            planckLike_->setCosmoParams(params);
        if(useBAO_)
            BAOLike_->setCosmoParams(params);
    }

    double likelihood()
    {
        double lnLike = 0; // This is -2*ln(likelihood)
        if(usePlanck_)
            lnLike = lnLike + planckLike_->likelihood();
        if(useBAO_)
            lnLike = lnLike + BAOLike_->likelihood();
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
        int extraPar = 0;
        if(usePlanck_)
            extraPar = 1;
        check(nPar == nModel + extraPar, "wrong number of model params");
    
        // Set all the parameters in vModel_ to the values in params
        for(int i = 0; i < nModel; ++i)
            vModel_[i] = params[i];
    
        // Sets the parameters in modelParams_ to the values in vModel_
        modelParams_->setAllParameters(vModel_);
        // Set the cosmological parameters to modelParams_.
        setCosmoParams(*modelParams_);

        if(usePlanck_)
            planckLike_->setAPlanck(params[nModel]);
    
        return likelihood();
    }

private:
    const CosmologicalParams* params_; // Cosmological parameters for initialization
    
    CosmologicalParams* modelParams_; // Cosmological parameters in test models
    std::vector<double> vModel_; // modelParams_ as a vector

    // Specifies which likelihoods to include
    bool usePlanck_, useBAO_;

    // Likelihood objects
    PlanckLikelihood* planckLike_;
    BAOLikelihood* BAOLike_;
};
