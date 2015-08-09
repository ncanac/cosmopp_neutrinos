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
    CombinedLikelihood(bool usePlanck, bool useBAO)
    {
        if(usePlanck)
            likeVector_.push_back(PlanckLikelihood(true, true, true, true, true, true, true));

        if(useBAO)
            likeVector_.push_back(BAOLikelihood());
    }

    ~CombinedLikelihood() { }

    void setCosmoParams(const CosmologicalParams& params)
    {
        params_ = &params;
        for(int i = 0; i < likeVector_.size(); ++i)
            likeVector_[i].setCosmoParams(params_);
    }

    double likelihood()
    {
        double lnLike = 0; // This is actually -2*ln(likelihood)
        for(int i = 0; i < likeVector_.size(); ++i)
            lnLike = lnLike + likeVector_[i].likelihood;
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

    const CosmologicalParams* params_; // Cosmological parameters for initialization
    
    CosmologicalParams* modelParams_; // Cosmological parameters in test models
    std::vector<double> vModel_; // modelParams_ as a vector

    // Vector of likelihood objects
    std::vector<Math::LikelihoodFunction> likeVector_;
};
