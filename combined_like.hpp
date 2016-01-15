#pragma once
#ifdef COSMO_MPI
#include <mpi.h>
#endif

#include <string>
#include <fstream>

#include <macros.hpp>
#include <planck_like.hpp>
#include <wmap9_like.hpp>
#include <bao_like.hpp>
#include <lrgdr7_like.hpp>
#include <numerics.hpp>
#include <cosmo_likelihood.hpp>
#include <wigglez_like.hpp>

class CombinedLikelihood : public Math::CosmoLikelihood
{
public:
    CombinedLikelihood(std::string datapath, Cosmo& cosmo, bool usePlanck, bool useWMAP, bool useBAO, bool useLRG, bool useWiggleZ) : usePlanck_(usePlanck), useWMAP_(useWMAP), useBAO_(useBAO), useLRG_(useLRG), useWiggleZ_(useWiggleZ)
    {
        nLikes_ = 0;
        //check(!(usePlanck_ && useWMAP_), "Both Planck and WMAP likelihoods should not be used at the same time.");
        //check(!(useBAO_ && useLRG_), "Both BAO and LRG likelihoods should not be used at the same time.");
        if(usePlanck_)
            planckLike_ = new PlanckLikelihood(true, true, true, false, true, false, false, false, 100);
        if(useWMAP_)
            wmapLike_ = new WMAP9Likelihood(true, true, true, true, true, true);
        if(useBAO_)
        {
            likes_.push_back(new BAOLikelihood(cosmo));
            ++nLikes_;
        }
        if(useLRG_)
        {
            likes_.push_back(new LRGDR7Likelihood(datapath, cosmo));
            ++nLikes_;
        }
        if(useWiggleZ_)
        {
            likes_.push_back(new WiggleZLikelihood(datapath, cosmo, 'a')); 
            likes_.push_back(new WiggleZLikelihood(datapath, cosmo, 'b')); 
            likes_.push_back(new WiggleZLikelihood(datapath, cosmo, 'c')); 
            likes_.push_back(new WiggleZLikelihood(datapath, cosmo, 'd')); 
            nLikes_ += 4;
        }
    }

    ~CombinedLikelihood()
    {
        if(usePlanck_)
            delete planckLike_;
        if(useWMAP_)
            delete wmapLike_;
    }

    void setCosmoParams(const CosmologicalParams& params)
    {
        params_ = &params;
        if(usePlanck_)
            planckLike_->setCosmoParams(params);
        if(useWMAP_)
        {
            wmapLike_->setCosmoParams(params);
            wmapLike_->calculateCls();
        }
        for(int i = 0; i < nLikes_; ++i)
            likes_[i]->setCosmoParams(params);
    }

    double likelihood()
    {
        double lnLike = 0; // This is -2*ln(likelihood)
        if(usePlanck_)
            lnLike = lnLike + planckLike_->likelihood();
        if(useWMAP_)
            lnLike = lnLike + wmapLike_->likelihood();
        for(int i = 0; i < nLikes_; ++i)
            lnLike += likes_[i]->likelihood();
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
    bool usePlanck_, useWMAP_, useBAO_, useLRG_, useWiggleZ_;

    // Likelihood objects
    PlanckLikelihood* planckLike_;
    WMAP9Likelihood* wmapLike_;
    std::vector<Math::CosmoLikelihood*> likes_;
    //BAOLikelihood* BAOLike_;

    // Number of likelihood objects excluding Planck and WMAP
    int nLikes_;
};
