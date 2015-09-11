#pragma once
#include <likelihood_function.hpp>
#include <cosmological_params.hpp>

namespace Math
{

/// An abstract class for likelihood functions which calculate their likelihood
/// by first setting up a CosmologicalParams object.
class CosmoLikelihood : public LikelihoodFunction
{
public:
    virtual ~CosmoLikelihood() {}

    /// This function sets the cosmological parameters member object of the
    /// likelihood object to be equal to params and also initializes the Class
    /// object. Should be called before calling likelihood().
    /// Example call: setCosmoParams(params)
    /// \param params The cosmological parameters object.
    virtual void setCosmoParams(const CosmologicalParams& params) = 0;

    /// This function calculates the likelihood after setting the cosmological
    /// parameters through setCosmoParams().
    virtual double likelihood() = 0;

    /// This function sets the cosmological model parameters member object of
    /// the likelihood function to be equal to params. It also sets the member
    /// vector of cosmological parameters for use with the calculate function.
    /// Should be called before running one of the parameter sampling
    /// algorithm, i.e. before calling calculate.
    /// Example call: setModelCosmoParams(&params)
    /// \params A pointer to the model cosmological parameters object
    virtual void setModelCosmoParams(CosmologicalParams *params) = 0;
};

} // namespace Math
