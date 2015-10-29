#pragma once
#ifdef COSMO_MPI
#include <mpi.h>
#endif

#include <string>
#include <fstream>
#include <cmath>
#include <iostream>
#include <iomanip>

#include <macros.hpp>
#include <cosmo.hpp>
#include <cosmo_likelihood.hpp>
#include <matrix_impl.hpp>
#include <cubic_spline.hpp>

#include <class.h>

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>

class LRGDR7Likelihood : public Math::CosmoLikelihood
{
public:
    LRGDR7Likelihood(bool primordialInitialize) 
    {
        cosmo_ = new Cosmo;
        lMax_ = 3000;
        cosmo_->preInitialize(lMax_, false, primordialInitialize, false, 0);

        // Number of points and kbands in the input files
        num_mpk_points_full_ = 45;
        num_mpk_kbands_full_ = 250;

        // Decide which bandpowers to use, min to max
        min_mpk_points_use_ = 1;
        max_mpk_points_use_ = 45;
        min_mpk_kbands_use_ = 1;
        max_mpk_kbands_use_ = 250;

        k_size_ = max_mpk_kbands_use_ - min_mpk_kbands_use_ + 1;
        k_.resize(k_size_);
        kh_.resize(k_size_);

        // Set parameters associated with nuisance parameters
        k1_ = 0.1;
        k2_ = 0.2;
        s1_ = 0.04;
        s2_ = 0.1;
        a1maxval_ = 1.1482;
        nptsa1_ = 41;
        nptsa2_ = 41;
        nptstot_ = 325;

        // Read in data file containing kbands
        root_ = "/Volumes/Data1/ncanac/cosmopp_neutrinos/data/LRGDR7/";
        std::ifstream datafile(root_ + "data/lrgdr7_kbands.txt");
        for(int i = 0; i < num_mpk_kbands_full_; ++i)
            if((i+2 > min_mpk_kbands_use_) && (i < max_mpk_kbands_use_))
                datafile >> kh_[i-min_mpk_kbands_use_+1];
        datafile.close();

        // Read in window functions
        n_size_ = max_mpk_points_use_ - min_mpk_points_use_ + 1; // 45
        window_.resize(n_size_, k_size_, 0);
        //std::vector< std::vector<double> > window(n_size_, std::vector<double>(k_size_));
        datafile.open(root_ + "data/lrgdr7_windows.txt");
        for(int i = 0; i < num_mpk_points_full_; ++i)
            for(int j = 0; j < k_size_; ++j)
                datafile >> window_(i, j);
        datafile.close();

        zerowindowfxn_.resize(k_size_, 1, 0);
        datafile.open(root_ + "data/lrgdr7_zerowindowfxn.txt");
        for(int i = 0; i < k_size_; ++i)
            datafile >> zerowindowfxn_(i, 0);
        datafile.close();

        zerowindowfxnsubdat_.resize(n_size_, 1, 0);
        datafile.open(root_ + "data/lrgdr7_zerowindowfxnsubtractdat.txt");
        datafile >> zerowindowfxnsubdatnorm_;
        for(int i = 0; i < n_size_; ++i)
            datafile >> zerowindowfxnsubdat_(i, 0);
        datafile.close();

        // Read in measurements
        std::string line;
        P_obs_.resize(n_size_, 1, 0);
        P_err_.resize(n_size_, 1, 0);
        datafile.open(root_ + "data/lrgdr7_ccmeasurements.txt");
        // Skip first two lines
        std::getline(datafile, line);
        std::getline(datafile, line);
        for(int i = 0; i < num_mpk_points_full_; ++i)
        {
            std::getline(datafile, line);
            if((i+2 > min_mpk_points_use_) && (i < max_mpk_points_use_))
            {
                std::istringstream iss(line);
                double kdum, klodum, khidum, pdum, errdum, dum;
                iss >> kdum >> klodum >> khidum >> pdum >> errdum >> dum;
                P_obs_(i-min_mpk_points_use_+1, 0) = pdum;
                P_err_(i-min_mpk_points_use_+1, 0) = errdum;
            }
        }
        datafile.close();

        // Read in inverse covariance matrix
        invcov_.resize(n_size_, n_size_, 0);
        datafile.open(root_ + "data/lrgdr7_invcov.txt");
        for(int i = 0; i < num_mpk_points_full_; ++i)
            if((i+2 > min_mpk_points_use_) && (i < max_mpk_points_use_))
                for(int j = 0; j < num_mpk_points_full_; ++j)
                    if((j+2 > min_mpk_points_use_) && (j < max_mpk_points_use_))
                        datafile >> invcov_(i, j);
        datafile.close();
    }

    ~LRGDR7Likelihood()
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
        modelParams_->setAllParameters(vModel_);
        // Set the cosmological parameters to modelParams_.
        setCosmoParams(*modelParams_);
    
        return likelihood();
    }

    double likelihood()
    {
        Math::Matrix<double> mpk_raw(k_size_, 1, 0);
        Math::Matrix<double> mpk_Pth(k_size_, 1, 0);
        Math::Matrix<double> mpk_Pth_k(k_size_, 1, 0);
        Math::Matrix<double> mpk_Pth_k2(k_size_, 1, 0);
        Math::Matrix<double> mpk_WPth(n_size_, 1, 0);
        Math::Matrix<double> mpk_WPth_k(n_size_, 1, 0);
        Math::Matrix<double> mpk_WPth_k2(n_size_, 1, 0);
        Math::Matrix<double> kh_scaled(k_size_, 1, 0); // h^-1 Mpc^-1

        // Compute scaling factor
        const double zeffDR7 = 0.312782; // redshift at which a_scl is evaluated.
        const double dvfid = 1211.884207; // value of dv for fiducial model
        double da = cosmo_->getAngularDistance(zeffDR7);
        double dr = zeffDR7 / cosmo_->getHubble(zeffDR7);
        double dv = pow(da*da*(1 + zeffDR7)*(1 + zeffDR7)*dr,1.0/3.0);
        double a_scl = dv/dvfid; // Value in BR09 example model = 1.012937

        // Initialize halopowerlrgtheory
        Math::TableFunction<double, double> halopowerlrgtheory;
        cosmo_->getLRGHaloPs(&halopowerlrgtheory);

        // Calculate kh_scaled and mpk_raw, which is just halopowerlrgtheory evaluated at kh_scaled*h
        // mpk_raw is in units of h^3 Mpc^3
        for(int i = 0; i < k_size_; ++i)
        {
            kh_scaled(i, 0) = a_scl * kh_[i];
            mpk_raw(i, 0) = halopowerlrgtheory.evaluate(kh_scaled(i, 0) / pow(a_scl, 3.0));
        }

        // Initialize
        mpk_Pth = mpk_raw;
        for(int i = 0; i < k_size_; ++i)
        {
            mpk_Pth_k(i, 0) = mpk_Pth(i, 0) * kh_scaled(i, 0);
            mpk_Pth_k2(i, 0) = mpk_Pth(i, 0) * pow(kh_scaled(i, 0), 2.0);
        }

        Math::Matrix<double>::multiplyMatrices(window_, mpk_Pth, &mpk_WPth);
        Math::Matrix<double>::multiplyMatrices(window_, mpk_Pth_k, &mpk_WPth_k);
        Math::Matrix<double>::multiplyMatrices(window_, mpk_Pth_k2, &mpk_WPth_k2);

        Math::Matrix<double> tempMat;
        Math::Matrix<double> zerowindowfxn_transpose = zerowindowfxn_.getTranspose();
        Math::Matrix<double>::multiplyMatrices(zerowindowfxn_transpose, mpk_Pth, &tempMat);
        double sumzerow_Pth = tempMat(0, 0) / zerowindowfxnsubdatnorm_;
        Math::Matrix<double>::multiplyMatrices(zerowindowfxn_transpose, mpk_Pth_k, &tempMat);
        double sumzerow_Pth_k = tempMat(0, 0) / zerowindowfxnsubdatnorm_;
        Math::Matrix<double>::multiplyMatrices(zerowindowfxn_transpose, mpk_Pth_k2, &tempMat);
        double sumzerow_Pth_k2 = tempMat(0, 0) / zerowindowfxnsubdatnorm_;

        Math::Matrix<double> covdat(n_size_, 1, 0);
        Math::Matrix<double> covth(n_size_, 1, 0);
        Math::Matrix<double> covth_k(n_size_, 1, 0);
        Math::Matrix<double> covth_k2(n_size_, 1, 0);
        Math::Matrix<double> covth_zerowin(n_size_, 1, 0);
        Math::Matrix<double>::multiplyMatrices(invcov_, P_obs_, &covdat);
        Math::Matrix<double>::multiplyMatrices(invcov_, mpk_WPth, &covth);
        Math::Matrix<double>::multiplyMatrices(invcov_, mpk_WPth_k, &covth_k);
        Math::Matrix<double>::multiplyMatrices(invcov_, mpk_WPth_k2, &covth_k2);
        Math::Matrix<double>::multiplyMatrices(invcov_, zerowindowfxnsubdat_, &covth_zerowin);

        Math::Matrix<double> P_obs_transpose = P_obs_.getTranspose();
        Math::Matrix<double>::multiplyMatrices(P_obs_transpose, covdat, &tempMat);
        double sumDD = tempMat(0, 0);
        Math::Matrix<double>::multiplyMatrices(P_obs_transpose, covth, &tempMat);
        double sumDT = tempMat(0, 0);
        Math::Matrix<double>::multiplyMatrices(P_obs_transpose, covth_k, &tempMat);
        double sumDT_k = tempMat(0, 0);
        Math::Matrix<double>::multiplyMatrices(P_obs_transpose, covth_k2, &tempMat);
        double sumDT_k2 = tempMat(0, 0);
        Math::Matrix<double>::multiplyMatrices(P_obs_transpose, covth_zerowin, &tempMat);
        double sumDT_zerowin = tempMat(0, 0);

        Math::Matrix<double> mpk_WPth_transpose = mpk_WPth.getTranspose();
        Math::Matrix<double> mpk_WPth_k_transpose = mpk_WPth_k.getTranspose();
        Math::Matrix<double> mpk_WPth_k2_transpose = mpk_WPth_k2.getTranspose();
        Math::Matrix<double>::multiplyMatrices(mpk_WPth_transpose, covth, &tempMat);
        double sumTT = tempMat(0, 0);
        Math::Matrix<double>::multiplyMatrices(mpk_WPth_transpose, covth_k, &tempMat);
        double sumTT_k = tempMat(0, 0);
        Math::Matrix<double>::multiplyMatrices(mpk_WPth_transpose, covth_k2, &tempMat);
        double sumTT_k2 = tempMat(0, 0);
        Math::Matrix<double>::multiplyMatrices(mpk_WPth_k_transpose, covth_k, &tempMat);
        double sumTT_k_k = tempMat(0, 0);
        Math::Matrix<double>::multiplyMatrices(mpk_WPth_k_transpose, covth_k2, &tempMat);
        double sumTT_k_k2 = tempMat(0, 0);
        Math::Matrix<double>::multiplyMatrices(mpk_WPth_k2_transpose, covth_k2, &tempMat);
        double sumTT_k2_k2 = tempMat(0, 0);
        Math::Matrix<double>::multiplyMatrices(mpk_WPth_transpose, covth_zerowin, &tempMat);
        double sumTT_zerowin = tempMat(0, 0);
        Math::Matrix<double>::multiplyMatrices(mpk_WPth_k_transpose, covth_zerowin, &tempMat);
        double sumTT_k_zerowin = tempMat(0, 0);
        Math::Matrix<double>::multiplyMatrices(mpk_WPth_k2_transpose, covth_zerowin, &tempMat);
        double sumTT_k2_zerowin = tempMat(0, 0);
        Math::Matrix<double>::multiplyMatrices(zerowindowfxnsubdat_.getTranspose(), covth_zerowin, &tempMat);
        double sumTT_zerowin_zerowin = tempMat(0, 0);

        // Nuisance parameter integration
        std::vector<double> chisq(nptstot_);
        std::vector<double> chisqmarg(nptstot_);
        double minchisq, maxchisq, deltaL;
        int myminchisqindx;
        double currminchisq, currminchisqmarg, minchisqtheoryamp, chisqnonuis;
        double minchisqtheoryampnonuis, minchisqtheoryampminnuis;
        double a1val, a2val, zerowinsub;
        std::vector<double> a1list(nptstot_);
        std::vector<double> a2list(nptstot_);
        nuisance_init_(a1list, a2list);

        double sumDT_tot, sumTT_tot;
        currminchisq = 1000.0;
        for(int i = 0; i < nptstot_; ++i)
        {
            a1val = a1list[i];
            a2val = a2list[i];
            zerowinsub = -1.0*(sumzerow_Pth + a1val*sumzerow_Pth_k + a2val*sumzerow_Pth_k2);

            sumDT_tot = sumDT + a1val*sumDT_k + a2val*sumDT_k2 + zerowinsub*sumDT_zerowin;
            sumTT_tot = sumTT + pow(a1val, 2.0)*sumTT_k_k + pow(a2val, 2.0)*sumTT_k2_k2 +
                        pow(zerowinsub, 2.0)*sumTT_zerowin_zerowin +
                        2.0*a1val*sumTT_k + 2.0*a2val*sumTT_k2 + 2.0*a1val*a2val*sumTT_k_k2 +
                        2.0*zerowinsub*sumTT_zerowin + 2.0*zerowinsub*a1val*sumTT_k_zerowin +
                        2.0*zerowinsub*a2val*sumTT_k2_zerowin;
            minchisqtheoryamp = sumDT_tot/sumTT_tot;
            chisq[i] = sumDD - 2.0*minchisqtheoryamp*sumDT_tot + pow(minchisqtheoryamp, 2.0)*sumTT_tot;
            chisqmarg[i] = sumDD - pow(sumDT_tot, 2.0)/sumTT_tot + std::log(sumTT_tot) -
                           2.0*std::log(1.0 + erf(sumDT_tot/2.0/sqrt(sumTT_tot)));
            
            if(i == 0 || chisq[i] < currminchisq)
            {
                myminchisqindx = i;
                currminchisq = chisq[i];
                currminchisqmarg = chisqmarg[i];
                minchisqtheoryampminnuis = minchisqtheoryamp;
            }

            if(i == (nptstot_ / 2))
            {
                chisqnonuis = chisq[i];
                minchisqtheoryampnonuis = minchisqtheoryamp;
                check(std::abs(a2val) <= 0.001 && std::abs(a2val) <= 0.001, "a1 or a2 > 0.001 failure");
            }
        }

        // Numerically marginalize over a1, a2 now using values stored in chisq
        minchisq = *std::min_element(std::begin(chisqmarg), std::end(chisqmarg));
        maxchisq = *std::max_element(std::begin(chisqmarg), std::end(chisqmarg));

        double LnLike = 0;
        for(int i = 0; i < nptstot_; ++i)
            LnLike += std::exp(-1.0*(chisqmarg[i]-minchisq)/2.0);
        LnLike = LnLike / (float)nptstot_;
        check(LnLike != 0, "LRG LnLike LogZero error");
        LnLike = -1.0*std::log(LnLike) + minchisq/2.0;
        //deltaL = (maxchisq - minchisq) * 0.5;

        return 2.0*(LnLike);
    }

    void setCosmoParams(const CosmologicalParams& params)
    {
        params_ = &params;
        // Does wantT need to be true?
        cosmo_->initialize(params, true, false, false, true, 0.5);
    }

    void setModelCosmoParams(CosmologicalParams *params)
    {
        modelParams_ = params;
        modelParams_->getAllParameters(vModel_);
    }

private:
    void nuisance_init_(std::vector<double>& a1list, std::vector<double>& a2list)
    {
        a1list.resize(nptstot_);
        a2list.resize(nptstot_);
        double a1val, a2val;
        double da1, da2;
        int countcheck = 0;
        
        da1 = a1maxval_ / (nptsa1_ / 2);
        da2 = a2maxpos_(-1.0*a1maxval_) / (nptsa2_ / 2);
        for(int i = -1 * nptsa1_ / 2; i <= nptsa1_ / 2; ++i)
        {
            for(int j = -1 * nptsa2_ / 2; j <= nptsa2_ / 2; ++j)
            {
                a1val = da1*i;
                a2val = da2*j;
                if((a2val >= 0.0 && a2val <= a2maxpos_(a1val) && a2val >= a2minfinalpos_(a1val)) ||
                   (a2val <= 0.0 && a2val <= a2maxfinalneg_(a1val) && a2val >= a2minneg_(a1val)))
                {
                    check(testa1a2_(a1val, a2val), "a1, a2 failure");
                    check(countcheck <  nptstot_, "countcheck >= nptstot failure");
                    a1list[countcheck] = a1val;
                    a2list[countcheck] = a2val;
                    ++countcheck;
                }
            }
        }
        check(countcheck == nptstot_, "countcheck failure");
    }

    double a2maxpos_(double a1val)
    {
        double a2max = -1.0;
        if(a1val <= std::min(s1_/k1_, s2_/k2_))
            a2max = std::min(s1_/pow(k1_, 2.0) - a1val/k1_, s2_/pow(k2_, 2.0) - a1val/k2_);
        return a2max;
    }

    double a2min1pos_(double a1val)
    {
        double a2min1 = 0.0;
        if(a1val <= 0.0)
            a2min1 = std::max(std::max(-1.0*s1_/pow(k1_, 2.0) - a1val/k1_, -1.0*s2_/pow(k2_, 2.0) - a1val/k2_), 0.0);
        return a2min1;
    }

    double a2min2pos_(double a1val)
    {
        double a2min2 = 0.0;
        if(std::abs(a1val) >= 2.0*s1_/k1_ && a1val <= 0.0)
            a2min2 = pow(a1val, 2.0)/s1_*0.25;
        return a2min2;
    }

    double a2min3pos_(double a1val)
    {
        double a2min3 = 0.0;
        if(std::abs(a1val) >= 2.0*s2_/k2_ && a1val <= 0.0)
            a2min3 = pow(a1val, 2.0)/s2_*0.25;
        return a2min3;
    }

    double a2minfinalpos_(double a1val)
    {
        return std::max(std::max(a2min1pos_(a1val), a2min2pos_(a1val)), a2min3pos_(a1val));
    }

    double a2minneg_(double a1val)
    {
        double a2min;
        if(a1val >= std::max(-1.0*s1_/k1_, -1.0*s2_/k2_))
            a2min = std::max(-1.0*s1_/pow(k1_, 2.0) - a1val/k1_, -1.0*s2_/pow(k2_, 2.0) - a1val/k2_);
        else
            a2min = 1.0;
        return a2min;
    }

    double a2max1neg_(double a1val)
    {
        double a2max1;
        if(a1val >= 0.0)
            a2max1 = std::min(std::min(s1_/pow(k1_, 2.0) - a1val/k1_, s2_/pow(k2_, 2.0) - a1val/k2_), 0.0);
        else
            a2max1 = 0.0;
        return a2max1;
    }

    double a2max2neg_(double a1val)
    {
        double a2max2 = 0.0;
        if(std::abs(a1val) >= 2.0*s1_/k1_ && a1val >= 0.0)
            a2max2 = -1.0*pow(a1val, 2.0)/s1_*0.25;
        return a2max2;
    }

    double a2max3neg_(double a1val)
    {
        double a2max3 = 0.0;
        if(std::abs(a1val) >= 2.0*s2_/k2_ && a1val >= 0.0)
            a2max3 = -1.0*pow(a1val, 2.0)/s2_*0.25;
        return a2max3;
    }

    double a2maxfinalneg_(double a1val)
    {
        return std::min(std::min(a2max1neg_(a1val), a2max2neg_(a1val)), a2max3neg_(a1val));
    }

    double testa1a2_(double a1val, double a2val)
    {
        bool testresult = true;

        double kext = -1.0*a1val/2.0/a2val;
        double diffval = std::abs(a1val*kext + a2val*pow(kext, 2.0));

        if(kext > 0.0 && kext <= k1_ && diffval > s1_)
            testresult = false;
        if(kext > 0.0 && kext <= k2_ && diffval > s2_)
            testresult = false;
        if (std::abs(a1val*k1_ + a2val*pow(k1_, 2.0)) > s1_)
            testresult = false;
        if (std::abs(a1val*k2_ + a2val*pow(k2_, 2.0)) > s2_)
            testresult = false;

        return testresult;
    }

    int lMax_;
    Cosmo* cosmo_;

    const CosmologicalParams* params_;

    CosmologicalParams* modelParams_;
    std::vector<double> vModel_;

    int num_mpk_points_full_;
    int num_mpk_kbands_full_;

    int min_mpk_points_use_;
    int max_mpk_points_use_;
    int min_mpk_kbands_use_;
    int max_mpk_kbands_use_;

    int k_size_;
    int n_size_;
    double zerowindowfxnsubdatnorm_;

    // Data vectors
    std::vector<double> kh_;
    std::vector<double> k_;
    Math::Matrix<double> window_;
    Math::Matrix<double> zerowindowfxn_;
    Math::Matrix<double> zerowindowfxnsubdat_;
    Math::Matrix<double> P_obs_;
    Math::Matrix<double> P_err_;
    Math::Matrix<double> invcov_;

    // Location of data/model files
    std::string root_;

    // Variables associated with nuisance parameters
    double k1_, k2_, s1_, s2_, a1maxval_;
    int nptsa1_, nptsa2_, nptstot_;
};
