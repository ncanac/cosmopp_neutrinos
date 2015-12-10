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

class WiggleZLikelihood : public Math::CosmoLikelihood
{
public:
    WiggleZLikelihood(std::string path, char redshift_bin) 
    {
        cosmo_ = new Cosmo;
        lMax_ = 3000;
        cosmo_->preInitialize(lMax_, false, false, false, lMax_);

        bool Q_marge = false;
        double Q_mid = 4.0;
        double Q_sigma = 1.5;
        double Ag = 1.4;
        bool Use_jennings = false;
        bool Use-simpledamp = false;
        int giggleZ_fidpk_size = 500;
        if(redshift_bin == 'a')
        {
            double d_angular_fid = 736.293;
            double d_radial_fid = 3876.81;
            double sigmav = 360; // or int?
            std::vector<double> giggleZ_fidpoly {4.61900,-13.7787,58.9410,-175.240,284.321,-187.284};
            double redshift = 0.22;
        }
        if(redshift_bin == 'b')
        {
            double d_angular_fid = 1134.87;
            double d_radial_fid = 3511.96;
            double sigmav = 308; // or int?
            std::vector<double> giggleZ_fidpoly {4.63079, -12.6293, 42.9265, -91.8068, 97.808, -37.633};
            double redshift = 0.41;
        }
        if(redshift_bin == 'c')
        {
            double d_angular_fid = 1396.05;
            double d_radial_fid = 3160.38;
            double sigmav = 325; // or int?
            std::vector<double> giggleZ_fidpoly {4.69659, -12.7287, 42.5681, -89.5578, 96.664, -41.2564};
            double redshift = 0.60;
        }
        if(redshift_bin == 'd')
        {
            double d_angular_fid = 1558.68;
            double d_radial_fid = 2852.95;
            double sigmav = 212; // or int?
            std::vector<double> giggleZ_fidpoly {4.6849, -13.4747, 53.7172, -145.832, 216.638, -132.782};
            double redshift = 0.78;
        }

        // Number of points and kbands in the input files
        num_mpk_points_full_ = 50;
        num_mpk_kbands_full_ = 100;

        // Decide which bandpowers to use, min to max
        min_mpk_points_use_ = 3;
        max_mpk_points_use_ = 30;
        min_mpk_kbands_use_ = 1;
        max_mpk_kbands_use_ = 100;

        k_size_ = max_mpk_kbands_use_ - min_mpk_kbands_use_ + 1;
        int mu_size = 1;
        k_.resize(k_size_);
        kh_.resize(k_size_);

        // Read in data file containing kbands
        root_ = path + "data/WiggleZ/";
        std::ifstream datafile(root_ + "wigglez_nov11_kbands.txt");
        for(int i = 0; i < num_mpk_kbands_full_; ++i)
            if((i+2 > min_mpk_kbands_use_) && (i < max_mpk_kbands_use_))
                datafile >> kh_[i-min_mpk_kbands_use_+1];
        datafile.close();

        double khmax = kh_[k_size_-1];

        // Read in data file containing fiducial power spectrum
        // to determine k_fid_size and ifid_discard. These are
        // basically parameters that select the values of k in
        // the fiducial model that match with the values of k
        // from the data.
        int k_fid_size;
        int ifid_discard;
        datafile.open(root_ + giggleZ_fidpk_file "gigglezfiducialmodel_matterpower_" + redshift_bin + ".dat");
        // TODO: Start from line 1238, need to read in data files correctly.
        //std::string line;
        double kdum, dum1, dum2, dum3;
        //std::getline(datafile, line);
        //std::istringstream iss(line);
        datafile >> kdum >> dum1 >> dum2 >> dum3;
        int line_number = 1;
        while(kdum < kh_[0])
        {
            check(line_numer <= 500, "read too many lines");
            //std::getline(datafile, line);
            //std::istringstream iss(line);
            datafile >> kdum >> dum1 >> dum2 >> dum3;
            line_number = line_number + 1;
        }
        ifid_discard = line_number - 2;
        while(kdum < khmax)
        {
            check(line_numer <= 500, "read too many lines");
            //std::getline(datafile, line);
            //std::istringstream iss(line);
            datafile >> kdum >> dum1 >> dum2 >> dum3;
            line_number = line_number + 1;
        }
        datafile.close();
        k_fid_size = line_number - ifid_discard + 1;
        khmax = kdum;

        //khmax = khmax * 2; // If using halofit. Why?

        //int num_regions = 1;
        //int num_regions_used = 1;
        //bool used_region = true;

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
        P_obs_.resize(n_size_, 1, 0);
        P_err_.resize(n_size_, 1, 0);
        //std::vector<double> P_obs(n_size_);
        //std::vector<double> P_err(n_size_);
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
        //std::vector< std::vector<double> > invcov(n_size_, std::vector<double>(n_size_));
        datafile.open(root_ + "data/lrgdr7_invcov.txt");
        for(int i = 0; i < num_mpk_points_full_; ++i)
            if((i+2 > min_mpk_points_use_) && (i < max_mpk_points_use_))
                for(int j = 0; j < num_mpk_points_full_; ++j)
                    if((j+2 > min_mpk_points_use_) && (j < max_mpk_points_use_))
                        datafile >> invcov_(i, j);
        datafile.close();

        // Read in fiducial model again, but different from first time
        k_fid_.resize(k_fid_size, 1, 0);
        Plin_fid_.resize(k_fid_size, 1, 0);
        Psmooth_fid_.resize(k_fid_size, 1, 0);
        rationwhalofit_.resize(k_fid_size, 1, 0);
        //std::vector<double> k_fid(k_fid_size);
        //std::vector<double> Plin_fid(k_fid_size); 
        //std::vector<double> Psmooth_fid(k_fid_size);
        //std::vector<double> rationwhalofit_fid(k_fid_size);
        datafile.open(root_ + "models/lrgdr7fiducialmodel_matterpowerzNEAR.dat"); // Need to do this for NEAR, MID, and FAR? What about z0?
        //Skip first ifid_discard lines
        for(int i = 0; i < ifid_discard; ++i)
            std::getline(datafile, line);
        for(int i = 0; i < k_fid_size; ++i)
        {
            std::getline(datafile, line);
            std::istringstream iss(line);
            double kdummy, plindummy, psmoothdummy, ratiodummy;
            iss >> kdummy, plindummy, psmoothdummy, ratiodummy; 
            k_fid_(i, 0) = kdummy;
            Plin_fid_(i, 0) = plindummy;
            Psmooth_fid_(i, 0) = psmoothdummy;
            rationwhalofit_(i, 0) = ratiodummy;
        }
        datafile.close();
    }

    ~WiggleZLikelihood()
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
    
        return likelihoodBR09();
    }

    double likelihood()
    {
        check(false, "Don't run me!");
        /*
        Data structures used:
        P_lin = a k_size_ x 1 column vector containing the linear power spectrum computed from Class
        P_data = a n_size_ x 1 column vector containing
        W_P_th = a n_size_ x 1 column vector containing
        cov_dat = a n_size_ x 1 column vector containing
        cov_th = a n_size_ x 1 column vector containing
        P_th = a k_size_ x 1 column vector containing
        window_ = an 
        */

        double h = params_->getH();
        // Rescale k_ and P_obs with correct factors of h
        for(int i = 0; i < k_size_; ++i)
            k_[i] = kh_[i]/h;
        for(int i = 0; i < n_size_; ++i)
            P_obs_(i, 0) = P_obs_(i, 0)*h*h*h;

        // Create a vector containing linear matter power spectrum
        // evaluated at values of k given in k_.
        Math::TableFunction<double, double> P_lin_function;
        //cosmo_->getMatterPs(0.3, &P_lin_function);
        cosmo_->getLRGHaloPs(&P_lin_function);

        // Writing P_lin_function to file
        //std::ofstream outTest(root_ + "models/test.txt");
        //for(std::map<double, double>::iterator ps_it = P_lin_function.begin(); ps_it != P_lin_function.end(); ++ps_it)
        //    outTest << ps_it->first << ' ' << ps_it->second << std::endl;
        //for(int i = 0; i < k_size_; ++i)
        //    outTest << k_[i] << std::endl;
        //outTest.close();

        Math::Matrix<double> P_lin(k_size_, 1, 0);
        //std::vector<double> P_lin(k_size_);
        for(int i = 0; i < k_size_; ++i)
            P_lin(i, 0) = P_lin_function.evaluate(k_[i]);

        // Initialize P_nw (linear power spectrum with "no wiggles") based on spline method in BR09
        // Temporary arrays to pass to function dopksmoothbspline_
        //double kvals[k_size_];
        //double lnP_lin[k_size_];
        //double lnP_nw[k_size_];
        // initialize arrays with values in k_ and P_lin
        //for(int i = 0; i < k_size_; ++i)
        //{
        //    kvals[i] = k_[i];
        //    lnP_lin[i] = log(P_lin(i,0));
        //}
        //dopksmoothbspline_(kvals, lnP_lin, lnP_nw, k_size_);
        // Store result of smoothing in P_nw
        //Math::Matrix<double> P_nw(k_size_, 1, 0);
        //for(int i = 0; i < k_size_; ++i)
        //    P_nw(i, 0) = exp(lnP_nw[i]);

        // Writing P_lin and P_nw to file
        //std::ofstream outTest("test.txt");
        //for(int i = 0; i < k_size_; ++i)
        //    outTest << k_[i] << " " << P_lin(i,0) << " " << P_nw(i,0) << std::endl;
        //outTest.close();

        // Compute P_damp according to eq. 10 in BR09
        //const double sigma2BAONEAR = 86.9988, sigma2BAOMID = 85.1374, sigma2BAOFAR = 84.5958 ;
        //Math::Matrix<double> P_damp(k_size_, 1, 0);
        //for(int i = 0; i < k_size_; ++i)
        //    P_damp[i] = P_lin * exp(-1.0 * pow(k_[i], 2) * sigma2BAONEAR * 0.5)
        //                + P_nw * exp(-1.0 * pow(k_[i], 2) * sigma2BAONEAR * 0.5)

        // Compute P_halofit,nw
        //P_halofitnw
        //nonlinear_halofit(P_nw);
        //Math::Matrix<double> Phalofitnw(k_size_, 1, 0);

        // TODO: I think some rescaling goes here. Use fiducial model.

        // Power spectrum data
        Math::Matrix<double> P_data(n_size_, 1, 0);
        //std::vector<double> P_data(n_size_);
        // Windowed matter power spectrum
        Math::Matrix<double> W_P_th(n_size_, 1, 0);
        //std::vector<double> W_P_th(n_size_);
        // Something covariance?
        Math::Matrix<double> cov_dat(n_size_, 1, 0);
        Math::Matrix<double> cov_th(n_size_, 1, 0);
        //std::vector<double> cov_dat(n_size_);
        //std::vector<double> cov_th(n_size_);
        
        Math::Matrix<double> P_th = P_lin;
        //std::vector<double> P_th(k_size_);
        //P_th = P_lin; // Can we just replace P_lin with P_th earlier?

        // Writing P_th to file
        //std::ofstream out(root_ + "models/test_pth.txt");
        //for(int i = 0; i < k_size_; ++i)
        //{
        //    out << k_[i] << " " << P_th(i, 0) << std::endl;
        //}
        //out.close();

        // Fill in theoretical windowed power spectrum
        Math::Matrix<double>::multiplyMatrices(window_, P_th, &W_P_th);
        
        // Set covariance matrices
        Math::Matrix<double>::multiplyMatrices(invcov_, P_obs_, &cov_dat);
        Math::Matrix<double>::multiplyMatrices(invcov_, W_P_th, &cov_th);
        
        // Calculate normV
        double normV = 0;
        Math::Matrix<double> tempMat;
        Math::Matrix<double>::multiplyMatrices(W_P_th.getTranspose(), cov_th, &tempMat);
        normV = normV + tempMat(0, 0);

        // Calculate bias factor
        //double b_out = 0;
        //Math::Matrix<double>::multiplyMatrices(W_P_th.getTranspose(), cov_dat, &tempMat);
        //b_out = b_out + tempMat(0, 0);
        //Math::Matrix<double>::multiplyMatrices(W_P_th.getTranspose(), cov_th, &tempMat);
        //b_out = b_out / tempMat(0, 0);

        // Calculate chi-squared
        //double chisq = 0;
        //Math::Matrix<double> delta(n_size_, 1, 0);
        //for(int i = 0; i < n_size_; ++i)
        //    delta(i, 0) = P_obs_(i, 0) - W_P_th(i, 0);
        //Math::Matrix<double>::multiplyMatrices(delta.getTranspose(), invcov_, &tempMat);
        //Math::Matrix<double> tempMat2;
        //Math::Matrix<double>::multiplyMatrices(tempMat, delta, &tempMat2);
        //chisq = tempMat2(0, 0); 
        double chisq = 0;
        Math::Matrix<double>::multiplyMatrices(P_obs_.getTranspose(), cov_dat, &tempMat);
        chisq = tempMat(0, 0);
        Math::Matrix<double>::multiplyMatrices(W_P_th.getTranspose(), cov_dat, &tempMat); 
        chisq = chisq - pow(tempMat(0, 0), 2.0) / normV;

        // Return -2ln(L)
        return chisq;
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
    //std::vector< std::vector<double> > window_;
    Math::Matrix<double> P_obs_;
    //std::vector<double> P_obs_;
    Math::Matrix<double> P_err_;
    //std::vector<double> P_err_;
    Math::Matrix<double> invcov_;
    //std::vector< std::vector<double> > invcov_;
    Math::Matrix<double> k_fid_;
    //std::vector<double> k_fid_;
    Math::Matrix<double> Plin_fid_;
    //std::vector<double> Plin_fid_;
    Math::Matrix<double> Psmooth_fid_;
    //std::vector<double> Psmooth_fid_;
    Math::Matrix<double> rationwhalofit_;
    //std::vector<double> rationwhalofit_;

    // Location of data/model files
    std::string root_;

    // Variables associated with nuisance parameters
    double k1_, k2_, s1_, s2_, a1maxval_;
    int nptsa1_, nptsa2_, nptstot_;
};
