#ifdef COSMO_MPI
#include <mpi.h>
#endif

#include <string>
#include <fstream>
#include <cmath>

#include <macros.hpp>
#include <cosmo.hpp>
#include <likelihood_function.hpp>
#include <matrix_impl.hpp>
#include <cubic_spline.hpp>

#include <class.h>

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>

class LRGDR7Likelihood : public Math::LikelihoodFunction
{
public:
    LRGDR7Likelihood() 
    {
        cosmo_ = new Cosmo;
        lMax_ = 3000;
        cosmo_->preInitialize(lMax_, false, false, false, lMax_);

        // TODO: set redshift
        redshift_ = 0;

        // Number of points and kbands in the input files
        num_mpk_points_full_ = 45;
        num_mpk_kbands_full_ = 250;

        // Decide which bandpowers to use, min to max
        min_mpk_points_use_ = 1;
        max_mpk_points_use_ = 45;
        min_mpk_kbands_use_ = 1;
        max_mpk_kbands_use_ = 250;

        k_size_ = max_mpk_kbands_use_ - min_mpk_kbands_use_ + 1;
        int mu_size = 1;
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

        double khmax = kh_[k_size_-1];

        // Read in data file containing fiducial power spectrum
        // to determine k_fid_size and ifid_discard. These are
        // basically parameters that select the values of k in
        // the fiducial model that match with the values of k
        // from the data.
        int k_fid_size;
        int ifid_discard;
        datafile.open(root_ + "models/lrgdr7fiducialmodel_matterpowerzNEAR.dat"); // Need to do this for NEAR, MID, and FAR? What about z0?
        // Skip first line
        std::string line;
        std::getline(datafile, line);
        double kdum, dum1, dum2, dum3;
        // Read first "real" line
        std::getline(datafile, line);
        std::istringstream iss(line);
        datafile >> kdum >> dum1 >> dum2 >> dum3;
        int line_number = 1;
        while(kdum < kh_[0])
        {
            std::getline(datafile, line);
            std::istringstream iss(line);
            datafile >> kdum >> dum1 >> dum2 >> dum3;
            line_number = line_number + 1;
        }
        ifid_discard = line_number - 2;
        while(kdum < khmax)
        {
            std::getline(datafile, line);
            std::istringstream iss(line);
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

        output_screen("Checkpoint 0" << std::endl);

        double h = params_->getH();
        // Rescale k_ and P_obs with correct factors of h
        for(int i = 0; i < k_size_; ++i)
            k_[i] = kh_[i]/h;
        for(int i = 0; i < n_size_; ++i)
            P_obs_(i, 0) = P_obs_(i, 0)*h*h*h;

        // Create a vector containing linear matter power spectrum
        // evaluated at values of k given in k_.
        Math::TableFunction<double, double> P_lin_function;
        //cosmo_->getMatterPs(redshift_, &P_lin_function);
        cosmo_->getLRGHaloPs(redshift_, &P_lin_function);

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

        output_screen("Checkpoint 1" << std::endl);

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

        output_screen("Checkpoint 2" << std::endl);

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

        output_screen("Checkpoint 3" << std::endl);

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

        output_screen("Checkpoint 4" << std::endl);

        // Return -2ln(L)
        return chisq;
    }

    double likelihoodBR09()
    {
        //output_screen("Checkpoint 0" << std::endl);
        // num_mpk_points_use = n_size_ = 45
        // num_mpk_kbands_use = k_size_ = 250
        double h = params_->getH();

        Math::Matrix<double> mpk_raw(k_size_, 1, 0);
        Math::Matrix<double> mpk_Pth(k_size_, 1, 0);
        Math::Matrix<double> mpk_Pth_k(k_size_, 1, 0);
        Math::Matrix<double> mpk_Pth_k2(k_size_, 1, 0);
        Math::Matrix<double> mpk_WPth(n_size_, 1, 0);
        Math::Matrix<double> mpk_WPth_k(n_size_, 1, 0);
        Math::Matrix<double> mpk_WPth_k2(n_size_, 1, 0);
        Math::Matrix<double> kh_scaled(k_size_, 1, 0); // h^-1 Mpc^-1

        // TODO: Compute scaling factor
        double a_scl = 1.012937; // Hard coded for now for fiducial model

        // Initialize k_ in units of 1/Mpc
        for(int i = 0; i < k_size_; ++i)
            k_[i] = kh_[i]*h;

        // Initialize halopowerlrgtheory
        //Math::TableFunction<double, double> halopowerlrgtheory;
        // zNEAR = 0.235, zMID = 0.342, zFAR = 0.421  
        // aNEAR = 0.809717d0, aMID = 0.745156d0, aFAR = 0.70373d0
        // zeffDR7 = 0.312782  !! redshift at which a_scl is evaluated.
        //cosmo_->getLRGHaloPs(0.312782, &halopowerlrgtheory);

        //*************************************************************
        // Debugging code

        int numlrg = 300;
        std::vector<double> halokvec(300);
        std::vector<double> halopowervec(300);
        std::ifstream datafile(root_ + "models/halopowerlrgtheory.txt");
        for(int i = 0; i < numlrg; ++i)
        {
            datafile >> halokvec[i];
            datafile >> halopowervec[i];
        }
        datafile.close();
        Math::CubicSpline halopowerlrgtheory(halokvec, halopowervec);

        //std::ifstream datafile(root_ + "models/mpk_raw.txt");
        //for(int i = 0; i < k_size_; ++i)
        //{
        //    datafile >> kh_scaled(i, 0);
        //    datafile >> mpk_raw(i, 0);
        //}
        //datafile.close();

        //*************************************************************

        // TODO: Do this like in BR09
        // Calculate kh_scaled and mpk_raw, which is just halopowerlrgtheory evaluated at kh_scaled*h
        // mpk_raw is in units of h^3 Mpc^3
        for(int i = 0; i < k_size_; ++i)
        {
            kh_scaled(i, 0) = a_scl * kh_[i];
            //mpk_raw(i, 0) = halopowerlrgtheory.evaluate(kh_scaled(i, 0) * h / pow(a_scl, 3.0)) * pow(h, 3.0);
            mpk_raw(i, 0) = halopowerlrgtheory.evaluate(kh_scaled(i, 0)) / pow(a_scl, 3.0);
        }

        
        //output_screen("Checkpoint 1" << std::endl);

        // Initialize
        mpk_Pth = mpk_raw;
        for(int i = 0; i < k_size_; ++i)
        {
            mpk_Pth_k(i, 0) = mpk_Pth(i, 0) * kh_scaled(i, 0);
            mpk_Pth_k2(i, 0) = mpk_Pth(i, 0) * pow(kh_scaled(i, 0), 2.0);
        }

        std::ofstream outTest("test.txt");
        outTest << zerowindowfxnsubdatnorm_ << std::endl;
        for(int i = 0; i < n_size_; ++i)
            outTest << zerowindowfxnsubdat_(i, 0) << std::endl;
            //outTest << k_[i] << " " << " " << zerowindowfxn_(i, 0) << " " <<  mpk_Pth(i,0) << std::endl;
        outTest.close();

        //output_screen("Checkpoint 2" << std::endl);
        Math::Matrix<double>::multiplyMatrices(window_, mpk_Pth, &mpk_WPth);
        Math::Matrix<double>::multiplyMatrices(window_, mpk_Pth_k, &mpk_WPth_k);
        Math::Matrix<double>::multiplyMatrices(window_, mpk_Pth_k2, &mpk_WPth_k2);

        //output_screen("Checkpoint 3" << std::endl);
        Math::Matrix<double> tempMat;
        Math::Matrix<double> zerowindowfxn_transpose = zerowindowfxn_.getTranspose();
        Math::Matrix<double>::multiplyMatrices(zerowindowfxn_transpose, mpk_Pth, &tempMat);
        double sumzerow_Pth = tempMat(0, 0) / zerowindowfxnsubdatnorm_;
        Math::Matrix<double>::multiplyMatrices(zerowindowfxn_transpose, mpk_Pth_k, &tempMat);
        double sumzerow_Pth_k = tempMat(0, 0) / zerowindowfxnsubdatnorm_;
        Math::Matrix<double>::multiplyMatrices(zerowindowfxn_transpose, mpk_Pth_k2, &tempMat);
        double sumzerow_Pth_k2 = tempMat(0, 0) / zerowindowfxnsubdatnorm_;

        //output_screen("Checkpoint 4" << std::endl);
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

        //output_screen("Checkpoint 5" << std::endl);
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

        //output_screen("Checkpoint 6" << std::endl);
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

        //output_screen("Checkpoint 7" << std::endl);
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
        //output_screen("Checkpoint 7.1" << std::endl);
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

        //output_screen(chisq[0] << " " << chisq[150] << " " <<chisq[300] << std::endl);
        //output_screen(chisqmarg[0] << " " << chisqmarg[150] << " " <<chisqmarg[300] << std::endl);

        // Numerically marginalize over a1, a2 now using values stored in chisq
        minchisq = *std::min_element(std::begin(chisqmarg), std::end(chisqmarg));
        maxchisq = *std::max_element(std::begin(chisqmarg), std::end(chisqmarg));

        //output_screen(minchisq << std::endl);
        //output_screen(maxchisq << std::endl);

        double LnLike = 0;
        for(int i = 0; i < nptstot_; ++i)
            LnLike += std::exp(-1.0*(chisqmarg[i]-minchisq)/2.0);
        LnLike = LnLike / (float)nptstot_;
        check(LnLike != 0, "LRG LnLike LogZero error");
        LnLike = -1.0*std::log(LnLike) + minchisq/2.0;
        deltaL = (maxchisq - minchisq) * 0.5;

        return LnLike;
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
    // Copied from BR09 LRG likelihood code, which was copied from http://www.gnu.org/software/gsl/manual/html_node/Basis-Splines.html
    // Procedure:
    // -Fit a cubic b-spline to Plin*k^1.5
    // -Eight equally spaced nodes starting at k = 0.0175 Mpc^-1 to 0.262 Mpc^-1.
    // -One additional node at k = 0.0007 Mpc^-1
    // Inputs:
    // kvals is an array containing the values of k to evaluate the power spectra at
    // lnpklinear is the natural logarithm of the linear power spectrum
    // npts is the length of the input (and output) arrays
    // Outputs:
    // lnpksmooth is the natural logarithm of the smoothed linear power spectrum
    void dopksmoothbspline_(double *kvals, double *lnpklinear, double *lnpksmooth, int npts)
    {
	    double kmaxsuppress = 0.01*0.7;
	    size_t n, ncoeffs, nbreak;
	    gsl_bspline_workspace *bw;
	    gsl_vector *B;
	    gsl_vector *c, *w, *x, *y;
	    gsl_matrix *X, *cov;
	    gsl_multifit_linear_workspace *mw;
	    double deltak,lastk;
	    int i,j,countkeep;

	    nbreak = 9;
	    gsl_vector *mybreaks = gsl_vector_alloc(nbreak);
	    gsl_vector_set(mybreaks,0,(0.001*0.7));
	    gsl_vector_set(mybreaks,1,(0.025*0.7));
	    gsl_vector_set(mybreaks,2,(0.075*0.7));
	    gsl_vector_set(mybreaks,3,(0.125*0.7));
	    gsl_vector_set(mybreaks,4,(0.175*0.7));
	    gsl_vector_set(mybreaks,5,(0.225*0.7));
	    gsl_vector_set(mybreaks,6,(0.275*0.7));
	    gsl_vector_set(mybreaks,7,(0.325*0.7));
	    gsl_vector_set(mybreaks,8,(0.375*0.7));

	    countkeep = 0;
	    for(i=0;i<npts;i++)  {
	    	if((kvals[i]) >= gsl_vector_get(mybreaks,0) && (kvals[i]) <= gsl_vector_get(mybreaks,nbreak-1)) {
	    		countkeep += 1;
	    		}
	    	}
	    n = countkeep;
	    ncoeffs = nbreak + 2;

	    /* allocate a cubic bspline workspace (k = 4) */
	    bw = gsl_bspline_alloc(4, nbreak);
	    B = gsl_vector_alloc(ncoeffs);     
	    x = gsl_vector_alloc(n);
	    y = gsl_vector_alloc(n);
	    X = gsl_matrix_alloc(n, ncoeffs);
	    c = gsl_vector_alloc(ncoeffs);
	    w = gsl_vector_alloc(n);
	    cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
	    mw = gsl_multifit_linear_alloc(n, ncoeffs);
	    i=0;
	    for(j=0;j<npts;j++)  {
	    	if((kvals[j]) >= gsl_vector_get(mybreaks,0) && (kvals[j]) <= gsl_vector_get(mybreaks,nbreak-1)) {
	    		gsl_vector_set(x,i,(kvals[j]));
	    		gsl_vector_set(y,i,exp(lnpklinear[j])*pow(kvals[j],1.5));
	    		if(j>0)  {
	    			deltak = kvals[j] - kvals[j-1];
	    			}
	    		else {
	    			deltak = kvals[0];
	    			if(kvals[1] - kvals[0] < deltak)  {
	    				deltak = kvals[1]-kvals[0];
	    				}
	    			}
	    		gsl_vector_set(w,i,deltak);
	    		i+=1;
	    		}
	    	}
	    gsl_bspline_knots(mybreaks,bw);
	    for(i=0;i<n;i++)  {
	    	double xi = gsl_vector_get(x,i);
	    	gsl_bspline_eval(xi,B,bw);
	    	for(j=0;j<ncoeffs;j++)  {
	    		double Bj = gsl_vector_get(B,j);
	    		gsl_matrix_set(X,i,j,Bj);
	    		}
	    	}
	    //do fit
	    double yi,yierr,chisq;
	    gsl_multifit_wlinear(X,w,y,c,cov,&chisq,mw);
	    i = 0;
	    for(j=0;j<npts;j++)  {
	    	if((kvals[j]) >= gsl_vector_get(mybreaks,0) && (kvals[j]) <= gsl_vector_get(mybreaks,nbreak-1)) {
	    		gsl_bspline_eval(gsl_vector_get(x,i),B,bw);
	    		gsl_multifit_linear_est(B,c,cov,&yi,&yierr);
	    		lnpksmooth[j] = log(yi*pow(kvals[j],-1.5));
	    		i += 1;
	    		}
	    	else {
	    		lnpksmooth[j] = lnpklinear[j];
	    		}
	    	//spline is wacky at small k -- suppress difference at k < 0.01
	    	if(kvals[j] < kmaxsuppress)  {
	    		lnpksmooth[j] = lnpklinear[j];
	    		}
	    	}
	    //assert(i==n);
	    gsl_bspline_free(bw);
	    gsl_vector_free(B);
	    gsl_vector_free(x);
	    gsl_vector_free(y);
	    gsl_vector_free(mybreaks);
	    gsl_matrix_free(X);
	    gsl_vector_free(c);
	    gsl_vector_free(w);
	    gsl_matrix_free(cov);
	    gsl_multifit_linear_free(mw);
    }

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
    double redshift_; // TODO need to initialize this... is it even used though?
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
