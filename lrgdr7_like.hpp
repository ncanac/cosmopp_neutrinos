#ifdef COSMO_MPI
#include <mpi.h>
#endif

#include <string>
#include <fstream>

#include <macros.hpp>
#include <cosmo.hpp>
#include <likelihood_function.hpp>

class LRGDR7Likelihood : public Math::LikelihoodFunction
{
    public:
    LRGDR7Likelihood() 
    {
        cosmo_ = new Cosmo;
        lMax_ = 2;
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

        // Read in data file containing kbands
        std::ifstream datafile("lrgdr7_kbands.txt");
        for(int i = 0; i < num_mpk_kbands_full_; ++i)
            if((i+2 > min_mpk_kbands_use_) && (i < self.max_mpk_kbands_use_))
                datafile >> kh[i-min_mpk_kbands_use_+1];
        datafile.close();

        double khmax = kh[k_size_-1];

        // Read in data file containing fiducial power spectrum
        // to determine k_fid_size and ifid_discard. These are
        // basically parameters that select the values of k in
        // the fiducial model that match with the values of k
        // from the data.
        int k_fid_size;
        int ifid_discard;
        datafile.open("lrgdr7fiducialmodel_matterpowerzNEAR.dat"); // Need to do this for NEAR, MID, and FAR? What about z0?
        // Skip first line
        std::getline(datafile, line);
        double kdum, dum1, dum2, dum3;
        // Read first "real" line
        std::getline(datafile, line);
        std::istringstream iss(line);
        datafile >> kdum >> dum1 >> dum2 >> dum3;
        int line_number = 1;
        while(kdum < kh[0])
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
        n_size_ = max_mpk_points_use_ - min_mpk_points_use_ + 1;
        window_.resize(n_size_, k_size_, 0);
        //std::vector< std::vector<double> > window(n_size_, std::vector<double>(k_size_));
        datafile.open("lrgdr7_windows.txt");
        for(int i = 0; i < num_mpk_points_use_; ++i)
            for(int j = 0; j < k_size_; ++j)
                datafile >> window_(i, j);
        datafile.close()

        // Read in measurements
        P_obs_.resize(n_size_, 1, 0);
        P_err_.resize(n_size_, 1, 0);
        //std::vector<double> P_obs(n_size_);
        //std::vector<double> P_err(n_size_);
        datafile.open("lrgdr7_ccmeasurements.txt");
        std::string line;
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
                P_obs(i-min_mpk_points_use_+1, 0) = pdum;
                P_err(i-min_mpk_points_use_+1, 0) = errdum;
            }
        }
        datafile.close();

        // Read in inverse covariance matrix
        invcov_.resize(n_size_, n_size_, 0);
        //std::vector< std::vector<double> > invcov(n_size_, std::vector<double>(n_size_));
        datafile.open("lrgdr7_invcov.txt");
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
        datafile.open("lrgdr7fiducialmodel_matterpowerzNEAR.dat"); // Need to do this for NEAR, MID, and FAR? What about z0?
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
        datafile.close()
    }

    ~LRGDR7Likelihood()
    {
        delete cosmo_; 
    }

    double calculate(double* params, int nPar)
    {
        
    }

    double likelihood()
    {
        /*
        Data structures defined in this function:
        P_lin = a k_size_ x 1 column vector containing the linear power spectrum computed from Class
        P_data = a n_size_ x 1 column vector containing
        W_P_th = a n_size_ x 1 column vector containing
        cov_dat = a n_size_ x 1 column vector containing
        cov_th = a n_size_ x 1 column vector containing
        P_th = a k_size_ x 1 column vector containing
        window_ = an 
        */
        double h = modelParams_.getH();

        // Create a vector containing linear matter power spectrum
        // evaluated at values of k given in k_.
        Math::TableFunction<double, double> P_lin_function;
        cosmo.getMatterPs(redshift_, &P_lin_function);
        Math::Matrix<double> P_lin(k_size_, 1, 0);
        //std::vector<double> P_lin(k_size_);
        for(i = 0; i < k_size_; ++i)
            P_lin(i, 0) = P_lin_function.evaluate(k_[i]);

        // TODO: I think rescaling goes here. Use fiducial model.

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
        
        Math::Matrix<double> P_th(k_size_, 1, 0);
        //std::vector<double> P_th(k_size_);
        P_th = P_lin; // Can we just replace P_lin with P_th earlier?

        // Fill in theoretical windowed power spectrum
        Math::Matrix<double>::multiplyMatrices(window_, P_th, &W_P_th);
        
        // Set covariance matrices
        Math::Matrix::multiplyMatrices(invcov_, P_obs_, &cov_dat);
        Math::Matrix::multiplyMatrices(invcov_, W_P_th, &cov_th);
        
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
        double chisq = 0;
        Math::Matrix<double>::multiplyMatrices(P_obs_.getTranspose(), cov_dat, &tempMat);
        chisq = tempMat(0, 0);
        Math::Matrix<double>::multiplyMatrices(W_P_th, cov_dat, &tempMat); 
        chisq = chisq - pow(tempMat(0, 0), 2.0) / normV;

        return chisq;
    }

    void setCosmoParams(const CosmologicalParams& params)
    {
        params_ = &params;
        // Does wantT need to be true?
        cosmo_->initialize(params, true, false, false, true);
    }

    void setModelCosmoParams(CosmologicalParams *params)
    {
        modelParams_ = params;
        modelParams_->getAllParameters(vModel_);
    }

private;
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
    double redshift_; // TODO need to initialize this

    // Data vectors
    std::vector<double> kh_;
    std::vector<double> k_;
    Math::Matrix<double> window_;
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
};
