#pragma once
#include <fstream>
#include <sstream>
#include <cmath>

#include <cmb.hpp>
#include <cubic_spline.hpp>
#include <class.h>

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>

class Cosmo : public CMB
{
public:
    // Constructor
    Cosmo() {}

    // Destructor
    ~Cosmo() {}

    double getrsdrag()
    {
        check(init_, "need to initialize first");
        return th_->rs_d;   
    } 

    double getAngularDistance(double z)
    {
        double tau;
        int last_index;
        double* pvecback;

        pvecback = (double*) calloc(br_->bg_size, sizeof(double));
        background_tau_of_z(br_, z, &tau);
        background_at_tau(br_, tau, br_->long_info, br_->inter_normal, &last_index, pvecback);

        const double res = pvecback[br_->index_bg_ang_distance];
        free(pvecback);
        return res;
    }

    double getHubble(double z)
    {
        double tau;
        int last_index;
        double* pvecback;

        pvecback = (double*) calloc(br_->bg_size, sizeof(double));
        background_tau_of_z(br_, z, &tau);
        background_at_tau(br_, tau, br_->long_info, br_->inter_normal, &last_index, pvecback);
        
        const double res = pvecback[br_->index_bg_H];
        free(pvecback);
        return res;
    }

    std::vector<double> z_of_r(std::vector<double>& z_array)
    {
        double tau = 0.0;
        int last_index = 0;
        double* pvecback;
        std::vector<double> r(z_array.size(), 0);
        std::vector<double> dzdr(z_array.size(), 0);
        std::vector<double> out(2*z_array.size(), 0);

        pvecback = (double*) calloc(br_->bg_size, sizeof(double));

        for(int i = 0; i < z_array.size(); ++i)
        {
            background_tau_of_z(br_, z_array[i], &tau);

            background_at_tau(br_, tau, br_->long_info, br_->inter_normal, &last_index, pvecback);

            // store r
            r[i] = pvecback[br_->index_bg_conf_distance];
            // store dz/dr = H
            dzdr[i] = pvecback[br_->index_bg_H];
        }

        for(int i = 0; i < z_array.size(); ++i)
        {
            out[i] = r[i];
            out[i+z_array.size()] = dzdr[i];
        }

        free(pvecback);
        return out;
    }

    void getMatterPsNL(double z, Math::TableFunction<double, double>* ps)
    {
        StandardException exc;
        check(init_, "need to initialize first");
        check(pt_->has_pk_matter, "matter ps not requested");
        check(z >= 0 && z <= sp_->z_max_pk, "invalid z = " << z);
        const int kSize = sp_->ln_k_size;
        check(kSize > 0, "");
    
        // To be done better
        check(kSize < 10000, "");
        double outTot[100000];
        double outIc[100000];
        double outTotNL[100000];
        double tau;
        double k_nl;
    
        if(spectra_pk_at_z(br_, sp_, linear, z, outTot, outIc) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: spectra_pk_at_z failed!" << std::endl << sp_->error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        if(background_tau_of_z(br_, z, &tau) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: background_tau_of_z failed!" << std::endl << br_->error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        if(nonlinear_halofit(pr_, br_, pm_, nl_, tau, outTot, outTotNL, &k_nl) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: nonlinear_halofit failed!" << std::endl << nl_->error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        ps->clear();
    
        for(int i = 0; i < kSize; ++i)
        {
            const double k = std::exp(sp_->ln_k[i]);
            (*ps)[k] = outTotNL[i];
        }
    }

    void getMatterPsNL2(double z, Math::TableFunction<double, double>* ps)
    {
        StandardException exc;
        check(init_, "need to initialize first");
        check(pt_->has_pk_matter, "matter ps not requested");
        check(z >= 0 && z <= sp_->z_max_pk, "invalid z = " << z);
        const int kSize = sp_->ln_k_size;
        check(kSize > 0, "");
    
        // To be done better
        check(kSize < 10000, "");
        double outTot[100000];
    
        if(spectra_pk_nl_at_z(br_, sp_, linear, z, outTot) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: spectra_pk_nl_at_z failed!" << std::endl << sp_->error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        ps->clear();
    
        for(int i = 0; i < kSize; ++i)
        {
            const double k = std::exp(sp_->ln_k[i]);
            (*ps)[k] = outTot[i];
        }
    }

    double getPkNLatk(double k, double z)
    {
        StandardException exc;
        check(init_, "need to initialize first");
        check(pt_->has_pk_matter, "matter ps not requested");
        check(z >= 0 && z <= sp_->z_max_pk, "invalid z = " << z);

        double pk;
    
        if(spectra_pk_nl_at_k_and_z(br_, pm_, sp_, k, z, &pk) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: spectra_pk_nl_at_k_and_z failed!" << std::endl << sp_->error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        return pk;
    }

    bool getLRGHaloPs(std::string root, Math::TableFunction<double, double>* ps)
    {
        StandardException exc;
        check(init_, "need to initialize first");
        check(pt_->has_pk_matter, "matter ps not requested");

        double h = params_->getH();

        const double zNEAR = 0.235, zMID = 0.342, zFAR = 0.421;
        //const double aNEAR = 0.809717, aMID = 0.745156, aFAR = 0.70373;

        // Get linear matter power spectrum for NEAR, MID, and FAR
        // Multiply by h^3 to be equivalent to CAMB output
        Math::TableFunction<double, double> P_lin_funcNEAR;
        Math::TableFunction<double, double> P_lin_funcMID;
        Math::TableFunction<double, double> P_lin_funcFAR;
        Math::TableFunction<double, double> P_lin_funcz0;
        getMatterPs(zNEAR, &P_lin_funcNEAR);
        getMatterPs(zMID, &P_lin_funcMID);
        getMatterPs(zFAR, &P_lin_funcFAR);
        getMatterPs(0.0, &P_lin_funcz0);

        // Create some arrays to store the various power spectra and initialize them
        int n = P_lin_funcNEAR.size();
        double kvals[n];
        double P_linNEAR[n];
        double P_linMID[n];
        double P_linFAR[n];
        double P_linz0[n];
        double lnP_linNEAR[n];
        double lnP_linMID[n];
        double lnP_linFAR[n];
        int itemp = 0;
        for(auto const &point : P_lin_funcNEAR)
        {
            kvals[itemp] = point.first;
            double Pval = point.second;
            P_linNEAR[itemp] = Pval;
            lnP_linNEAR[itemp] = std::log(Pval);
            ++itemp;
        }
        check(n == itemp, "number of k values mismatch");
        itemp = 0;
        for(auto const &point : P_lin_funcMID)
        {
            check(std::abs(kvals[itemp] - point.first) < 0.1*kvals[itemp], "kvals should be identical");
            double Pval = point.second;
            P_linMID[itemp] = Pval;
            lnP_linMID[itemp] = std::log(Pval);
            ++itemp;
        }
        check(n == itemp, "number of k values mismatch");
        itemp = 0;
        for(auto const &point : P_lin_funcFAR)
        {
            check(std::abs(kvals[itemp] - point.first) < 0.1*kvals[itemp], "kvals should be identical");
            double Pval = point.second;
            P_linFAR[itemp] = Pval;
            lnP_linFAR[itemp] = std::log(Pval);
            ++itemp;
        }
        check(n == itemp, "number of k values mismatch");
        itemp = 0;
        for(auto const &point : P_lin_funcz0)
        {
            check(std::abs(kvals[itemp] - point.first) < 0.1*kvals[itemp], "kvals should be identical");
            double Pval = point.second;
            P_linz0[itemp] = Pval;
            ++itemp;
        }
        check(n == itemp, "number of k values mismatch");

        // Populate sparse vectors with k values from original BR09 LRG code
        // Read in k values
        int n_sparse = 106;
        double kvals_sparse[n];
        std::ifstream infile("/Users/ncanac/lrgdr7like/models/kvals_sparse.txt");
        for(int i = 0; i < 106; ++i)
        {
            infile >> kvals_sparse[i];
        }
        infile.close();

        check(kvals_sparse[0] > kvals[0], "kmin outside range");
        check(kvals_sparse[n_sparse-1] < kvals[n-1], "max outside range");

        // Initialize linear power spectrum vectors
        double P_linNEAR_sparse[n_sparse];
        double P_linMID_sparse[n_sparse];
        double P_linFAR_sparse[n_sparse];
        double P_linz0_sparse[n_sparse];
        double lnP_linNEAR_sparse[n_sparse];
        double lnP_linMID_sparse[n_sparse];
        double lnP_linFAR_sparse[n_sparse];
        for(int i = 0; i < n_sparse; ++i)
        {
            P_linNEAR_sparse[i] = P_lin_funcNEAR.evaluate(kvals_sparse[i]);
            lnP_linNEAR_sparse[i] = std::log(P_linNEAR_sparse[i]);
            P_linMID_sparse[i] = P_lin_funcMID.evaluate(kvals_sparse[i]);
            lnP_linMID_sparse[i] = std::log(P_linMID_sparse[i]);
            P_linFAR_sparse[i] = P_lin_funcFAR.evaluate(kvals_sparse[i]);
            lnP_linFAR_sparse[i] = std::log(P_linFAR_sparse[i]);
            P_linz0_sparse[i] = P_lin_funcz0.evaluate(kvals_sparse[i]);
        }

        // extract getabstransferscale for NEAR, MID, and FAR
        const double khmindata = 0.02;
        std::vector<double> getabstransferscale(4);
        itemp = 0;
        while(kvals_sparse[itemp]/h < khmindata && itemp < n_sparse)
            ++itemp;
        check(itemp < n_sparse, "khmindata > kvals");
        getabstransferscale[0] = sqrt(P_linNEAR_sparse[itemp]*pow(h, 3.0));        
        getabstransferscale[1] = sqrt(P_linMID_sparse[itemp]*pow(h, 3.0));        
        getabstransferscale[2] = sqrt(P_linFAR_sparse[itemp]*pow(h, 3.0));        
        getabstransferscale[3] = sqrt(P_linz0_sparse[itemp]*pow(h, 3.0));        

        // Compute smoothed power spectrum using basis spline
        double lnP_nwNEAR_sparse[n_sparse];
        double lnP_nwMID_sparse[n_sparse];
        double lnP_nwFAR_sparse[n_sparse];
        double P_nwNEAR_sparse[n_sparse];
        double P_nwMID_sparse[n_sparse];
        double P_nwFAR_sparse[n_sparse];
        //for(int i = 0; i < n; ++i)
        //    output_screen(kvals[i] << " " << lnP_linNEAR[i] << std::endl);
        dopksmoothbspline_(kvals_sparse, lnP_linNEAR_sparse, lnP_nwNEAR_sparse, n_sparse);
        //for(int i = 0; i < n; ++i)
        //    output_screen(kvals[i] << " " << lnP_nwNEAR[i] << std::endl);
        dopksmoothbspline_(kvals_sparse, lnP_linMID_sparse, lnP_nwMID_sparse, n_sparse);
        dopksmoothbspline_(kvals_sparse, lnP_linFAR_sparse, lnP_nwFAR_sparse, n_sparse);
        for(int i = 0; i < n; ++i)
        {
            P_nwNEAR_sparse[i] = std::exp(lnP_nwNEAR_sparse[i]);
            P_nwMID_sparse[i] = std::exp(lnP_nwMID_sparse[i]);
            P_nwFAR_sparse[i] = std::exp(lnP_nwFAR_sparse[i]);
        }

        std::ofstream outfile("P_lin_nw_sparsezNEAR.txt");
        for(int i = 0; i < n_sparse; ++i)
            outfile << kvals_sparse[i] << " " << P_linNEAR_sparse[i] << " " << P_nwNEAR_sparse[i] << std::endl;
        outfile.close();
        outfile.open("P_lin_nw_sparsezMID.txt");
        for(int i = 0; i < n_sparse; ++i)
            outfile << kvals_sparse[i] << " " << P_linMID_sparse[i] << " " << P_nwMID_sparse[i] << std::endl;
        outfile.close();
        outfile.open("P_lin_nw_sparsezFAR.txt");
        for(int i = 0; i < n_sparse; ++i)
            outfile << kvals_sparse[i] << " " << P_linFAR_sparse[i] << " " << P_nwFAR_sparse[i] << std::endl;
        outfile.close();

        // Initialize interpolating function for P_nw
        Math::TableFunction<double, double> lnP_nwNEAR_sparse_spline;
        Math::TableFunction<double, double> lnP_nwMID_sparse_spline;
        Math::TableFunction<double, double> lnP_nwFAR_sparse_spline;
        for(int i = 0; i < n_sparse; ++i)
        {
            lnP_nwNEAR_sparse_spline[std::log(kvals_sparse[i])] = lnP_nwNEAR_sparse[i];
            lnP_nwMID_sparse_spline[std::log(kvals_sparse[i])] = lnP_nwMID_sparse[i];
            lnP_nwFAR_sparse_spline[std::log(kvals_sparse[i])] = lnP_nwFAR_sparse[i];
        }

        // Now initialize full smoothed power spectrum by interpolating
        // For values outside interpolating range, just use linear power spectrum
        //double lnP_nwNEAR[n];
        //double lnP_nwMID[n];
        //double lnP_nwFAR[n];
        double P_nwNEAR[n];
        double P_nwMID[n];
        double P_nwFAR[n];
        for(int i = 0; i < n; ++i)
        {
            if(kvals[i] > kvals_sparse[0] && kvals[i] < kvals_sparse[n_sparse-1])
            {
                P_nwNEAR[i] = std::exp(lnP_nwNEAR_sparse_spline[std::log(kvals[i])]);
                P_nwMID[i] = std::exp(lnP_nwMID_sparse_spline[std::log(kvals[i])]);
                P_nwFAR[i] = std::exp(lnP_nwFAR_sparse_spline[std::log(kvals[i])]);
            }
            else
            {
                P_nwNEAR[i] = P_linNEAR[i];
                P_nwMID[i] = P_linMID[i];
                P_nwFAR[i] = P_linFAR[i];
            }
        }

        // Debugging: Read in P_lin from CAMB for fiducial model
        //int n = 106;
        //double kvals[n];
        //double P_linNEAR[n];
        //double P_linMID[n];
        //double P_linFAR[n];
        //double P_linz0[n];
        //double lnP_linNEAR[n];
        //double lnP_linMID[n];
        //double lnP_linFAR[n];
        //// Read in NEAR model
        //std::ifstream infile("/Users/ncanac/lrgdr7like/models/plinzNEAR.txt");
        //double ktemp, ptemp;
        //int itemp = 0;
        //while(infile >> ktemp >> ptemp)
        //{
        //    kvals[itemp] = ktemp;
        //    P_linNEAR[itemp] = ptemp/pow(h, 3);
        //    lnP_linNEAR[itemp] = std::log(ptemp/pow(h, 3));
        //    ++itemp;
        //}
        //infile.close();
        //check(itemp == n, "Read in wrong number of kvals.");
        //// Read in MID model
        //infile.open("/Users/ncanac/lrgdr7like/models/plinzMID.txt");
        //itemp = 0;
        //while(infile >> ktemp >> ptemp)
        //{
        //    P_linMID[itemp] = ptemp/pow(h, 3);
        //    lnP_linMID[itemp] = std::log(ptemp/pow(h, 3));
        //    ++itemp;
        //}
        //infile.close();
        //// Read in FAR model
        //infile.open("/Users/ncanac/lrgdr7like/models/plinzFAR.txt");
        //itemp = 0;
        //while(infile >> ktemp >> ptemp)
        //{
        //    P_linFAR[itemp] = ptemp/pow(h, 3);
        //    lnP_linFAR[itemp] = std::log(ptemp/pow(h, 3));
        //    ++itemp;
        //}
        //infile.close();
        //// Read in z=0 model
        //infile.open("/Users/ncanac/lrgdr7like/models/plinz0.txt");
        //itemp = 0;
        //while(infile >> ktemp >> ptemp)
        //{
        //    P_linz0[itemp] = ptemp/pow(h, 3);
        //    ++itemp;
        //}
        //infile.close();

        ////////////////////////////////////////////////////////

        // extract getabstransferscale for NEAR, MID, and FAR
        //const double khmindata = 0.0221168;
        //const double khmindata = 0.02;
        //std::vector<double> getabstransferscale(4);
        //double getabstransferscaleNEAR, getabstransferscaleMID, getabstransferscaleFAR;
        //itemp = 0;
        //while(kvals[itemp]/h < khmindata && itemp < n)
        //    ++itemp;
        //check(itemp < n, "khmindata > kvals");
        //getabstransferscale[0] = sqrt(P_linNEAR[itemp]*pow(h, 3.0));        
        //getabstransferscale[1] = sqrt(P_linMID[itemp]*pow(h, 3.0));        
        //getabstransferscale[2] = sqrt(P_linFAR[itemp]*pow(h, 3.0));        
        //getabstransferscale[3] = sqrt(P_linz0[itemp]*pow(h, 3.0));        

        ////std::ofstream outfile("P_linFAR.txt");
        ////for(int i = 0; i < n; ++i)
        ////    outfile << kvals[i] << " " << P_linFAR[i]*pow(h, 3) << std::endl;
        ////outfile.close();
       
        //// Initialize P_nw (linear power spectrum with "no wiggles") based on spline method used in BR09        
        //double lnP_nwNEAR[n];
        //double lnP_nwMID[n];
        //double lnP_nwFAR[n];
        //double P_nwNEAR[n];
        //double P_nwMID[n];
        //double P_nwFAR[n];
        ////for(int i = 0; i < n; ++i)
        ////    output_screen(kvals[i] << " " << lnP_linNEAR[i] << std::endl);
        //dopksmoothbspline_(kvals, lnP_linNEAR, lnP_nwNEAR, n);
        ////for(int i = 0; i < n; ++i)
        ////    output_screen(kvals[i] << " " << lnP_nwNEAR[i] << std::endl);
        //dopksmoothbspline_(kvals, lnP_linMID, lnP_nwMID, n);
        //dopksmoothbspline_(kvals, lnP_linFAR, lnP_nwFAR, n);
        //for(int i = 0; i < n; ++i)
        //{
        //    P_nwNEAR[i] = std::exp(lnP_nwNEAR[i]);
        //    P_nwMID[i] = std::exp(lnP_nwMID[i]);
        //    P_nwFAR[i] = std::exp(lnP_nwFAR[i]);
        //}

        //for(int i = 0; i < n; ++i)
        //    output_screen(kvals[i] << " " << P_nwNEAR[i] << std::endl);

        // Apply halofit model for nonlinear structure growth to P_nw to generate P_halofitnw
        // Tau is the conformal time at z
        // k_nl is the value of k where power spectrum becomes nonlinear
        double P_halofitnwNEAR[n];
        double P_halofitnwMID[n];
        double P_halofitnwFAR[n];
        double tau;
        double k_nl;
        bool halofit_success = true;
        // Apply halofit model to P_nw for NEAR
        //nonlinear_k_nl_at_z(br_, nl_, zNEAR, &k_nl);
        //background_tau_of_z(br_, zNEAR, &tau);
        //nonlinear_halofit(pr_, br_, pm_, nl_, tau, P_nwNEAR, P_halofitnwNEAR, &k_nl);
        if(background_tau_of_z(br_, zNEAR, &tau) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: background_tau_of_z failed!" << std::endl << br_->error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        //for(int i = 0; i < n; ++i)
        //    output_screen(kvals[i] << " " << P_nwNEAR[i] << std::endl);

        if(nonlinear_halofit(pr_, br_, pm_, nl_, tau, P_nwNEAR, P_halofitnwNEAR, &k_nl) == _FAILURE_)
        {
            halofit_success = false;
            for(int i = 0; i < n; ++i)
                P_halofitnwNEAR[i] = P_nwNEAR[i];
        }

        for(int i = 0; i < n; ++i)
            output_screen(kvals[i]/h << " " << P_halofitnwNEAR[i] << std::endl);

        // Apply halofit model to P_nw for MID 
        //nonlinear_k_nl_at_z(br_, nl_, zMID, &k_nl);
        //background_tau_of_z(br_, zMID, &tau);
        //nonlinear_halofit(pr_, br_, pm_, nl_, tau, P_nwMID, P_halofitnwMID, &k_nl);
        if(background_tau_of_z(br_, zMID, &tau) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: background_tau_of_z failed!" << std::endl << br_->error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }
        if(nonlinear_halofit(pr_, br_, pm_, nl_, tau, P_nwMID, P_halofitnwMID, &k_nl) == _FAILURE_)
        {
            halofit_success = false;
            for(int i = 0; i < n; ++i)
                P_halofitnwMID[i] = P_nwMID[i];
        }
        // Apply halofit model to P_nw for FAR 
        //nonlinear_k_nl_at_z(br_, nl_, zFAR, &k_nl);
        //background_tau_of_z(br_, zFAR, &tau);
        //nonlinear_halofit(pr_, br_, pm_, nl_, tau, P_nwFAR, P_halofitnwFAR, &k_nl);
        if(background_tau_of_z(br_, zFAR, &tau) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: background_tau_of_z failed!" << std::endl << br_->error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }
        if(nonlinear_halofit(pr_, br_, pm_, nl_, tau, P_nwFAR, P_halofitnwFAR, &k_nl) == _FAILURE_)
        {
            halofit_success = false;
            for(int i = 0; i < n; ++i)
                P_halofitnwFAR[i] = P_nwFAR[i];
        }

        if(!halofit_success)
        {
            output_screen("WARNING: nonlinear_halofit failed!" << std::endl);
        }

        // Make sure that P_halofitnw is positive everywhere
        //for(int i = 0; i < n; ++i)
        //{
        //    if(P_halofitnwNEAR[i] <= 0 || P_halofitnwMID[i] <= 0 || P_halofitnwFAR[i] <= 0)
        //    {
        //        output_screen("BADHALOFIT" << std::endl);
        //        return false;
        //    }
        //}
        //output_screen("GOODHALOFIT" << std::endl);

        // Calculate factor r_halofit, the ratio of P_halofitnw to P_nw
        std::vector<double> r_halofitNEAR(n);
        std::vector<double> r_halofitMID(n);
        std::vector<double> r_halofitFAR(n);
        for(int i = 0; i < n; ++i)
        {
            r_halofitNEAR[i] = P_halofitnwNEAR[i]/P_nwNEAR[i];
            r_halofitMID[i] = P_halofitnwMID[i]/P_nwMID[i];
            r_halofitFAR[i] = P_halofitnwFAR[i]/P_nwFAR[i];
            //r_halofitNEAR[i] = getPkNLatk(kvals[i], zNEAR)/P_nwNEAR[i];
            //r_halofitMID[i] = getPkNLatk(kvals[i], zMID)/P_nwMID[i];
            //r_halofitFAR[i] = getPkNLatk(kvals[i], zFAR)/P_nwFAR[i];
        }

        // Initialize all arrays to be passed to LRGTheory_()
        std::vector<double> kh(n);
        for(int i = 0; i < n; ++i)
        {
            kh[i] = kvals[i] / h;
        }
        std::vector< std::vector<double> > P_lin(3, std::vector<double>(n));
        std::vector< std::vector<double> > P_nw(3, std::vector<double>(n));
        std::vector< std::vector<double> > r_nwhalofit(3, std::vector<double>(n));
        for(int i = 0; i < n; ++i)
        {
            P_lin[0][i] = P_linNEAR[i] * pow(h, 3.0);
            P_nw[0][i] = P_nwNEAR[i] * pow(h, 3.0);
            r_nwhalofit[0][i] = r_halofitNEAR[i];
        }
        for(int i = 0; i < n; ++i)
        {
            P_lin[1][i] = P_linMID[i] * pow(h, 3.0);
            P_nw[1][i] = P_nwMID[i] * pow(h, 3.0);
            r_nwhalofit[1][i] = r_halofitMID[i];
        }
        for(int i = 0; i < n; ++i)
        {
            P_lin[2][i] = P_linFAR[i] * pow(h, 3.0);
            P_nw[2][i] = P_nwFAR[i] * pow(h, 3.0);
            r_nwhalofit[2][i] = r_halofitFAR[i];
        }

        std::vector<double> kh_fid;
        std::vector<double> P_halo;
        // Calculate P_halo by calling LRGTheory_()
        LRGTheory_(root, kh, P_lin, P_nw, r_nwhalofit, getabstransferscale, kh_fid, P_halo);

        // Store LRG power spectrum in table function
        ps->clear();
        for(int i = 0; i < kh_fid.size(); ++i)
            (*ps)[kh_fid[i]] = P_halo[i];

        return true;
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

    // Copied form BR09 likelihood code
    // Hard coding of polynomial fits to near, mid, and far subsamples.
    // k is really k/h
    void LRGtoICsmooth_(double k, std::vector<double>& fidpolys)
    {
        double fidNEAR, fidMID, fidFAR;
        fidpolys.resize(3);

        if(k < 0.194055)
            fidNEAR = (1.0 - 0.680886*k + 6.48151*k*k);
        else
            fidNEAR = (1.0 - 2.13627*k + 21.0537*k*k - 50.1167*pow(k, 3) + 36.8155*pow(k,4)) * 1.04482;

        if(k < 0.19431)
            fidMID = (1.0 - 0.530799*k + 6.31822*k*k);
        else
            fidMID = (1.0 - 1.97873*k + 20.8551*k*k - 50.0376*pow(k, 3) + 36.4056*pow(k, 4)) * 1.04384;

        if(k < 0.19148)
            fidFAR = (1.0 - 0.475028*k + 6.69004*k*k);
        else
            fidFAR = (1.0 - 1.84891*k + 21.3479*k*k - 52.4846*pow(k, 3) + 38.9541*pow(k, 4)) * 1.03753;

        fidpolys[0] = fidNEAR;
        fidpolys[1] = fidMID;
        fidpolys[2] = fidFAR;
    }

    // Calculates the theoretical halo power spectrum
    void LRGTheory_(std::string root, const std::vector<double>& kh, const std::vector< std::vector<double> >& P_lin,
                    const std::vector< std::vector<double> >& P_nw, const std::vector< std::vector<double> >& r_nwhalofit,
                    const std::vector<double>& getabstransferscale, std::vector<double>& kh_fid, std::vector<double>& P_halo)
    {
        const int k_size = 300;
        std::vector< std::vector<double> > r_fid(3, std::vector<double>(k_size, 0.0));
        kh_fid.resize(k_size);
        // Read in NEAR model
        //std::ifstream datafile(root + "models/lrgdr7fiducialmodel_matterpowerzNEAR.dat");
        std::ifstream datafile(root + "models/lrgdr7fiducialmodel_matterpowerzNEAR_new.dat");
        // Skip first line
        std::string line;
        std::getline(datafile, line);
        // Read in rest of data file
        for(int i = 0; i < k_size; ++i)
        {
            std::getline(datafile, line);
            std::istringstream iss(line);
            double kdummy, plindummy, psmoothdummy, ratiodummy;
            iss >> kdummy >> plindummy >> psmoothdummy >> ratiodummy; 
            kh_fid[i] = kdummy; // Only need to do this once
            r_fid[0][i] = ratiodummy;
        }
        datafile.close();
        // Read in MID model
        //datafile.open(root + "models/lrgdr7fiducialmodel_matterpowerzMID.dat");
        datafile.open(root + "models/lrgdr7fiducialmodel_matterpowerzMID_new.dat");
        // Skip first line
        std::getline(datafile, line);
        // Read in rest of data file
        for(int i = 0; i < k_size; ++i)
        {
            std::getline(datafile, line);
            std::istringstream iss(line);
            double kdummy, plindummy, psmoothdummy, ratiodummy;
            iss >> kdummy >> plindummy >> psmoothdummy >> ratiodummy; 
            r_fid[1][i] = ratiodummy;
        }
        datafile.close();
        // Read in FAR model
        //datafile.open(root + "models/lrgdr7fiducialmodel_matterpowerzFAR.dat");
        datafile.open(root + "models/lrgdr7fiducialmodel_matterpowerzFAR_new.dat");
        // Skip first line
        std::getline(datafile, line);
        // Read in rest of data file
        for(int i = 0; i < k_size; ++i)
        {
            std::getline(datafile, line);
            std::istringstream iss(line);
            double kdummy, plindummy, psmoothdummy, ratiodummy;
            iss >> kdummy >> plindummy >> psmoothdummy >> ratiodummy; 
            r_fid[2][i] = ratiodummy;
        }
        datafile.close();

        // Evaluate P_lin, P_nw, and r_nwhalofit as values of k_fid using spline in log space
        int n = P_lin[0].size();
        std::vector<double> lnP_lin(n);
        std::vector<double> lnP_nw(n);
        std::vector<double> lnkh(n);
        for(int i = 0; i < n; ++i)
            lnkh[i] = std::log(kh[i]);
        std::vector< std::vector<double> > P_lin_atfid(3, std::vector<double>(k_size));
        std::vector< std::vector<double> > P_nw_atfid(3, std::vector<double>(k_size));
        std::vector< std::vector<double> > r_nwhalofit_atfid(3, std::vector<double>(k_size));
        for(int i = 0; i < 3; ++i)
        {
            for(int j = 0; j < n; ++j)
            {
                lnP_lin[j] = std::log(P_lin[i][j]);
                lnP_nw[j] = std::log(P_nw[i][j]);
            }
            //Math::CubicSpline P_lin_spline(kh, P_lin[i]);
            //Math::CubicSpline P_nw_spline(kh, P_nw[i]);
            //Math::CubicSpline lnP_lin_spline(lnkh, lnP_lin);
            //Math::CubicSpline lnP_nw_spline(lnkh, lnP_nw);
            //Math::CubicSpline r_nwhalofit_spline(lnkh, r_nwhalofit[i]);
            Math::TableFunction<double, double> lnP_lin_spline;
            Math::TableFunction<double, double> lnP_nw_spline;
            Math::TableFunction<double, double> r_nwhalofit_spline;
            for(int j = 0; j < n; ++j)
            {
                lnP_lin_spline[lnkh[j]] = lnP_lin[j];
                lnP_nw_spline[lnkh[j]] = lnP_nw[j];
                r_nwhalofit_spline[lnkh[j]] = r_nwhalofit[i][j];
            }
            for(int j = 0; j < k_size; ++j)
            {
                //P_lin_atfid[i][j] = P_lin_spline.evaluate(kh_fid[j]);
                P_lin_atfid[i][j] = std::exp(lnP_lin_spline.evaluate(std::log(kh_fid[j])));
                //P_nw_atfid[i][j] = P_nw_spline.evaluate(kh_fid[j]);
                P_nw_atfid[i][j] = std::exp(lnP_nw_spline.evaluate(std::log(kh_fid[j])));
                r_nwhalofit_atfid[i][j] = r_nwhalofit_spline.evaluate(std::log(kh_fid[j]));
            }
        }

        // Debugging code
        //std::vector<double> getabstransferscale2(4);
        //// Read in P_lin, P_nw, and r_nwhalofit_atfid from file
        //datafile.open("/Users/ncanac/lrgdr7like/models/lrgdr7model_matterpowerzNEAR.dat");
        //// Skip first line
        ////std::getline(datafile, line);
        //// Read in rest of data file
        ////for(int i = 0; i < k_size; ++i)
        //for(int i = 0; i < k_size+1; ++i)
        //{
        //    std::getline(datafile, line);
        //    std::istringstream iss(line);
        //    double kdummy, plindummy, psmoothdummy, ratiodummy;
        //    iss >> kdummy >> plindummy >> psmoothdummy >> ratiodummy; 
        //    if(i == 0)
        //    {
        //        getabstransferscale2[0] = kdummy; 
        //    }
        //    else
        //    {
        //        P_lin_atfid[0][i-1] = plindummy;
        //        P_nw_atfid[0][i-1] = psmoothdummy;
        //        r_nwhalofit_atfid[0][i-1] = ratiodummy;
        //    }
        //}
        //datafile.close();
        //// Read in MID model
        //datafile.open("/Users/ncanac/lrgdr7like/models/lrgdr7model_matterpowerzMID.dat");
        //// Skip first line
        ////std::getline(datafile, line);
        //// Read in rest of data file
        ////for(int i = 0; i < k_size; ++i)
        //for(int i = 0; i < k_size+1; ++i)
        //{
        //    std::getline(datafile, line);
        //    std::istringstream iss(line);
        //    double kdummy, plindummy, psmoothdummy, ratiodummy;
        //    iss >> kdummy >> plindummy >> psmoothdummy >> ratiodummy; 
        //    if(i == 0)
        //    {
        //        getabstransferscale2[1] = kdummy; 
        //    }
        //    else
        //    {
        //        P_lin_atfid[1][i-1] = plindummy;
        //        P_nw_atfid[1][i-1] = psmoothdummy;
        //        r_nwhalofit_atfid[1][i-1] = ratiodummy;
        //    }
        //}
        //datafile.close();
        //// Read in FAR model
        //datafile.open("/Users/ncanac/lrgdr7like/models/lrgdr7model_matterpowerzFAR.dat");
        //// Skip first line
        ////std::getline(datafile, line);
        //// Read in rest of data file
        ////for(int i = 0; i < k_size; ++i)
        //for(int i = 0; i < k_size+1; ++i)
        //{
        //    std::getline(datafile, line);
        //    std::istringstream iss(line);
        //    double kdummy, plindummy, psmoothdummy, ratiodummy;
        //    iss >> kdummy >> plindummy >> psmoothdummy >> ratiodummy; 
        //    if(i == 0)
        //    {
        //        getabstransferscale2[2] = kdummy; 
        //    }
        //    else
        //    {
        //        P_lin_atfid[2][i-1] = plindummy;
        //        P_nw_atfid[2][i-1] = psmoothdummy;
        //        r_nwhalofit_atfid[2][i-1] = ratiodummy;
        //        //output_screen(plindummy << std::endl);
        //    }
        //}
        //datafile.close();
        // // Read in z=0 model
        //datafile.open("/Users/ncanac/lrgdr7like/models/lrgdr7model_matterpowerz0.dat");
        //std::getline(datafile, line);
        //std::istringstream iss(line);
        //double kdummy, plindummy, psmoothdummy, ratiodummy;
        //iss >> kdummy >> plindummy >> psmoothdummy >> ratiodummy; 
        //getabstransferscale2[3] = kdummy; 
        //datafile.close();
      
        //output P_lin, P_nw, and r_nwhalofit to files
        std::ofstream outfile("lrgPsNEAR.txt");
        for(int i = 0; i < k_size; ++i)
            outfile << kh_fid[i] << " " << P_lin_atfid[0][i] << " " << P_nw_atfid[0][i] << " " << r_nwhalofit_atfid[0][i] << std::endl;
        outfile.close();
        outfile.open("lrgPsMID.txt");
        for(int i = 0; i < k_size; ++i)
            outfile << kh_fid[i] << " " << P_lin_atfid[1][i] << " " << P_nw_atfid[1][i] << " " << r_nwhalofit_atfid[1][i] << std::endl;
        outfile.close();
        outfile.open("lrgPsFAR.txt");
        for(int i = 0; i < k_size; ++i)
            outfile << kh_fid[i] << " " << P_lin_atfid[2][i] << " " << P_nw_atfid[2][i] << " " << r_nwhalofit_atfid[2][i] << std::endl;
        outfile.close();

        // End debugging code

        // Set weights to do weighted sum for P_halo as in eq. 17 of BR09
        std::vector<double> zweight {0.395, 0.355, 0.250};
        // Set values of sigma
        std::vector<double> sigma2BAO {86.9988, 85.1374, 84.5958};
        std::vector<double> powerscaletoz0(3, 0.0);
        powerscaletoz0[0] = pow(getabstransferscale[3], 2.0)/pow(getabstransferscale[0], 2.0);
        powerscaletoz0[1] = pow(getabstransferscale[3], 2.0)/pow(getabstransferscale[1], 2.0);
        powerscaletoz0[2] = pow(getabstransferscale[3], 2.0)/pow(getabstransferscale[2], 2.0);

        // Debugging: TODO remove later
        //powerscaletoz0[0] = pow(getabstransferscale2[3], 2.0)/pow(getabstransferscale2[0], 2.0);
        //powerscaletoz0[1] = pow(getabstransferscale2[3], 2.0)/pow(getabstransferscale2[1], 2.0);
        //powerscaletoz0[2] = pow(getabstransferscale2[3], 2.0)/pow(getabstransferscale2[2], 2.0);
        ///////////////////////////////

        output_screen("getabstransferscalez0: " << getabstransferscale[3] << std::endl);
        output_screen("getabstransferscaleNEAR: " << getabstransferscale[0] << std::endl);
        output_screen("getabstransferscaleMID: " << getabstransferscale[1] << std::endl);
        output_screen("getabstransferscaleFAR: " << getabstransferscale[2] << std::endl);

        // Calculate P_halo for NEAR, MID, and FAR using eq. 10 from BR09
        P_halo.resize(k_size);
        for(int i = 0; i < k_size; ++i)
        {
            double khval = kh_fid[i];
            for(int j = 0; j < 3; ++j)
            {
                double expval = std::exp(-1.0*pow(khval, 2.0)*sigma2BAO[j]*0.5);
                double psmear = (P_lin_atfid[j][i])*expval + (P_nw_atfid[j][i])*(1.0-expval);
                psmear = psmear*powerscaletoz0[j];
                double nlrat = r_nwhalofit_atfid[j][i]/r_fid[j][i];
                std::vector<double> fidpolys(3, 0.0);
                LRGtoICsmooth_(khval, fidpolys);
                P_halo[i] = P_halo[i] + zweight[j]*psmear*nlrat*fidpolys[j];
            }
        }
        outfile.open("phalo.txt");
        for(int i = 0; i < k_size; ++i)
            outfile << kh_fid[i] << " " << P_halo[i] << std::endl;
        outfile.close();
    }
};
