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

    void setkMax(double kMax, bool iskh=false)
    {
        double currentkMax = pt_->k_max_for_pk;
        if(iskh)
            pt_->k_max_for_pk = std::max(currentkMax, kMax*params_->getH());
        else
            pt_->k_max_for_pk = std::max(currentkMax, kMax);
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

    bool getLRGHaloPs(std::string root, Math::TableFunction<double, double>* ps, double khmax)
    {
        output_screen("getLRGHaloPs Checkpoint 1" << std::endl);
        check(init_, "need to initialize first");
        check(pt_->has_pk_matter, "matter ps not requested");

        double h = params_->getH();

        const double zNEAR = 0.235, zMID = 0.342, zFAR = 0.421;
        //const double aNEAR = 0.809717, aMID = 0.745156, aFAR = 0.70373;

        // Get linear matter power spectrum for NEAR, MID, and FAR
        Math::TableFunction<double, double> P_lin_funcNEAR;
        Math::TableFunction<double, double> P_lin_funcMID;
        Math::TableFunction<double, double> P_lin_funcFAR;
        Math::TableFunction<double, double> P_lin_funcz0;
        getMatterPs(zNEAR, &P_lin_funcNEAR);
        getMatterPs(zMID, &P_lin_funcMID);
        getMatterPs(zFAR, &P_lin_funcFAR);
        getMatterPs(0.0, &P_lin_funcz0);

        output_screen("Checkpoint 2" << std::endl);
        // Get value of k max and size of arrays
        int n = 0;
        Math::TableFunction<double, double>::iterator point;
        point = P_lin_funcNEAR.begin();
        while(point->first <= khmax*h && point != P_lin_funcNEAR.end())
        {
            ++point;
            ++n;
        }
        ++n;
        output_screen("point->first: " << point->first << std::endl);
        check(point->first > khmax*h, "getMatterPs kmax too small");
        // Create some arrays to store the various power spectra and initialize them
        double kvals[n];
        double P_linNEAR[n];
        double P_linMID[n];
        double P_linFAR[n];
        double P_linz0[n];
        double lnP_linNEAR[n];
        double lnP_linMID[n];
        double lnP_linFAR[n];
        //int itemp = 0;
        //for(auto const &point : P_lin_funcNEAR)
        point = P_lin_funcNEAR.begin();
        kvals[0] = point->first;
        P_linNEAR[0] = point->second;
        lnP_linNEAR[0] = std::log(point->second);
        for(int i = 1; i < n; ++i)
        {
            ++point;
            kvals[i] = point->first;
            double Pval = point->second;
            P_linNEAR[i] = Pval;
            lnP_linNEAR[i] = std::log(Pval);
            //++itemp;
        }
        output_screen("max kval = " << kvals[n-1] << std::endl);
        //check(n == itemp, "number of k values mismatch");
        check(kvals[n-1] > khmax*h, "k values mismatch");
        //itemp = 0;
        //for(auto const &point : P_lin_funcMID)
        point = P_lin_funcMID.begin();
        P_linMID[0] = point->second;
        lnP_linMID[0] = std::log(point->second);
        for(int i = 1; i < n; ++i)
        {
            ++point;
            check(std::abs(kvals[i] - point->first) < 0.01*kvals[i], "kvals should be identical");
            double Pval = point->second;
            P_linMID[i] = Pval;
            lnP_linMID[i] = std::log(Pval);
            //++itemp;
        }
        //check(n == itemp, "number of k values mismatch");
        //itemp = 0;
        //for(auto const &point : P_lin_funcFAR)
        point = P_lin_funcFAR.begin();
        P_linFAR[0] = point->second;
        lnP_linFAR[0] = std::log(point->second);
        for(int i = 1; i < n; ++i)
        {
            ++point;
            check(std::abs(kvals[i] - point->first) < 0.01*kvals[i], "kvals should be identical");
            double Pval = point->second;
            P_linFAR[i] = Pval;
            lnP_linFAR[i] = std::log(Pval);
            //++itemp;
        }
        //check(n == itemp, "number of k values mismatch");
        //itemp = 0;
        //for(auto const &point : P_lin_funcz0)
        point = P_lin_funcz0.begin();
        P_linz0[0] = point->second;
        for(int i = 1; i < n; ++i)
        {
            ++point;
            check(std::abs(kvals[i] - point->first) < 0.01*kvals[i], "kvals should be identical");
            double Pval = point->second;
            P_linz0[i] = Pval;
            //++itemp;
        }
        //check(n == itemp, "number of k values mismatch");

        // Check that all arrays are initialized
        bool success = true;
        for(int i = 0; i < n; ++i)
        {
            if(P_linNEAR[i] <= 0 || lnP_linNEAR[i] <= 0 || P_linMID[i] <= 0 || lnP_linMID[i] <= 0 ||
                P_linFAR[i] <= 0 || lnP_linFAR[i] <= 0 || P_linz0[i] <= 0)
                success = false;
        }
        check(success, "P_lin arrays not initialized properly");

        // extract getabstransferscale for NEAR, MID, and FAR
        const double khmindata = 0.0221168;
        std::vector<double> getabstransferscale(4);
        double getabstransferscaleNEAR, getabstransferscaleMID, getabstransferscaleFAR;
        int itemp = 0;
        while(kvals[itemp]/h < khmindata && itemp < n)
            ++itemp;
        check(itemp < n, "khmindata > kvals");
        getabstransferscale[0] = sqrt(P_linNEAR[itemp]*pow(h, 3.0));        
        getabstransferscale[1] = sqrt(P_linMID[itemp]*pow(h, 3.0));        
        getabstransferscale[2] = sqrt(P_linFAR[itemp]*pow(h, 3.0));        
        getabstransferscale[3] = sqrt(P_linz0[itemp]*pow(h, 3.0));        
       
        // Initialize P_nw (linear power spectrum with "no wiggles") based on spline method used in BR09        
        double lnP_nwNEAR[n];
        double lnP_nwMID[n];
        double lnP_nwFAR[n];
        double P_nwNEAR[n];
        double P_nwMID[n];
        double P_nwFAR[n];
        dopksmoothbspline_(kvals, lnP_linNEAR, lnP_nwNEAR, n);
        dopksmoothbspline_(kvals, lnP_linMID, lnP_nwMID, n);
        dopksmoothbspline_(kvals, lnP_linFAR, lnP_nwFAR, n);
        for(int i = 0; i < n; ++i)
        {
            P_nwNEAR[i] = std::exp(lnP_nwNEAR[i]);
            P_nwMID[i] = std::exp(lnP_nwMID[i]);
            P_nwFAR[i] = std::exp(lnP_nwFAR[i]);
        }

        // Check that P_nw arrays are initialized
        for(int i = 0; i < n; ++i)
        {
            if(P_nwNEAR[i] <= 0 || lnP_nwNEAR[i] <= 0 || P_nwMID[i] <= 0 || lnP_nwMID[i] <= 0 ||
                P_nwFAR[i] <= 0 || lnP_nwFAR[i] <= 0)
                success = false;
        }
        check(success, "P_nw arrays not initialized properly");

        // Apply halofit model for nonlinear structure growth to P_nw to generate P_halofitnw
        // Tau is the conformal time at z
        // k_nl is the value of k where power spectrum becomes nonlinear
        double P_halofitnwNEAR[n];
        double P_halofitnwMID[n];
        double P_halofitnwFAR[n];
        double tau;
        double k_nl;
        // Apply halofit model to P_nw for NEAR
        StandardException exc;
        bool halofit_success = true;
        //nonlinear_k_nl_at_z(br_, nl_, zNEAR, &k_nl);
        if(background_tau_of_z(br_, zNEAR, &tau) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: background_tau_of_z failed!" << std::endl << br_->error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }
        if(nonlinear_halofit(pr_, br_, pm_, nl_, tau, P_nwNEAR, P_halofitnwNEAR, &k_nl) == _FAILURE_)
        {
            halofit_success = false;
            for(int i = 0; i < n; ++i)
                P_halofitnwNEAR[i] = P_nwNEAR[i];
        }
        // Apply halofit model to P_nw for MID 
        //nonlinear_k_nl_at_z(br_, nl_, zMID, &k_nl);
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

        // Check that P_nw arrays are initialized
        for(int i = 0; i < n; ++i)
        {
            if(P_halofitnwNEAR[i] <= 0 || P_halofitnwMID[i] <= 0 || P_halofitnwFAR[i] <= 0)
                success = false;
        }
        check(success, "P_halofitnw arrays not initialized properly");

        // Make sure that P_halofitnw is positive everywhere
        //for(int i = 0; i < n; ++i)
        //{
        //    if(P_halofitnwNEAR[i] <= 0 || P_halofitnwMID[i] <= 0 || P_halofitnwFAR[i] <= 0)
        //    {
        //        output_screen("BADHALOFIT" << std::endl);
        //        return false;
        //    }
        //}
        if(halofit_success)
            output_screen("BADHALOFIT" << std::endl);
        output_screen("GOODHALOFIT" << std::endl);

        // Calculate factor r_halofit, the ratio of P_halofitnw to P_nw
        std::vector<double> r_halofitNEAR(n);
        std::vector<double> r_halofitMID(n);
        std::vector<double> r_halofitFAR(n);
        for(int i = 0; i < n; ++i)
        {
            r_halofitNEAR[i] = P_halofitnwNEAR[i]/P_nwNEAR[i];
            r_halofitMID[i] = P_halofitnwMID[i]/P_nwMID[i];
            r_halofitFAR[i] = P_halofitnwFAR[i]/P_nwFAR[i];
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

        for(int i = 0; i < P_halo.size(); ++i)
        {
            if(kh_fid[i] <= 0 || P_halo[i] <= 0)
                success = false; 
        }
        check(success, "P_halo initialization error" << std::endl);

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
    // I think input k is really k/h
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
        std::ifstream datafile(root + "models/lrgdr7fiducialmodel_matterpowerzNEAR.dat");
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
        check(kh_fid[k_size-1] < kh[kh.size()-1], "kh_fid > kh" << std::endl);
        // Read in MID model
        datafile.open(root + "models/lrgdr7fiducialmodel_matterpowerzMID.dat");
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
        datafile.open(root + "models/lrgdr7fiducialmodel_matterpowerzFAR.dat");
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

        // Evaluate P_lin, P_nw, and r_nwhalofit at values of k_fid using cubic spline in log space
        // Note: Doing this in linear space doesn't seem to make a difference
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
                double lnkhj = std::log(kh_fid[j]);
                check(lnkhj > lnkh[0] && lnkhj < lnkh[n-1], "lnkh: (" << lnkh[0] << ", " << lnkh[n-1] << "), requested value: " << lnkhj);
                P_lin_atfid[i][j] = std::exp(lnP_lin_spline.evaluate(lnkhj));
                //P_lin_atfid[i][j] = std::exp(lnP_lin_spline.evaluate(std::log(kh_fid[j])));
                //P_nw_atfid[i][j] = P_nw_spline.evaluate(kh_fid[j]);
                P_nw_atfid[i][j] = std::exp(lnP_nw_spline.evaluate(lnkhj));
                r_nwhalofit_atfid[i][j] = r_nwhalofit_spline.evaluate(lnkhj);
                //P_nw_atfid[i][j] = std::exp(lnP_nw_spline.evaluate(std::log(kh_fid[j])));
                //r_nwhalofit_atfid[i][j] = r_nwhalofit_spline.evaluate(std::log(kh_fid[j]));
            }
        }

        // Set weights to do weighted sum for P_halo as in eq. 17 of BR09
        std::vector<double> zweight {0.395, 0.355, 0.250};
        // Set values of sigma
        std::vector<double> sigma2BAO {86.9988, 85.1374, 84.5958};
        std::vector<double> powerscaletoz0(3, 0.0);
        powerscaletoz0[0] = pow(getabstransferscale[3], 2.0)/pow(getabstransferscale[0], 2.0);
        powerscaletoz0[1] = pow(getabstransferscale[3], 2.0)/pow(getabstransferscale[1], 2.0);
        powerscaletoz0[2] = pow(getabstransferscale[3], 2.0)/pow(getabstransferscale[2], 2.0);

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
    }
};
