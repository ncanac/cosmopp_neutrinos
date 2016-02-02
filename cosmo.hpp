#pragma once
#include <fstream>
#include <sstream>
#include <cmath>
#include <cfloat>

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

        free(pvecback);

        return pvecback[br_->index_bg_ang_distance];
    }

    //double getSigma8()
    //{
    //    check(init_, "need to initialize first");
    //    return sp_->sigma8;
    //}

    double getHubble(double z)
    {
        double tau;
        int last_index;
        double* pvecback;

        pvecback = (double*) calloc(br_->bg_size, sizeof(double));
        background_tau_of_z(br_, z, &tau);
        background_at_tau(br_, tau, br_->long_info, br_->inter_normal, &last_index, pvecback);
        
        free(pvecback);
        
        return pvecback[br_->index_bg_H];
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

        return out;
    }

    bool getLRGHaloPs(std::string root, Math::TableFunction<double, double>* ps)
    {
        check(init_, "need to initialize first");
        check(pt_->has_pk_matter, "matter ps not requested");

        double h = params_->getH();
        
        //std::vector<double> v;
        //params_->getAllParameters(v);

        //for(int i = 0; i < v.size(); ++i)
        //    output_screen(v[i] << std::endl);

        const double zNEAR = 0.235, zMID = 0.342, zFAR = 0.421;
        //const double aNEAR = 0.809717, aMID = 0.745156, aFAR = 0.70373;

        // Get linear matter power spectrum for NEAR, MID, and FAR
        //Math::TableFunction<double, double> P_lrg;
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
            check(std::isnan(Pval) == 0, "NEAR power is nan");
            //if(itemp == 0)
            //    output_screen("P_lin_NEAR[0]: " << Pval << std::endl);
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
            check(std::isnan(Pval) == 0, "MID power is nan");
            //if(itemp == 0)
            //    output_screen("P_lin_MID[0]: " << Pval << std::endl);
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
            check(std::isnan(Pval) == 0, "FAR power is nan");
            //if(itemp == 0)
            //    output_screen("P_lin_FAR[0]: " << Pval << std::endl);
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
            check(std::isnan(Pval) == 0, "z0 power is nan");
            //if(itemp == 0)
            //    output_screen("P_lin_z0[0]: " << Pval << std::endl);
            P_linz0[itemp] = Pval;
            ++itemp;
        }
        check(n == itemp, "number of k values mismatch");

        // extract getabstransferscale for NEAR, MID, and FAR
        const double khmindata = 0.0221168;
        std::vector<double> getabstransferscale(4);
        //double getabstransferscaleNEAR, getabstransferscaleMID, getabstransferscaleFAR;
        itemp = 0;
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
            check(std::isnan(lnP_nwNEAR[i]) == 0 || std::isnan(lnP_nwMID[i]) == 0 || std::isnan(lnP_nwFAR[i]) == 0, "No wiggles power is nan");
            //if(i == 0)
            //{
            //    output_screen("lnP_nwNEAR[0]: " << lnP_nwNEAR[i] << std::endl);
            //    output_screen("lnP_nwMID[0]: " << lnP_nwMID[i] << std::endl);
            //    output_screen("lnP_nwFAR[0]: " << lnP_nwFAR[i] << std::endl);
            //}
            P_nwNEAR[i] = std::exp(lnP_nwNEAR[i]);
            P_nwMID[i] = std::exp(lnP_nwMID[i]);
            P_nwFAR[i] = std::exp(lnP_nwFAR[i]);
        }

        //output_screen("P_nwNEAR: " << P_nwNEAR[0] << std::endl);
        //output_screen("P_nwMID: " << P_nwMID[0] << std::endl);
        //output_screen("P_nwFAR: " << P_nwFAR[0] << std::endl);

        // Apply halofit model for nonlinear structure growth to P_nw to generate P_halofitnw
        // Tau is the conformal time at z
        // k_nl is the value of k where power spectrum becomes nonlinear
        double P_halofitnwNEAR[n];
        double P_halofitnwMID[n];
        double P_halofitnwFAR[n];
        double tau;
        double k_nl;
        // Apply halofit model to P_nw for NEAR
        nonlinear_k_nl_at_z(br_, nl_, zNEAR, &k_nl);
        //if(k_nl < pr_->halofit_min_k_nonlinear)
        //    k_nl = 0.3; 
        background_tau_of_z(br_, zNEAR, &tau);
        //output_screen("NEAR tau, k_nl: " << tau << ", " << k_nl << std::endl);
        nonlinear_halofit(pr_, br_, pm_, nl_, tau, P_nwNEAR, P_halofitnwNEAR, &k_nl);
        // Apply halofit model to P_nw for MID 
        nonlinear_k_nl_at_z(br_, nl_, zMID, &k_nl);
        background_tau_of_z(br_, zMID, &tau);
        //output_screen("MID tau, k_nl: " << tau << ", " << k_nl << std::endl);
        nonlinear_halofit(pr_, br_, pm_, nl_, tau, P_nwMID, P_halofitnwMID, &k_nl);
        // Apply halofit model to P_nw for FAR 
        nonlinear_k_nl_at_z(br_, nl_, zFAR, &k_nl);
        background_tau_of_z(br_, zFAR, &tau);
        //output_screen("FAR tau, k_nl: " << tau << ", " << k_nl << std::endl);
        nonlinear_halofit(pr_, br_, pm_, nl_, tau, P_nwFAR, P_halofitnwFAR, &k_nl);

        // Make sure that P_halofitnw is positive everywhere
        for(int i = 0; i < n; ++i)
        {
            if(P_halofitnwNEAR[i] <= 0 || P_halofitnwMID[i] <= 0 || P_halofitnwFAR[i] <= 0)
            {
                output_screen("BADHALOFIT" << std::endl);
                return false;
            }
        }
        //double sumNEAR = 0;
        //double sumMID = 0;
        //double sumFAR = 0;
        //for(int i = 0; i < n; ++i)
        //{
        //    sumNEAR += P_halofitnwNEAR[i];
        //    sumMID += P_halofitnwMID[i];
        //    sumFAR += P_halofitnwFAR[i];
        //}
        //if(sumNEAR < 1 || sumMID < 1 || sumFAR < 1)
        //{
        //    output_screen("BADHALOFIT" << std::endl);
        //    return false;
        //}
        output_screen("GOODHALOFIT" << std::endl);

        // Calculate factor r_halofit, the ratio of P_halofitnw to P_nw
        std::vector<double> r_halofitNEAR(n);
        std::vector<double> r_halofitMID(n);
        std::vector<double> r_halofitFAR(n);
        for(int i = 0; i < n; ++i)
        {
            //if(i == 0)
            //{
            //    output_screen("P_halofitnwNEAR, P_nwNEAR: " << P_halofitnwNEAR[i] << ", " << P_nwNEAR[i] << std::endl);
            //    output_screen("P_halofitnwMID, P_nwMID: " << P_halofitnwMID[i] << ", " << P_nwMID[i] << std::endl);
            //    output_screen("P_halofitnwFAR, P_nwFAR: " << P_halofitnwFAR[i] << ", " << P_nwFAR[i] << std::endl);
            //}
            r_halofitNEAR[i] = P_halofitnwNEAR[i]/P_nwNEAR[i];
            r_halofitMID[i] = P_halofitnwMID[i]/P_nwMID[i];
            r_halofitFAR[i] = P_halofitnwFAR[i]/P_nwFAR[i];
            check(std::isnan(r_halofitNEAR[i]) == 0 || std::isnan(r_halofitMID[i]) == 0 || std::isnan(r_halofitFAR[i]) == 0, "halofit power is nan");
            //if(i == 0)
            //{
            //    output_screen("r_halofitNEAR[0]: " << r_halofitNEAR[i] << std::endl);
            //    output_screen("r_halofitMID[0]: " << r_halofitMID[i] << std::endl);
            //    output_screen("r_halofitFAR[0]: " << r_halofitFAR[i] << std::endl);
            //}
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
        //output_screen(kh[0] << " " << P_lin[0][0] << " " << P_nw[0][0] << " " << r_nwhalofit[0][0] << std::endl);
        //output_screen(kh[0] << " " << P_lin[1][0] << " " << P_nw[1][0] << " " << r_nwhalofit[1][0] << std::endl);
        //output_screen(kh[0] << " " << P_lin[2][0] << " " << P_nw[2][0] << " " << r_nwhalofit[2][0] << std::endl);
        //output_screen(getabstransferscale[0] << " " << getabstransferscale[1] << " " << getabstransferscale[2] << " " << getabstransferscale[3] << std::endl);
        LRGTheory_(root, kh, P_lin, P_nw, r_nwhalofit, getabstransferscale, kh_fid, P_halo);
        //output_screen("kh, P_halo: " << kh_fid[0] << ", " << P_halo[0] << std::endl);

        // Store LRG power spectrum in table function
        ps->clear();
        for(int i = 0; i < kh_fid.size(); ++i)
        {
            check(!std::isnan(P_halo[i]), "P_halo is nan");
            //if(i == 0)
                //output_screen("P_halo[0]: " << P_halo[i] << std::endl);
            const double k = kh_fid[i];
            (*ps)[k] = P_halo[i];
        }

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
        //check(root + "models/lrgdr7fiducialmodel_matterpowerzNEAR.dat" == "/Volumes/Data1/ncanac/cosmopp_neutrinos/data/LRGDR7/models/lrgdr7fiducialmodel_matterpowerzNEAR.dat", "NEAR fiducial model file name error: " << root + "models/lrgdr7fiducialmodel_matterpowerzNEAR.dat");
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
        check(kh_fid[0] == 0.99900E-04 && kh_fid[k_size-1] == 0.78557, "NEAR fiducial model read error!");
        datafile.close();
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
        //output_screen("lnkh: (" << lnkh[0] << ", " << lnkh[n-1] << ")" << std::endl);
        //output_screen("lnkhfid: (" << std::log(kh_fid[0]) << ", " << std::log(kh_fid[k_size-1]) << ")" << std::endl);
        std::vector< std::vector<double> > P_lin_atfid(3, std::vector<double>(k_size));
        std::vector< std::vector<double> > P_nw_atfid(3, std::vector<double>(k_size));
        std::vector< std::vector<double> > r_nwhalofit_atfid(3, std::vector<double>(k_size));
        //output_screen("r_nwhalofit[0][0]: " << r_nwhalofit[0][0] << std::endl);
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
                //P_nw_atfid[i][j] = P_nw_spline.evaluate(kh_fid[j]);
                P_nw_atfid[i][j] = std::exp(lnP_nw_spline.evaluate(lnkhj));
                r_nwhalofit_atfid[i][j] = r_nwhalofit_spline.evaluate(lnkhj);
            }
        }

        //output_screen("P_lin_atfid/P_nw_atfid/r_nwhalofit_atfid:" << std::endl);
        //output_screen(P_lin_atfid[0][0] << " " << P_nw_atfid[0][0] << " " << r_nwhalofit_atfid[0][0] << std::endl);
        //output_screen(P_lin_atfid[1][0] << " " << P_nw_atfid[1][0] << " " << r_nwhalofit_atfid[1][0] << std::endl);
        //output_screen(P_lin_atfid[2][0] << " " << P_nw_atfid[2][0] << " " << r_nwhalofit_atfid[2][0] << std::endl);

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
        //output_screen("P_halo_init: " << P_halo[0] << std::endl);
        for(int i = 0; i < k_size; ++i)
        {
            double khval = kh_fid[i];
            for(int j = 0; j < 3; ++j)
            {
                double expval = std::exp(-1.0*pow(khval, 2.0)*sigma2BAO[j]*0.5);
                double psmear = (P_lin_atfid[j][i])*expval + (P_nw_atfid[j][i])*(1.0-expval);
                psmear = psmear*powerscaletoz0[j];
                //if(i == 0)
                //    output_screen("r_nwhalofit_atfid/r_fid: " << r_nwhalofit_atfid[j][i] << " " << r_fid[j][i] << std::endl);
                double nlrat = r_nwhalofit_atfid[j][i]/r_fid[j][i];
                //if(i == 0)
                //    output_screen("psmear/nlrat: " << psmear << " " << nlrat << std::endl);
                std::vector<double> fidpolys(3, 0.0);
                LRGtoICsmooth_(khval, fidpolys);
                //if(i == 0)
                //    output_screen("fidpolys: " << fidpolys[0] << " " << fidpolys[1] << " " << fidpolys[2] << std::endl);
                P_halo[i] = P_halo[i] + zweight[j]*psmear*nlrat*fidpolys[j];
            }
        }
    }
};
