#include <fstream>
#include <cmath>

#include <cmb.hpp>
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

    double getrsdrag_fudge()
    {
        check(init_, "need to initialize first");
        return th_->rs_d_fudge;
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

    double getHubble(double z)
    {
        double tau;
        int last_index;
        double* pvecback;

        pvecback = (double*) calloc(br_->bg_size,sizeof(double));
        background_tau_of_z(br_, z, &tau);
        background_at_tau(br_, tau, br_->long_info, br_->inter_normal, &last_index, pvecback);
        
        free(pvecback);
        
        return pvecback[br_->index_bg_H];
    }

    void getLRGPs(double z, Math::TableFunction<double, double>* ps)
    {
        check(init_, "need to initialize first");
        check(pt_->has_pk_matter, "matter ps not requested");
        check(z >= 0 && z <= sp_->z_max_pk, "invalid z = " << z);

        // Get linear matter power spectrum
        Math::TableFunction<double, double> P_lin;
        getMatterPs(z, &P_lin);

        // Create some temporary arrays to pass to dopksmoothbspline
        int n = P_lin.size();
        double kvals[n];
        double lnP_lin[n];
        double lnP_nw[n];
        //for(std::map<double, double>::iterator ps_it = P_lin.begin(); ps_it != P_lin.end(); ++ps_it)
        int i = 0;
        for(auto const &point : P_lin)
        {
            kvals[i] = point.first;
            lnP_lin[i] = std::log(point.second);
            ++i;
        }
        
        dopksmoothbspline_(kvals, lnP_lin, lnP_nw, n);

        ps->clear();

        for(i = 0; i < n; ++i)
            (*ps)[kvals[i]] = std::exp(lnP_nw[i]);
    }

    double getR_NL()
    {
        double z, k_nl, k;
        double* pvecback;
        int junk;
        double r_nl;
        int index_tau; // Think I can get rid of this and replace last instance with nl_->tau[index_tau]

        output_screen("Checkpoint 1" << std::endl);

        pt_->k_max_for_pk = kMax_;
        nl_->method = nl_halofit;
        pr_->halofit_dz = 0.1;
        pr_->halofit_min_k_nonlinear = 0.0035;
        pr_->halofit_sigma_precision = 0.05;

        nonlinear_free(nl_);
        nonlinear_init(pr_, br_, th_, pt_, pm_, nl_);

        //nl_->method = nl_halofit;
        output_screen("Checkpoint 1.1" << std::endl);
        z = 0.0;
        nonlinear_k_nl_at_z(br_, nl_, z, &k_nl);

        output_screen("Checkpoint 2" << std::endl);

        pvecback = (double*) calloc(br_->bg_size,sizeof(double));

        for(index_tau = 0; index_tau < nl_->tau_size; ++index_tau)
        {
            background_at_tau(br_, nl_->tau[index_tau], br_->short_info,
                                br_->inter_normal, &junk, pvecback);
        }

        output_screen("Checkpoint 3" << std::endl);

        z = br_->a_today/pvecback[br_->index_bg_a]-1.0;
        
        std::ofstream out("test_nonlinear.txt");
        for(int index_k = 0; index_k < nl_->k_size; ++index_k)
        {
            k = nl_->k[index_k];
            r_nl = nl_->nl_corr_density[index_tau * nl_->k_size + index_k];
            out << z << " " << k << " " << r_nl << std::endl;
        }
        out.close();
        free(pvecback);
        return 0;
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
};
