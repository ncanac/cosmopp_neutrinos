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
        Math::TableFunction<double, double> P_lin_func;
        getMatterPs(z, &P_lin_func);

        // Create some arrays to store the various power spectra
        int n = P_lin_func.size();
        double kvals[n];
        double P_lin[n];
        double lnP_lin[n];
        //for(std::map<double, double>::iterator ps_it = P_lin_func.begin(); ps_it != P_lin_func.end(); ++ps_it)
        int i = 0;
        for(auto const &point : P_lin_func)
        {
            kvals[i] = point.first;
            P_lin[i] = point.second;
            lnP_lin[i] = std::log(point.second);
            ++i;
        }
        
        // Initialize P_nw (linear power spectrum with "no wiggles") based on spline method used in BR09        
        double lnP_nw[n];
        double P_nw[n];
        dopksmoothbspline_(kvals, lnP_lin, lnP_nw, n);
        for(int i = 0; i < n; ++i)
            P_nw[i] = std::exp(lnP_nw[i]);

        // Calculate P_damp for NEAR, MID, and FAR using eq. 10 from BR09
        std::vector< std::vector<double> > P_damp(3, std::vector<double>(n));
        std::vector<double> sigma2BAO {86.9988, 85.1374, 84.5958};
        std::vector<double> powerscaletoz0(3);
        double getabstransferscalez0 = 162.71; // TODO: Figure out what these are
        double getabstransferscaleNEAR = 144.31;
        double getabstransferscaleMID = 136.64;
        double getabstransferscaleFAR = 131.3;
        powerscaletoz0[0] = pow(getabstransferscalez0, 2.0)/pow(getabstransferscaleNEAR, 2.0);
        powerscaletoz0[1] = pow(getabstransferscalez0, 2.0)/pow(getabstransferscaleMID, 2.0);
        powerscaletoz0[2] = pow(getabstransferscalez0, 2.0)/pow(getabstransferscaleFAR, 2.0);
        for(int i = 0; i < 3; ++i)
        {
            for(int j = 0; j < n; ++j)
            {
                double expval = std::exp(-1.0 * pow(kvals[j], 2) * sigma2BAO[i] * 0.5);
                P_damp[i][j] = P_lin[j] * expval + P_nw[j] * (1.0 - expval);
                P_damp[i][j] = P_damp[i][j] * powerscaletoz0[i];
            }
        }
        output_screen("Checkpoint 1" << std::endl);
        // Initalize a cubic spline for each P_damp to be used for calculating P_halo later
        //std::vector<double> kvec(kvals, kvals + sizeof(kvals) / sizeof(double));
        //std::vector<Math::CubicSpline*> P_damp_spline(3);
        //for(int i = 0; i < 3; ++i)
        //{
        //    Math::CubicSpline tempSpline(kvec, P_damp[i]);
        //    P_damp_spline[i] = &tempSpline;
        //}
        std::vector<double> kvec(kvals, kvals + sizeof(kvals) / sizeof(double));
        Math::CubicSpline P_damp_splineNEAR(kvec, P_damp[0]);
        Math::CubicSpline P_damp_splineMID(kvec, P_damp[1]);
        Math::CubicSpline P_damp_splineFAR(kvec, P_damp[2]);
        
        output_screen("Checkpoint 2" << std::endl);

        // Apply halofit model for nonlinear structure growth to P_nw to generate P_halofitnw
        // Tau is the conformal time at z
        // k_nl is the value of k where power spectrum becomes nonlinear
        double P_halofitnw[n];
        double tau;
        double k_nl;
        nonlinear_k_nl_at_z(br_, nl_, z, &k_nl);
        background_tau_of_z(br_, z, &tau);
        nonlinear_halofit(pr_, br_, pm_, nl_, tau, P_nw, P_halofitnw, &k_nl);

        // Calculate factor r_halofit, the ratio of P_halofitnw to P_nw
        // Then P_DMhalofit = P_damp * r_halofit
        std::vector<double> r_halofit(n);
        std::vector<double> P_DMhalofit(n);
        double temp;
        for(int i = 0; i < n; ++i)
        {
            temp = P_halofitnw[i]/P_nw[i];
            r_halofit[i] = temp;
            //P_DMhalofit[i] = P_damp[i] * temp;
        }
        // Create a cubic spline for r_halofit for use later
        Math::CubicSpline r_halofit_spline(kvec, r_halofit);
        output_screen("Checkpoint 3" << std::endl);

        // Calculate r_DMdamp, model for ratio of nonlinear matter power spectrum
        // to damped linear power spectrum.
        // r_DMdamp = (r_halofit/r_halofitfid) * (P_DMfid/P_dampfid)
        // TODO: What is r_halofitfid, P_DMfid, and P_dampfid?

        // Now include halo bias
        // Calculate r_haloDMfid (P_halofid/P_damp) / (P_DMfid/P_damp)
        // P_halo = P_damp*r_DMdamp * r_haloDMfid * F_nuis
        // In BR09 likelihood code P_halo = psmear*nlrat*fidpolys
        // where psmear = P_damp, I think?
        // nlrat = outpowerrationwhalofit/ratio_power_nw_nl_fid **
        // fidpoly = hard coded polynomial fits to near, mid, and far subsamples calculated in LRGtoICsmooth
        // **
        // ratio_power_nw_nl_fid is rationwhalofit from lrgdr7fiducialmodel_matterpowerz*.dat
        // outpowerrationwhalofit is rationwhalofit from lrgdr7model_matterpowerz*.dat, maybe this is just r_halofit?
        // TODO: I don't think the above is true, but using it for now...
        std::vector<double> k_fid; // k values should be identical so only need one
        std::vector< std::vector<double> > r_fid(3); // two dimensional vector to hold ratio of nw to nl for NEAR, MID, and FAR
        std::string root = "/Volumes/Data1/ncanac/cosmopp_neutrinos/data/LRGDR7/";
        // Read in NEAR model
        std::ifstream datafile(root + "models/lrgdr7fiducialmodel_matterpowerzNEAR.dat");
        // Skip first line
        std::string line;
        std::getline(datafile, line);
        // Read in rest of data file
        //std::vector<double> temp_vec;
        while(!datafile.eof())
        {
            std::getline(datafile, line);
            std::istringstream iss(line);
            double kdummy, plindummy, psmoothdummy, ratiodummy;
            iss >> kdummy >> plindummy >> psmoothdummy >> ratiodummy; 
            k_fid.push_back(kdummy); // Only need to do this once
            //Plin_fid[i] = plindummy;
            //Psmooth_fid[i] = psmoothdummy;
            //temp_vec.push_back(ratiodummy);
            r_fid[0].push_back(ratiodummy);
        }
        //r_fid.push_back(temp_vec);
        datafile.close();
        // Read in MID model
        datafile.open(root + "models/lrgdr7fiducialmodel_matterpowerzMID.dat");
        // Skip first line
        //std::string line;
        std::getline(datafile, line);
        // Read in rest of data file
        //temp_vec.clear();
        while(!datafile.eof())
        {
            std::getline(datafile, line);
            std::istringstream iss(line);
            double kdummy, plindummy, psmoothdummy, ratiodummy;
            iss >> kdummy >> plindummy >> psmoothdummy >> ratiodummy; 
            //k_fid.push_back(kdummy);
            //Plin_fid[i] = plindummy;
            //Psmooth_fid[i] = psmoothdummy;
            //temp_vec.push_back(ratiodummy);
            r_fid[1].push_back(ratiodummy);
        }
        //r_fid.push_back(temp_vec);
        datafile.close();
        // Read in FAR model
        datafile.open(root + "models/lrgdr7fiducialmodel_matterpowerzFAR.dat");
        // Skip first line
        //std::string line;
        std::getline(datafile, line);
        // Read in rest of data file
        //temp_vec.clear()
        while(!datafile.eof())
        {
            std::getline(datafile, line);
            std::istringstream iss(line);
            double kdummy, plindummy, psmoothdummy, ratiodummy;
            iss >> kdummy >> plindummy >> psmoothdummy >> ratiodummy; 
            //k_fid.push_back(kdummy);
            //Plin_fid[i] = plindummy;
            //Psmooth_fid[i] = psmoothdummy;
            r_fid[2].push_back(ratiodummy);
        }
        //r_fid.push_back(temp_vec);
        datafile.close();

        // Hard coded for testing purposes. TODO: delete later
        // Read in data from model files for r_nwhalofit
        std::vector< std::vector<double> > r_nwhalofit(3);
        datafile.open(root + "models/lrgdr7model_matterpowerzNEAR.dat");
        // Skip first line
        //std::string line;
        std::getline(datafile, line);
        // Read in rest of data file
        //temp_vec.clear()
        while(!datafile.eof())
        {
            std::getline(datafile, line);
            std::istringstream iss(line);
            double kdummy, plindummy, psmoothdummy, ratiodummy;
            iss >> kdummy >> plindummy >> psmoothdummy >> ratiodummy; 
            //k_fid.push_back(kdummy);
            //Plin_fid[i] = plindummy;
            //Psmooth_fid[i] = psmoothdummy;
            r_nwhalofit[0].push_back(ratiodummy);
        }
        //r_fid.push_back(temp_vec);
        datafile.close();
        datafile.open(root + "models/lrgdr7model_matterpowerzMID.dat");
        // Skip first line
        //std::string line;
        std::getline(datafile, line);
        // Read in rest of data file
        //temp_vec.clear()
        while(!datafile.eof())
        {
            std::getline(datafile, line);
            std::istringstream iss(line);
            double kdummy, plindummy, psmoothdummy, ratiodummy;
            iss >> kdummy >> plindummy >> psmoothdummy >> ratiodummy; 
            //k_fid.push_back(kdummy);
            //Plin_fid[i] = plindummy;
            //Psmooth_fid[i] = psmoothdummy;
            r_nwhalofit[1].push_back(ratiodummy);
        }
        //r_fid.push_back(temp_vec);
        datafile.close();
        datafile.open(root + "models/lrgdr7model_matterpowerzFAR.dat");
        // Skip first line
        //std::string line;
        std::getline(datafile, line);
        // Read in rest of data file
        //temp_vec.clear()
        while(!datafile.eof())
        {
            std::getline(datafile, line);
            std::istringstream iss(line);
            double kdummy, plindummy, psmoothdummy, ratiodummy;
            iss >> kdummy >> plindummy >> psmoothdummy >> ratiodummy; 
            //k_fid.push_back(kdummy);
            //Plin_fid[i] = plindummy;
            //Psmooth_fid[i] = psmoothdummy;
            r_nwhalofit[2].push_back(ratiodummy);
        }
        //r_fid.push_back(temp_vec);
        datafile.close();

        output_screen("Checkpoint 4" << std::endl);
        
        // Set number of k values in fiducial model
        int k_size = k_fid.size();

        // Set weights to do weighted sum for P_halo as in eq. 17 of BR09
        std::vector<double> zweight {0.395, 0.355, 0.250};

        // Now calculate P_halo as P_damp * r_fid * fid_polys
        std::vector<double> P_halo(k_size, 0.0);
        //std::ofstream out("k_Pdamp_rfid_fidpoly.txt");
        //for(int i = 0; i < k_size; ++i)
        //{
        //    for(int j = 0; j < 3; ++j)
        //    {
        //        std::vector<double> fidpolys;
        //        LRGtoICsmooth(k_fid[i], fidpolys);
        //        //out << k_fid[i] << " " << P_damp_spline[j]->evaluate(k_fid[i]) << " " << r_fid[j][i] << " " << fidpolys[j] << std::endl; 
        //        P_halo[i] += zweight[j] * (P_damp_spline[j])->evaluate(k_fid[i]) * r_fid[j][i] * fidpolys[j];
        //        //P_halo[i] += zweight[j] * r_fid[j][i] * fidpolys[j];
        //    }
        //}
        //out.close();
        std::ofstream out("r_halofit.txt");
        for(int i = 0; i < k_size; ++i)
        {
            out << k_fid[i] << " " << r_halofit_spline.evaluate(k_fid[i]) << std::endl;
        }
        out.close();
        for(int i = 0; i < k_size; ++i)
        {
            std::vector<double> fidpolys;
            LRGtoICsmooth(k_fid[i], fidpolys);
            //double nlratNEAR = r_halofit_spline.evaluate(k_fid[i]) / r_fid[0][i];
            //double nlratMID = r_halofit_spline.evaluate(k_fid[i]) / r_fid[1][i];
            //double nlratFAR = r_halofit_spline.evaluate(k_fid[i]) / r_fid[2][i];
            double nlratNEAR = r_nwhalofit[0][i] / r_fid[0][i];
            double nlratMID = r_nwhalofit[1][i] / r_fid[1][i];
            double nlratFAR = r_nwhalofit[2][i] / r_fid[2][i];
            //P_halo[i] = zweight[0] * P_damp_splineNEAR.evaluate(k_fid[i]) * r_fid[0][i] * fidpolys[0]
            //            + zweight[1] * P_damp_splineMID.evaluate(k_fid[i]) * r_fid[1][i] * fidpolys[1]
            //            + zweight[2] * P_damp_splineFAR.evaluate(k_fid[i]) * r_fid[2][i] * fidpolys[2];
            P_halo[i] = zweight[0] * P_damp_splineNEAR.evaluate(k_fid[i]) * nlratNEAR * fidpolys[0]
                        + zweight[1] * P_damp_splineMID.evaluate(k_fid[i]) * nlratMID * fidpolys[1]
                        + zweight[2] * P_damp_splineFAR.evaluate(k_fid[i]) * nlratFAR * fidpolys[2];
            //P_halo[i] = P_damp_splineNEAR.evaluate(k_fid[i]) * r_fid[0][i] * fidpolys[0];
        }

        output_screen("Checkpoint 5" << std::endl);

        // Generate linear and nonlinear power spectra from Class for comparison
        //int ksize = sp_->ln_k_size;
        //double pk_nl[100000];
        //spectra_pk_nl_at_z(br_, sp_, linear, z, pk_nl);
        //std::ofstream out("pk_nl_class.txt");
        //for(int i = 0; i < nl_->k_size; ++i)
        //    out << std::exp(sp_->ln_k[i]) << " " << " " << pk_nl[i] << std::endl; 
        //out.close();

        // Store LRG power spectrum in table function
        ps->clear();
        for(i = 0; i < k_size; ++i)
            (*ps)[k_fid[i]] = P_halo[i];
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

    // Copied form BR09 likelihood code
    // Hard coding of polynomial fits to near, mid, and far subsamples.
    void LRGtoICsmooth(double k, std::vector<double>& fidpolys)
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
            fidMID = (1.0 - 1.97873*k + 20.8551*k*k - 50.0376*pow(k, 3) + 36.4056*(k, 4)) * 1.04384;

        if(k < 0.19148)
            fidFAR = (1.0 - 0.475028*k + 6.69004*k*k);
        else
            fidFAR = (1.0 - 1.84891*k + 21.3479*k*k - 52.4846*pow(k, 3) + 38.9541*pow(k, 4)) * 1.03753;

        fidpolys[0] = fidNEAR;
        fidpolys[1] = fidMID;
        fidpolys[2] = fidFAR;
    }
};
