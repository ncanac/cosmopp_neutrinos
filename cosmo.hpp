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

    void getLRGHaloPs(Math::TableFunction<double, double>* ps)
    {
        check(init_, "need to initialize first");
        check(pt_->has_pk_matter, "matter ps not requested");

        const double zNEAR = 0.235, zMID = 0.342, zFAR = 0.421;
        const double aNEAR = 0.809717, aMID = 0.745156, aFAR = 0.70373;
        const double zeffDR7 = 0.312782;  // redshift at which a_scl is evaluated.

        // Get linear matter power spectrum for NEAR, MID, and FAR
        Math::TableFunction<double, double> P_lin_funcNEAR;
        Math::TableFunction<double, double> P_lin_funcMID;
        Math::TableFunction<double, double> P_lin_funcFAR;
        getMatterPs(zNEAR, &P_lin_funcNEAR);
        getMatterPs(zMID, &P_lin_funcMID);
        getMatterPs(zFAR, &P_lin_funcFAR);

        // Create some arrays to store the various power spectra
        int n = P_lin_funcNEAR.size();
        double kvals[n];
        double P_linNEAR[n];
        double P_linMID[n];
        double P_linFAR[n];
        double lnP_linNEAR[n];
        double lnP_linMID[n];
        double lnP_linFAR[n];
        int itemp = 0;
        for(auto const &point : P_lin_funcNEAR)
        {
            kvals[itemp] = point.first; // Same for all
            P_linNEAR[itemp] = point.second;
            lnP_linNEAR[itemp] = std::log(point.second);
            ++itemp;
        }
        itemp = 0;
        for(auto const &point : P_lin_funcMID)
        {
            P_linMID[itemp] = point.second;
            lnP_linMID[itemp] = std::log(point.second);
            ++itemp;
        }
        itemp = 0;
        for(auto const &point : P_lin_funcFAR)
        {
            P_linFAR[itemp] = point.second;
            lnP_linFAR[itemp] = std::log(point.second);
            ++itemp;
        }
       
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

        // Calculate P_damp for NEAR, MID, and FAR using eq. 10 from BR09
        std::vector<double> P_dampNEAR(n);
        std::vector<double> P_dampMID(n);
        std::vector<double> P_dampFAR(n);
        double sigma2BAONEAR = 86.9988, sigma2BAOMID = 85.1374, sigma2BAOFAR = 84.5958;
        double powerscaletoz0NEAR, powerscaletoz0MID, powerscaletoz0FAR;
        double getabstransferscalez0 = 162.71; // TODO: Figure out what these are. Hard coded for fiducial model for now.
        double getabstransferscaleNEAR = 144.31;
        double getabstransferscaleMID = 136.64;
        double getabstransferscaleFAR = 131.3;
        powerscaletoz0NEAR = pow(getabstransferscalez0, 2.0)/pow(getabstransferscaleNEAR, 2.0);
        powerscaletoz0MID = pow(getabstransferscalez0, 2.0)/pow(getabstransferscaleMID, 2.0);
        powerscaletoz0FAR = pow(getabstransferscalez0, 2.0)/pow(getabstransferscaleFAR, 2.0);
        for(int i = 0; i < n; ++i)
        {
            double expval = std::exp(-1.0 * pow(kvals[i], 2) * sigma2BAONEAR * 0.5);
            P_dampNEAR[i] = (P_linNEAR[i] * expval + P_nwNEAR[i] * (1.0 - expval)) * powerscaletoz0NEAR;
            expval = std::exp(-1.0 * pow(kvals[i], 2) * sigma2BAOMID * 0.5);
            P_dampMID[i] = (P_linMID[i] * expval + P_nwMID[i] * (1.0 - expval)) * powerscaletoz0NEAR;
            expval = std::exp(-1.0 * pow(kvals[i], 2) * sigma2BAOFAR * 0.5);
            P_dampFAR[i] = (P_linFAR[i] * expval + P_nwFAR[i] * (1.0 - expval)) * powerscaletoz0NEAR;
        }
        //output_screen("Checkpoint 1" << std::endl);
        // Initalize a cubic spline for each P_damp to be used for calculating P_halo later
        //std::vector<double> kvec(kvals, kvals + sizeof(kvals) / sizeof(double));
        //std::vector<Math::CubicSpline*> P_damp_spline(3);
        //for(int i = 0; i < 3; ++i)
        //{
        //    Math::CubicSpline tempSpline(kvec, P_damp[i]);
        //    P_damp_spline[i] = &tempSpline;
        //}
        std::vector<double> kvec(kvals, kvals + sizeof(kvals) / sizeof(double));
        Math::CubicSpline P_damp_splineNEAR(kvec, P_dampNEAR);
        Math::CubicSpline P_damp_splineMID(kvec, P_dampMID);
        Math::CubicSpline P_damp_splineFAR(kvec, P_dampFAR);
        
        //output_screen("Checkpoint 2" << std::endl);

        // Apply halofit model for nonlinear structure growth to P_nw to generate P_halofitnw
        // Tau is the conformal time at z
        // k_nl is the value of k where power spectrum becomes nonlinear
        double P_halofitnwNEAR[n];
        double P_halofitnwMID[n];
        double P_halofitnwFAR[n];
        double tau;
        double k_nl;
        nonlinear_k_nl_at_z(br_, nl_, zNEAR, &k_nl);
        background_tau_of_z(br_, zNEAR, &tau);
        nonlinear_halofit(pr_, br_, pm_, nl_, tau, P_nwNEAR, P_halofitnwNEAR, &k_nl);
        nonlinear_k_nl_at_z(br_, nl_, zMID, &k_nl);
        background_tau_of_z(br_, zMID, &tau);
        nonlinear_halofit(pr_, br_, pm_, nl_, tau, P_nwMID, P_halofitnwMID, &k_nl);
        nonlinear_k_nl_at_z(br_, nl_, zFAR, &k_nl);
        background_tau_of_z(br_, zFAR, &tau);
        nonlinear_halofit(pr_, br_, pm_, nl_, tau, P_nwFAR, P_halofitnwFAR, &k_nl);

        // Calculate factor r_halofit, the ratio of P_halofitnw to P_nw
        // Then P_DMhalofit = P_damp * r_halofit
        std::vector<double> r_halofitNEAR(n);
        std::vector<double> r_halofitMID(n);
        std::vector<double> r_halofitFAR(n);
        //std::vector<double> P_DMhalofit(n);
        //double temp;
        for(int i = 0; i < n; ++i)
        {
            //temp = P_halofitnw[i]/P_nw[i];
            r_halofitNEAR[i] = P_halofitnwNEAR[i]/P_nwNEAR[i];
            r_halofitMID[i] = P_halofitnwMID[i]/P_nwMID[i];
            r_halofitFAR[i] = P_halofitnwFAR[i]/P_nwFAR[i];
            //P_DMhalofit[i] = P_damp[i] * temp;
        }
        // Create a cubic spline for r_halofit for use later
        Math::CubicSpline r_halofitNEAR_spline(kvec, r_halofitNEAR);
        Math::CubicSpline r_halofitMID_spline(kvec, r_halofitMID);
        Math::CubicSpline r_halofitFAR_spline(kvec, r_halofitFAR);
        //output_screen("Checkpoint 3" << std::endl);

        // Output k, P_lin, P_nw, and r_halofit for NEAR to a file for debugging
        double h = params_->getH();
        std::ofstream out("fidmodelNEAR.txt");
        for(int i = 0; i < n; ++i)
            out << kvals[i]/h << " " << P_linNEAR[i]*pow(h, 3.0) << " " << P_nwNEAR[i]*pow(h, 3.0) << " " << r_halofitNEAR[i] << std::endl;
        out.close();

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
        std::vector<double> k_fid; // k values should be identical so only need one. Is this k or k/h?
        std::vector<double> r_fidNEAR; // vector to hold ratio of nw to nl
        std::vector<double> r_fidMID; // vector to hold ratio of nw to nl
        std::vector<double> r_fidFAR; // vector to hold ratio of nw to nl
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
            r_fidNEAR.push_back(ratiodummy);
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
            r_fidMID.push_back(ratiodummy);
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
            r_fidFAR.push_back(ratiodummy);
        }
        //r_fid.push_back(temp_vec);
        datafile.close();

        // Hard coded for testing purposes. TODO: delete later
        // Read in data from model files for r_nwhalofit
        //std::vector< std::vector<double> > r_nwhalofit(3);
        //datafile.open(root + "models/lrgdr7model_matterpowerzNEAR.dat");
        //// Skip first line
        ////std::string line;
        //std::getline(datafile, line);
        //// Read in rest of data file
        ////temp_vec.clear()
        //while(!datafile.eof())
        //{
        //    std::getline(datafile, line);
        //    std::istringstream iss(line);
        //    double kdummy, plindummy, psmoothdummy, ratiodummy;
        //    iss >> kdummy >> plindummy >> psmoothdummy >> ratiodummy; 
        //    //k_fid.push_back(kdummy);
        //    //Plin_fid[i] = plindummy;
        //    //Psmooth_fid[i] = psmoothdummy;
        //    r_nwhalofit[0].push_back(ratiodummy);
        //}
        ////r_fid.push_back(temp_vec);
        //datafile.close();
        //datafile.open(root + "models/lrgdr7model_matterpowerzMID.dat");
        //// Skip first line
        ////std::string line;
        //std::getline(datafile, line);
        //// Read in rest of data file
        ////temp_vec.clear()
        //while(!datafile.eof())
        //{
        //    std::getline(datafile, line);
        //    std::istringstream iss(line);
        //    double kdummy, plindummy, psmoothdummy, ratiodummy;
        //    iss >> kdummy >> plindummy >> psmoothdummy >> ratiodummy; 
        //    //k_fid.push_back(kdummy);
        //    //Plin_fid[i] = plindummy;
        //    //Psmooth_fid[i] = psmoothdummy;
        //    r_nwhalofit[1].push_back(ratiodummy);
        //}
        ////r_fid.push_back(temp_vec);
        //datafile.close();
        //datafile.open(root + "models/lrgdr7model_matterpowerzFAR.dat");
        //// Skip first line
        ////std::string line;
        //std::getline(datafile, line);
        //// Read in rest of data file
        ////temp_vec.clear()
        //while(!datafile.eof())
        //{
        //    std::getline(datafile, line);
        //    std::istringstream iss(line);
        //    double kdummy, plindummy, psmoothdummy, ratiodummy;
        //    iss >> kdummy >> plindummy >> psmoothdummy >> ratiodummy; 
        //    //k_fid.push_back(kdummy);
        //    //Plin_fid[i] = plindummy;
        //    //Psmooth_fid[i] = psmoothdummy;
        //    r_nwhalofit[2].push_back(ratiodummy);
        //}
        ////r_fid.push_back(temp_vec);
        //datafile.close();

        //output_screen("Checkpoint 4" << std::endl);
        
        // Set number of k values in fiducial model
        int k_size = k_fid.size();

        // Set weights to do weighted sum for P_halo as in eq. 17 of BR09
        double zweightNEAR = 0.395;
        double zweightMID = 0.355;
        double zweightFAR = 0.250;

        // Now calculate P_halo as P_damp * r_fid * fid_polys
        std::vector<double> P_halo(k_size, 0.0);
        ////std::ofstream out("k_Pdamp_rfid_fidpoly.txt");
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
        //std::ofstream out("r_halofit.txt");
        //for(int i = 0; i < k_size; ++i)
        //{
        //    out << k_fid[i] << " " << r_halofit_spline.evaluate(k_fid[i]) << std::endl;
        //}
        //out.close();
        for(int i = 0; i < k_size; ++i)
        {
            double k = k_fid[i]*h;
            std::vector<double> fidpolys;
            LRGtoICsmooth(k_fid[i], fidpolys);
            //LRGtoICsmooth(k, fidpolys);
            double nlratNEAR = r_halofitNEAR_spline.evaluate(k) / r_fidNEAR[i];
            double nlratMID = r_halofitMID_spline.evaluate(k) / r_fidMID[i];
            double nlratFAR = r_halofitFAR_spline.evaluate(k) / r_fidFAR[i];
            //double nlratNEAR = r_nwhalofit[0][i] / r_fid[0][i];
            //double nlratMID = r_nwhalofit[1][i] / r_fid[1][i];
            //double nlratFAR = r_nwhalofit[2][i] / r_fid[2][i];
            //P_halo[i] = zweight[0] * P_damp_splineNEAR.evaluate(k_fid[i]) * r_fid[0][i] * fidpolys[0]
            //            + zweight[1] * P_damp_splineMID.evaluate(k_fid[i]) * r_fid[1][i] * fidpolys[1]
            //            + zweight[2] * P_damp_splineFAR.evaluate(k_fid[i]) * r_fid[2][i] * fidpolys[2];
            P_halo[i] = zweightNEAR * P_damp_splineNEAR.evaluate(k) * nlratNEAR * fidpolys[0]
                        + zweightMID * P_damp_splineMID.evaluate(k) * nlratMID * fidpolys[1]
                        + zweightFAR * P_damp_splineFAR.evaluate(k) * nlratFAR * fidpolys[2];
            //P_halo[i] = P_damp_splineNEAR.evaluate(k_fid[i]) * r_fid[0][i] * fidpolys[0];
        }

        //output_screen("Checkpoint 5" << std::endl);

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
        for(int i = 0; i < k_size; ++i)
            (*ps)[k_fid[i]] = P_halo[i];
    }

    double getR_NL()
    {
        double z, k_nl, k;
        double* pvecback;
        int junk;
        double r_nl;
        int index_tau; // Think I can get rid of this and replace last instance with nl_->tau[index_tau]


        pt_->k_max_for_pk = kMax_;
        nl_->method = nl_halofit;
        pr_->halofit_dz = 0.1;
        pr_->halofit_min_k_nonlinear = 0.0035;
        pr_->halofit_sigma_precision = 0.05;

        nonlinear_free(nl_);
        nonlinear_init(pr_, br_, th_, pt_, pm_, nl_);

        //nl_->method = nl_halofit;
        z = 0.0;
        nonlinear_k_nl_at_z(br_, nl_, z, &k_nl);


        pvecback = (double*) calloc(br_->bg_size,sizeof(double));

        for(index_tau = 0; index_tau < nl_->tau_size; ++index_tau)
        {
            background_at_tau(br_, nl_->tau[index_tau], br_->short_info,
                                br_->inter_normal, &junk, pvecback);
        }


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
