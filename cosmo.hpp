#include <fstream>
#include <cmb.hpp>
#include <class.h>

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
};
