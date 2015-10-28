#include <fstream>
#include <string>
#include <sstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <likelihood_function.hpp>

#include <table_function.hpp>

#include <cosmological_params.hpp>
#include <mcmc.hpp>
#include <mn_scanner.hpp>
#include <cmb.hpp>
#include <planck_like.hpp>

//#include "complete_likelihood.hpp"

const double kMin = 1e-6;
const double kMax = 1.0;
const double aMin = -2;
const double aMax = 4;

const double nT = 0.0;
const double pivotT = 0.05;

class RLike : public Math::RealFunction
{
public:
    RLike(const char* fileName)
    {
        std::ifstream in(fileName);
        StandardException exc;
        if(!in)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Cannot open input file " << fileName << ".";
            exc.set(exceptionStr.str());
            throw exc;
        }

        min_ = std::numeric_limits<double>::max();
        max_ = std::numeric_limits<double>::min();

        while(!in.eof())
        {
            std::string s;
            std::getline(in, s);
            if(s == "")
                break;

            double x, y;
            std::stringstream str(s);
            str >> x >> y;
            if(x < 0 || y < 0)
            {
                std::stringstream exceptionStr;
                exceptionStr << "Invalid input line: " << s;
                exc.set(exceptionStr.str());
                throw exc;
            }
            rLike_[x] = y;
            if(x < min_)
                min_ = x;
            if(x > max_)
                max_ = x;
        }

        output_screen("r likelihood initialized successfully between " << min_ << " and " << max_ << std::endl);
    }

    double min() const { return min_; }
    double max() const { return max_; }

    virtual double evaluate(double x) const
    {
        if(x < min_ || x > max_)
            return 1e-30;

        return rLike_.evaluate(x);
    }
private:
    Math::TableFunction<double, double> rLike_;
    double min_, max_;
};

class LinSplWithDegenerateNeutrinosParams : public LinearSplineParams
{
public:
    LinSplWithDegenerateNeutrinosParams(double omBH2, double omCH2, double h, double tau, double r, double nt, double pivot, const std::vector<double>& kVals, const std::vector<double>& amplitudes,  double nEff, int nMassive, double sumMNu) : LinearSplineParams(omBH2, omCH2, h, tau, kVals, amplitudes), r_(r), nEff_(nEff), nMassive_(nMassive), sumMNu_(sumMNu)
    {
        check(nEff > 0, "invalid nEff = " << nEff);
        check(sumMNu >= 0, "invalid sumMNu = " << sumMNu);
        check(nMassive >= 0, "number of massive neutrinos is negative: " << nMassive);
        check(nEff > nMassive, "nEff needs to be more than the number of massive neutrinos");

        //psTensor_ = new StandardPowerSpectrumTensor(powerSpectrum(), r, nt, pivot);
        psTensor_ = new StandardPowerSpectrumTensor(r * 2.1867e-9, nt, pivot);
        check(r_ >= 0, "invalid r");
    }

    ~LinSplWithDegenerateNeutrinosParams() { delete psTensor_; }

    virtual const Math::RealFunction& powerSpectrumTensor() const { return *psTensor_; }
    virtual double getR() const { return r_; }
    virtual double getNt() const { return psTensor_->getNt(); }

    virtual double getNEff() const { return nEff_ - nMassive_; }
    virtual int getNumNCDM() const { return nMassive_; }
    virtual double getNCDMParticleMass(int i) const
    {
        check(i >= 0 && i < nMassive_, "invalid index = " << i);
        return sumMNu_ / nMassive_;
    }

    virtual double getNCDMParticleTemp(int i) const
    {
        check(i >= 0 && i < nMassive_, "invalid index = " << i);
        //return 0.715985;
        return 0.713765855506013;
    }

private:
    double nEff_;
    int nMassive_;
    double sumMNu_;
    double r_;
    StandardPowerSpectrumTensor* psTensor_;
};

class CubSplWithDegenerateNeutrinosParams : public CubicSplineParams
{
public:
    CubSplWithDegenerateNeutrinosParams(double omBH2, double omCH2, double h, double tau, double r, double nt, double pivot, const std::vector<double>& kVals, const std::vector<double>& amplitudes,  double nEff, int nMassive, double sumMNu) : CubicSplineParams(omBH2, omCH2, h, tau, kVals, amplitudes), r_(r), nEff_(nEff), nMassive_(nMassive), sumMNu_(sumMNu)
    {
        check(nEff > 0, "invalid nEff = " << nEff);
        check(sumMNu >= 0, "invalid sumMNu = " << sumMNu);
        check(nMassive >= 0, "number of massive neutrinos is negative: " << nMassive);
        check(nEff > nMassive, "nEff needs to be more than the number of massive neutrinos");

        //psTensor_ = new StandardPowerSpectrumTensor(powerSpectrum(), r, nt, pivot);
        psTensor_ = new StandardPowerSpectrumTensor(r * 2.1867e-9, nt, pivot);
        check(r_ >= 0, "invalid r");
    }

    ~CubSplWithDegenerateNeutrinosParams() { delete psTensor_; }

    virtual const Math::RealFunction& powerSpectrumTensor() const { return *psTensor_; }
    virtual double getR() const { return r_; }
    virtual double getNt() const { return psTensor_->getNt(); }

    virtual double getNEff() const { return nEff_ - nMassive_; }
    virtual int getNumNCDM() const { return nMassive_; }
    virtual double getNCDMParticleMass(int i) const
    {
        check(i >= 0 && i < nMassive_, "invalid index = " << i);
        return sumMNu_ / nMassive_;
    }

    virtual double getNCDMParticleTemp(int i) const
    {
        check(i >= 0 && i < nMassive_, "invalid index = " << i);
        //return 0.715985;
        return 0.713765855506013;
    }

private:
    double nEff_;
    int nMassive_;
    double sumMNu_;
    double r_;
    StandardPowerSpectrumTensor* psTensor_;
};

// parameters come in the following order: 4 base params (ombh2, omch2, h, tau), nEff (if varied), sumMNu (if varied), r (if varied), background params (14 for planck, 2 for own like), ps params (2 for standard case, however many for spline)
class CMBLike : public Math::LikelihoodFunction
{
public:
    enum PSType{ STANDARD_PS = 0, LINEAR_SPLINE, CUBIC_SPLINE, CUTOFF, RUN, PSTYPE_MAX };
public:
    CMBLike(bool useClik, PSType psType = STANDARD_PS, int numberOfNodes = 0, bool varyStandard = true, bool varyBackground = true, bool varyNEff = false, bool varySumMNu = false, bool varyR = false, bool usePlanckPol = false) : lMax_(analysisLMax + 1000), pivot_(0.05), planckLike_(NULL), psType_(psType), numberOfNodes_(-1), varyStandard_(varyStandard), varyBackground_(varyBackground), varyNEff_(varyNEff), varySumMNu_(varySumMNu), varyR_(varyR)
    {
        check(psType_ >= STANDARD_PS && psType_ < PSTYPE_MAX, "invalid psType");

        nPar_ = 0;

        if(varyStandard_)
            nPar_ += 4;

        if(varyNEff)
            ++nPar_;

        if(varySumMNu)
            ++nPar_;

        if(varyR)
            ++nPar_;

        nParBase_ = nPar_;

        if(useClik)
        {
            planckLike_ = new PlanckLikelihood(true, true, true, usePlanckPol, false, varyR);
            if(varyBackground_)
                nPar_ += 14;
        }
        else
        {
            cmb_.preInitialize(lMax_, false, true, varyR, lMax_);
            if(varyBackground_)
                nPar_ += 2;
        }

        nParNoPs_ = nPar_;

        if(psType_ == STANDARD_PS)
            nPar_ += 2;
        else if(psType_ == CUTOFF || psType_ == RUN)
            nPar_ += 3;
        else
        {
            check(numberOfNodes >= 0, "");
            nPar_ += (2 + 2 * numberOfNodes);
            numberOfNodes_ = numberOfNodes;
        }

        /*
        if(varyR)
            rLike_ = new RLike("r_like.txt");
        */
    }

    ~CMBLike()
    {
        if(planckLike_)
            delete planckLike_;

        /*
        if(varyR_)
            delete rLike_;
        */
    }

    int getNPar() const { return nPar_; }

    virtual double calculate(double* p, int nParams)
    {
        check(nParams == nPar_, "");

        bool bigLike = false; // return a big number if the nodes are out of order, also if the primordial power spectrum is crazy

        CosmologicalParams* params;
        double as;
        std::vector<double> kVals, amplitudes;
        double kPrev = kMin;

        int nStandard = (varyStandard_ ? 4 : 0);

        // fixed values
        double omBH2 = (planckLike_ ? 0.0220877 : 0.0224029);
        double omCH2 = (planckLike_ ? 0.118108 : 0.115274);
        double h = (planckLike_ ? 0.684021 : 0.698577);
        double tau = (planckLike_ ? 0.0932567 : 0.103401);

        if(varyStandard_)
        {
            omBH2 = p[0];
            omCH2 = p[1];
            h = p[2];
            tau = p[3];
        }

        double nEff = 3.046;
        if(varyNEff_)
        {
            nEff = p[nStandard];
            ++nStandard;
        }

        int nMassive = 0;
        double sumMNu = 0;
        if(varySumMNu_)
        {
            nMassive = 1;
            sumMNu = p[nStandard];
            ++nStandard;
        }

        double r = 0;
        if(varyR_)
        {
            r = p[nStandard];
            ++nStandard;
        }

        switch(psType_)
        {
            case STANDARD_PS:
                as = std::exp(p[nParNoPs_ + 1]) / 1e10;
                params = new LCDMWithTensorAndDegenerateNeutrinosParams(omBH2, omCH2, h, tau, p[nParNoPs_], as, pivot_, r, nT, pivotT, nEff, nMassive, sumMNu);
                break;

            case LINEAR_SPLINE:
                kVals.push_back(kMin);

                for(int i = 0; i < numberOfNodes_; ++i)
                {
                    const double kValue = std::exp(p[nParNoPs_ + i]);
                    kVals.push_back(kValue);
                    if(kValue < kPrev)
                        bigLike = true;

                    kPrev = kValue;
                }
                kVals.push_back(kMax);

                for(int i = 0; i < numberOfNodes_ + 2; ++i)
                    amplitudes.push_back(std::exp(p[nParNoPs_ + numberOfNodes_ + i]) / 1e10);

                params = new LinSplWithDegenerateNeutrinosParams(omBH2, omCH2, h, tau, r, nT, pivotT, kVals, amplitudes, nEff, nMassive, sumMNu);
                break;
                
            case CUBIC_SPLINE:
                kVals.push_back(kMin);
                for(int i = 0; i < numberOfNodes_; ++i)
                {
                    const double kValue = std::exp(p[nParNoPs_ + i]);
                    kVals.push_back(kValue);
                    if(kValue < kPrev)
                        bigLike = true;

                    kPrev = kValue;
                }
                kVals.push_back(kMax);

                for(int i = 0; i < numberOfNodes_ + 2; ++i)
                    amplitudes.push_back(std::exp(p[nParNoPs_ + numberOfNodes_ + i]) / 1e10);

                params = new CubSplWithDegenerateNeutrinosParams(omBH2, omCH2, h, tau, r, nT, pivotT, kVals, amplitudes, nEff, nMassive, sumMNu);
                break;

            case CUTOFF:
                as = std::exp(p[nParNoPs_ + 1]) / 1e10;
                params = new LCDMWithCutoffTensorDegenerateNeutrinosParams(omBH2, omCH2, h, tau, std::exp(p[nParNoPs_ + 2]), p[nParNoPs_], as, pivot_, r, nT, pivotT, nEff, nMassive, sumMNu);
                break;

            case RUN:
                as = std::exp(p[nParNoPs_ + 1]) / 1e10;
                params = new LCDMWithTensorAndDegenerateNeutrinosParams(omBH2, omCH2, h, tau, p[nParNoPs_], as, pivot_, r, nT, pivotT, nEff, nMassive, sumMNu, p[nParNoPs_ + 2]);
                break;

            default:
                check(false, "");
        }

        if(crazyPS(params))
            bigLike = true;

        if(bigLike)
        {
            delete params;
            output_screen("Likelihood = " << 1e10 << std::endl);
            return 1e10;
        }

        double like;

        if(planckLike_)
        {
            check(nParNoPs_ - nParBase_ == (varyBackground_ ? 14 : 0), "");

            planckLike_->setCosmoParams(*params);
            if(varyBackground_)
                planckLike_->setCamspecExtraParams(p[nParBase_], p[nParBase_ + 1], p[nParBase_ + 2], p[nParBase_ + 3], p[nParBase_ + 4], p[nParBase_ + 5], p[nParBase_ + 6], p[nParBase_ + 7], p[nParBase_ + 8], p[nParBase_ + 9], p[nParBase_ + 10], p[nParBase_ + 11], p[nParBase_ + 12], p[nParBase_ + 13]);
            else
                planckLike_->setCamspecExtraParams(168, 66, 98, 3.5, 34.5, 4.9, 0.971, 0.37, 0.65, 1.0005, 0.996, 0.47, 1.5, 0.6);
            planckLike_->calculateCls();

            like = planckLike_->likelihood();
        }

        //else
        //{
        //    check(nParNoPs_ - nParBase_ == (varyBackground_ ? 2 : 0), "");

        //    cmb_.initialize(*params, true, false, true);
        //    std::vector<double> cl;
        //    cmb_.getLensedCl(&cl);
        //    check(cl.size() >= analysisLMax + 1, "");

        //    double A_ps = 6.2, A_cl = 51;

        //    if(varyBackground_)
        //    {
        //        A_ps = p[nParBase_] * 2.0 * Math::pi / (3000 * 3001);
        //        A_cl = p[nParBase_ + 1] / std::pow(3000.0, 0.8);
        //    }

        //    std::vector<double> clVec(analysisLMax + 1, 0);

        //    for(int l = 2; l <= analysisLMax; ++l)
        //        clVec[l] = cl[l] + A_ps + A_cl * std::pow(double(l), 0.8) * 2.0 * Math::pi / double(l * (l + 1));

        //    double lowChi2, lowLogDet, highL;

        //    like = CompleteLikelihood::create().calculate(clVec, lowChi2, lowLogDet, highL);
        //}

        delete params;

        /*
        if(varyR_)
            like -= 2 * std::log(rLike_->evaluate(r));
        */

        output_screen("Likelihood = " << like << std::endl);
        return like;
    }

private:
    bool crazyPS(const CosmologicalParams* params) const
    {
        if(psType_ == CUTOFF)
            return false;

        const double psMin = 1e-12, psMax = 1e-7;
        const int n = 1000;
        const double delta = (std::log(kMax) - std::log(kMin)) / n;
        for(int i = 0; i <= n; ++i)
        {
            const double k = std::exp(std::log(kMin) + delta * i);
            const double a = params->powerSpectrum().evaluate(k);
            if(a < psMin || a > psMax)
                return true;
        }
        return false;
    }

private:
    CMB cmb_;
    PlanckLikelihood* planckLike_;
    int lMax_;
    double pivot_;
    int nPar_;
    int nParNoPs_;
    int nParBase_;
    PSType psType_;
    int numberOfNodes_;
    bool varyStandard_;
    bool varyBackground_;
    bool varyNEff_;
    bool varySumMNu_;
    bool varyR_;

    //RLike* rLike_;
};

int main(int argc, char *argv[])
{
    try {
        StandardException exc;
        if(argc < 3)
        {
            std::string exceptionStr = "Likelihood type (planck or standard), power spectrum type (standard, linear_spline, cubic_spline, cutoff, or running), and number of nodes (optional, not used for standard ps or cutoff, used 0 for others if not specified) must be specified. If want to vary n_eff, give true after the number of nodes(optional, default is false). If want to vary sum(m_nu), give true after the last one (optional, default is false). If do not want to vary standard params, give false after the last one (optional, default is true). If do not want to vary background params, give false after the last one (optional, default is true). If want to vary r, give true after the last one (optional, default = false). If want to use polarization data with planck (only with planck, not standard), give true after the last one (optional, default = false).";
            exc.set(exceptionStr);
            throw exc;
        }

        bool useClik = true;
        if(argv[1][0] == 's')
            useClik = false;

        CMBLike::PSType psType = CMBLike::STANDARD_PS;
        if(argv[2][0] == 'l')
            psType = CMBLike::LINEAR_SPLINE;
        else if(argv[2][0] == 'r')
            psType = CMBLike::RUN;
        else if(argv[2][0] == 'c' && argv[2][2] == 'b')
            psType = CMBLike::CUBIC_SPLINE;
        else if(argv[2][0] == 'c' && argv[2][2] == 't')
            psType = CMBLike::CUTOFF;

        int numberOfNodes = 0;
        if(argc > 3)
        {
            std::stringstream str;
            str << argv[3];
            str >> numberOfNodes;

            if(numberOfNodes < 0)
            {
                std::stringstream exceptionStr;
                exceptionStr << "Invalid number of nodes " << numberOfNodes << ", must be non-negative.";
                exc.set(exceptionStr.str());
                throw exc;
            }
        }

        bool varyNEff = false, varySumMNu = false, varyStandard = true, varyBackground = true, varyR = false, usePlanckPol = false;

        if(argc > 4 && argv[4][0] == 't')
            varyNEff = true;

        if(argc > 5 && argv[5][0] == 't')
            varySumMNu = true;

        if(argc > 6 && argv[6][0] == 'f')
            varyStandard = false;

        if(argc > 7 && argv[7][0] == 'f')
            varyBackground = false;

        if(argc > 8 && argv[8][0] == 't')
            varyR = true;

        if(argc > 9 && argv[9][0] == 't')
            usePlanckPol = true;


        std::stringstream fileRoot;
        fileRoot << "mn_";
        if(useClik)
            fileRoot << "planck_";
        else
            fileRoot << "stdlike_";

        switch (psType)
        {
        case CMBLike::STANDARD_PS:
            fileRoot << "stdps_";
            break;

        case CMBLike::LINEAR_SPLINE:
            fileRoot << "linsplps_" << numberOfNodes << "_";
            break;

        case CMBLike::CUBIC_SPLINE:
            fileRoot << "cubsplps_" << numberOfNodes << "_";
            break;

        case CMBLike::CUTOFF:
            fileRoot << "cutoff_";
            break;

        case CMBLike::RUN:
            fileRoot << "run_";
            break;

        default:
            check(false, "");
            break;
        }

        if(!varyStandard)
            fileRoot << "fixedst_";

        if(!varyBackground)
            fileRoot << "fixedback_";

        int nParBase = (varyStandard ? 4 : 0);

        if(varyNEff)
        {
            ++nParBase;
            fileRoot << "neff_";
        }

        if(varySumMNu)
        {
            ++nParBase;
            fileRoot << "summnu_";
        }

        if(varyR)
        {
            ++nParBase;
            fileRoot << "r_";
        }

        if(usePlanckPol)
            fileRoot << "pol_";


        CMBLike like(useClik, psType, numberOfNodes, varyStandard, varyBackground, varyNEff, varySumMNu, varyR, usePlanckPol);
        MnScanner scanner(like.getNPar(), like, 300, fileRoot.str());
        //Math::MetropolisHastings scanner(like.getNPar(), like, fileRoot.str());

        int nStandard = 0;
        if(varyStandard)
        {
            scanner.setParam(0, std::string("ombh2"), 0.02, 0.025);
            scanner.setParam(1, std::string("omch2"), 0.1, 0.2);
            scanner.setParam(2, std::string("h"), 0.55, 0.85);

            if(useClik)
                scanner.setParam(3, std::string("tau"), 0.02, 0.20);
            else
                scanner.setParamGauss(3, std::string("tau"), 0.0851, 0.014);

            nStandard = 4;
        }

        if(varyNEff)
        {
            scanner.setParam(nStandard, std::string("n_eff"), 2.0, 5.0);
            ++nStandard;
        }

        if(varySumMNu)
        {
            scanner.setParam(nStandard, std::string("sum_mnu"), 0.01, 3.0);
            ++nStandard;
        }
        
        if(varyR)
        {
            RLike rLike("r_like.txt");
            scanner.setParamGeneral(nStandard, std::string("r"), rLike.min(), rLike.max(), rLike);
        }

        int nParNoPs = nParBase;
        if(varyBackground)
        {
            if(useClik)
            {
                scanner.setParam(nParBase, "A_ps_100", 0, 360);
                scanner.setParam(nParBase + 1, "A_ps_143", 0, 270);
                scanner.setParam(nParBase + 2, "A_ps_217", 0, 450);
                scanner.setParam(nParBase + 3, "A_cib_143", 0, 20);
                scanner.setParam(nParBase + 4, "A_cib_217", 0, 80);
                scanner.setParam(nParBase + 5, "A_sz", 0, 10);
                scanner.setParam(nParBase + 6, "r_ps", 0.0, 1.0);
                scanner.setParam(nParBase + 7, "r_cib", 0.0, 1.0);
                scanner.setParam(nParBase + 8, "n_Dl_cib", -2, 2);
                scanner.setParam(nParBase + 9, "cal_100", 0.98, 1.02);
                scanner.setParam(nParBase + 10, "cal_127", 0.95, 1.05);
                scanner.setParam(nParBase + 11, "xi_sz_cib", 0, 1);
                scanner.setParam(nParBase + 12, "A_ksz", 0, 10);
                scanner.setParam(nParBase + 13, "Bm_1_1", -20, 20);
                nParNoPs = nParBase + 14;
            }
            else
            {
                scanner.setParam(nParBase, std::string("A_ps"), 0, 200);
                scanner.setParam(nParBase + 1, std::string("A_cl"), 0, 60);
                nParNoPs = nParBase + 2;
            }
        }

        if(psType == CMBLike::STANDARD_PS)
        {
            scanner.setParam(nParNoPs, "ns", 0.9, 1.1);
            scanner.setParam(nParNoPs + 1, "as", 2.7, 3.5);
        }
        else if(psType == CMBLike::CUTOFF)
        {
            scanner.setParam(nParNoPs, "ns", 0.9, 1.1);
            scanner.setParam(nParNoPs + 1, "as", 2.7, 3.5);
            scanner.setParam(nParNoPs + 2, "cut", std::log(kMin), std::log(1e-2));
        }
        else if(psType == CMBLike::RUN)
        {
            scanner.setParam(nParNoPs, "ns", 0.9, 1.1);
            scanner.setParam(nParNoPs + 1, "as", 2.7, 3.5);
            scanner.setParam(nParNoPs + 2, "run", -0.1, 0.1);
        }
        else
        {
            for(int i = 0; i < numberOfNodes; ++i)
            {
                std::stringstream paramName;
                paramName << "k_node_" << i;
                scanner.setParam(nParNoPs + i, paramName.str().c_str(), std::log(kMin), std::log(kMax));
            }
            scanner.setParam(nParNoPs + numberOfNodes, "A_lower", aMin, aMax);
            for(int i = 0; i < numberOfNodes; ++i)
            {
                std::stringstream paramName;
                paramName << "A_node_" << i;
                scanner.setParam(nParNoPs + numberOfNodes + 1 + i, paramName.str().c_str(), aMin, aMax);
            }
            scanner.setParam(nParNoPs + 2 * numberOfNodes + 1, "A_upper", aMin, aMax);
        }

        scanner.run();
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}

