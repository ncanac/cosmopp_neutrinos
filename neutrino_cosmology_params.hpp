#include <cosmological_params.hpp>

class DegenerateNeutrinosParams : public LambdaCDMParams
{
public:
    DegenerateNeutrinosParams(double omBH2, double omCH2, double h, double tau, double ns, double as, double pivot, double nEff, int nMassive, double sumMNu) : LambdaCDMParams(omBH2, omCH2, h, tau, ns, as, pivot), nEff_(nEff), nMassive_(nMassive), sumMNu_(sumMNu)
    {
        check(nEff > 0, "invalid nEff = " << nEff);
        check(sumMNu >= 0, "invalid sumMNu = " << sumMNu);
        check(nMassive >= 0, "number of massive neutrinos is negative: " << nMassive);
        check(nEff > nMassive, "nEff needs to be more than the number of massive neutrinos");
    }

    ~DegenerateNeutrinosParams() {}

    virtual double getNEff() const { return nEff_; }
    virtual int getNumNCDM() const { return nMassive_; }
    double getSumMNu() const { return sumMNu_; }
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

    virtual std::string name() const { return "DegenerateNeutrinos"; }

    virtual void getAllParameters(std::vector<double>& v) const
    { 
        v.resize(8);
        v[0] = getOmBH2();
        v[1] = getOmCH2();
        v[2] = getH();
        v[3] = getTau();
        v[4] = getNs();
        v[5] = std::log(getAs() * 1e10);
        v[6] = getNEff();
        v[7] = getSumMNu();
    }

    virtual bool setAllParameters(const std::vector<double>& v, double *badLike = NULL)
    {
        check(v.size() >= 8, "");
        omBH2_ = v[0];
        omCH2_ = v[1];
        h_ = v[2];
        tau_ = v[3];
        ps_.setNs(v[4]);
        ps_.setAs(std::exp(v[5]) / 1e10);
        nEff_ = v[6];
        sumMNu_ = v[7]; 

        if(badLike)
            *badLike = 0;

        return true;
    }

private:
    double nEff_;
    int nMassive_;
    double sumMNu_;
};

class LinearSplineParams2 : public CosmologicalParams
{
    class DummyPS : public Math::RealFunction
    {
    public:
        DummyPS() { }
        double evaluate(double x) const { return 0.0; }
    };

public:
    LinearSplineParams2(double omBH2, double omCH2, double h, double tau, const std::vector<double>& kVals, const std::vector<double>& amplitudes) : CosmologicalParams(), omBH2_(omBH2), omCH2_(omCH2), h_(h), tau_(tau), ps_(kVals, amplitudes) {}
    ~LinearSplineParams2() {}

    virtual double getOmBH2() const { return omBH2_; }
    virtual double getOmCH2() const { return omCH2_; }
    virtual double getH() const { return h_; }
    virtual double getOmB() const { return omBH2_ / (h_ * h_); }
    virtual double getOmC() const { return omCH2_ / (h_ * h_); }
    virtual double getOmLambda() const { return 1.0 - CosmologicalParams::getOmM(); }
    virtual double getOmK() const { return 0.0; }
    virtual double getNs() const { return ps_.getNs(); }
    virtual double getAs() const { return ps_.getAs(); }
    virtual double getPivot() const { return ps_.getPivot(); }
    virtual double getTau() const { return tau_; }
    virtual double getNEff() const { return 3.046; }
    virtual int getNumNCDM() const { return 0; }
    virtual double getNCDMParticleMass(int i) const { check(false, ""); }
    virtual double getNCDMParticleTemp(int i) const { check(false, ""); }
    virtual double getYHe() const { return 0.0; } // means use BBN
    virtual double getR() const { return 0.0; }
    virtual double getNt() const { return 0.0; }

    virtual const Math::RealFunction& powerSpectrum() const { return ps_; }
    virtual const Math::RealFunction& powerSpectrumTensor() const { return psTensor_; }

    virtual std::string name() const { return "LinearSpline"; }
    virtual void getAllParameters(std::vector<double>& v) const { check(false, "not implemented"); }
    virtual void setAllParameters(const std::vector<double>& v) { check(false, "not implemented"); }

private:
    double omBH2_;
    double omCH2_;
    double h_;
    double tau_;
    LinearSplinePowerSpectrum ps_;
    DummyPS psTensor_;
};

class NeutrinosAndLinearSplineParams : public LinearSplineParams2
{
public:
    NeutrinosAndLinearSplineParams(double omBH2, double omCH2, double h, double tau, const std::vector<double>& kVals, const std::vector<double>& amplitudes, double nEff, int nMassive, double sumMNu) : LinearSplineParams2(omBH2, omCH2, h, tau, kVals, amplitudes), nEff_(nEff), nMassive_(nMassive), sumMNu_(sumMNu)
    {
        check(nEff > 0, "invalid nEff = " << nEff);
        check(sumMNu >= 0, "invalid sumMNu = " << sumMNu);
        check(nMassive >= 0, "number of massive neutrinos is negative: " << nMassive);
        check(nEff > nMassive, "nEff needs to be more than the number of massive neutrinos");
    }

    ~NeutrinosAndLinearSplineParams() {}

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

    virtual std::string name() const { return "NeutrinosAndLinearSpline"; }

private:
    double nEff_;
    int nMassive_;
    double sumMNu_;
};

class CubicSplineParams2 : public CosmologicalParams
{
    class DummyPS : public Math::RealFunction
    {
    public:
        DummyPS() { }
        double evaluate(double x) const { return 0.0; }
    };

public:
    CubicSplineParams2(double omBH2, double omCH2, double h, double tau, const std::vector<double>& kVals, const std::vector<double>& amplitudes) : CosmologicalParams(), omBH2_(omBH2), omCH2_(omCH2), h_(h), tau_(tau), ps_(kVals, amplitudes) {}
    ~CubicSplineParams2() {}

    virtual double getOmBH2() const { return omBH2_; }
    virtual double getOmCH2() const { return omCH2_; }
    virtual double getH() const { return h_; }
    virtual double getOmB() const { return omBH2_ / (h_ * h_); }
    virtual double getOmC() const { return omCH2_ / (h_ * h_); }
    virtual double getOmLambda() const { return 1.0 - CosmologicalParams::getOmM(); }
    virtual double getOmK() const { return 0.0; }
    virtual double getNs() const { return ps_.getNs(); }
    virtual double getAs() const { return ps_.getAs(); }
    virtual double getPivot() const { return ps_.getPivot(); }
    virtual double getTau() const { return tau_; }
    virtual double getNEff() const { return 3.046; }
    virtual int getNumNCDM() const { return 0; }
    virtual double getNCDMParticleMass(int i) const { check(false, ""); }
    virtual double getNCDMParticleTemp(int i) const { check(false, ""); }
    virtual double getYHe() const { return 0.0; } // means use BBN
    virtual double getR() const { return 0.0; }
    virtual double getNt() const { return 0.0; }

    virtual const Math::RealFunction& powerSpectrum() const { return ps_; }
    virtual const Math::RealFunction& powerSpectrumTensor() const { return psTensor_; }

    virtual std::string name() const { return "CubicSpline"; }
    virtual void getAllParameters(std::vector<double>& v) const { check(false, "not implemented"); }
    virtual void setAllParameters(const std::vector<double>& v) { check(false, "not implemented"); }

private:
    double omBH2_;
    double omCH2_;
    double h_;
    double tau_;
    CubicSplinePowerSpectrum ps_;
    DummyPS psTensor_;
};

class NeutrinosAndCubicSplineParams : public CubicSplineParams2
{
public:
    NeutrinosAndCubicSplineParams(double omBH2, double omCH2, double h, double tau, const std::vector<double>& kVals, const std::vector<double>& amplitudes, double nEff, int nMassive, double sumMNu) : CubicSplineParams2(omBH2, omCH2, h, tau, kVals, amplitudes), nEff_(nEff), nMassive_(nMassive), sumMNu_(sumMNu)
    {
        check(nEff > 0, "invalid nEff = " << nEff);
        check(sumMNu >= 0, "invalid sumMNu = " << sumMNu);
        check(nMassive >= 0, "number of massive neutrinos is negative: " << nMassive);
        check(nEff > nMassive, "nEff needs to be more than the number of massive neutrinos");
    }

    ~NeutrinosAndCubicSplineParams() {}

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

    virtual std::string name() const { return "NeutrinosAndCubicSpline"; }

private:
    double nEff_;
    int nMassive_;
    double sumMNu_;
};
