// Greens function class for 2d Green's Function for 2d annulus with radial and
// axial dependence. Inner boundary is radiative (rad) (reaction event), outer
// boundary is absorbing (abs) (escape event). Different "draw" functions
// provide a way to draw certain values from the Green's Function, e.g. an
// escape angle theta ("drawTheta" function).
//
// Written by Laurens Bossen. Adapted by Martijn Wehrens.
// FOM Institute AMOLF.
//
// Modified by Toru Niina.
// RIKEN BDR.

#ifndef GREENS_FUNCTION_2D_RAD_ABS_HPP
#define GREENS_FUNCTION_2D_RAD_ABS_HPP

#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/function.hpp>
#include <boost/array.hpp>
#include <boost/math/constants/constants.hpp>
#include <gsl/gsl_roots.h>
#include "PairGreensFunction.hpp"

namespace greens_functions
{

class GreensFunction2DRadAbs
    : public PairGreensFunction
{
  public:
    typedef std::vector<Real> RealVector;

  public:

    GreensFunction2DRadAbs(const Real D,     const Real kf, const Real r0,
                           const Real Sigma, const Real a)
        : PairGreensFunction(D, kf, r0, Sigma), a(a),
          h(kf / (D * 2.0 * boost::math::constants::pi<Real>() * Sigma)),
          estimated_alpha_root_distance_(
              // observed convergence of distance roots f_alpha().
              boost::math::constants::pi<Real>() / (a - Sigma))
    {
        if(a < this->getSigma()) // check (outer boundary > innter boundary)
        {
            throw std::invalid_argument((boost::format(
                "GreensFunction2DRadAbs: a >= sigma : a=%.16g, sigma=%.16g")
                % a % sigma).str());
        }
        GreensFunction2DRadAbs::clearAlphaTable();
    }
    virtual ~GreensFunction2DRadAbs() {}

    Real geth() const throw() {return this->h;}
    Real geta() const throw() {return this->a;}

    Real getestimated_alpha_root_distance_() const throw()
    {
        return this->estimated_alpha_root_distance_;
    }

    virtual EventKind drawEventType(const Real rnd, const Real t) const;

    virtual Real drawTime (const Real rnd) const;
    virtual Real drawR    (const Real rnd, const Real t) const;
    virtual Real drawTheta(const Real rnd, const Real r, const Real t) const;

    Real f_alpha0(const Real alpha) const;
    Real f_alpha (const Real alpha, const Integer n) const;


    Real p_survival      (const Real t) const;
    Real p_survival_table(const Real t, RealVector& table) const;


    Real leaves(const Real t) const;
    Real leavea(const Real t) const;

    Real p_m(const Integer n, const Real r, const Real t) const;

    Real dp_m_at_a(const Integer m, const Real t) const;

    Real p_m_alpha(const unsigned int n, const unsigned int m,
                   const Real r,         const Real t ) const;

    Real dp_m_alpha_at_a(const unsigned int n, const unsigned int m,
                         const Real t) const;

    const char* getName() const {return "GreensFunction2DRadAbs";}

// -------------------------------------------------------------------------- //
//  methods below are kept public for debugging purpose.
// -------------------------------------------------------------------------- //

    std::string dump() const
    {
        return (boost::format("GreensFunction2DRadAbs: D = %.16g, sigma = %.16g"
                ", a = %.16g, kf = %.16g, r0 = %.16g, h = %.16g") % this->getD()
                % this->getSigma() % this->geta() % this->getkf()
                % this->getr0()    % this->geth()).str();
    }

    void GiveRootInterval      (Real& low, Real& high, const Integer n) const;
    void GiveRootIntervalSimple(Real& low, Real& high, const Integer n,
                                const Real i) const;

    Real getAlphaRoot0(const Real low, const Real high) const;
    Real getAlphaRootN(const Real low, const Real high, const Integer n) const;
    Real getAlphaRoot (const Real high, const Real low, const Integer n) const;

    void decideOnMethod2(size_t n, RealVector::size_type i) const;

    void needToSwitchBackMethod1(size_t n, RealVector::size_type i) const;

    Real getAlpha(size_t n, RealVector::size_type i) const;

    Real p_survival_i(const Real alpha) const;
    Real calc_A_i_0  (const Real alpha) const;
    Real leaves_i    (const Real alpha) const;

    boost::tuple<Real,Real,Real>
    Y0J0J1_constants(const Real alpha, const Real t) const;

//    const Real getAlpha( const size_t n, const RealVector::size_type i ) const;
//    const Real getAlpha0( const RealVector::size_type i ) const;

    Real givePDFTheta(const Real theta, const Real r, const Real t) const;
    Real givePDFR    (const Real r,     const Real t) const;

    void dumpRoots(int n);

protected:

    void clearAlphaTable() const;

    RealVector& getAlphaTable( const size_t n ) const
    {
        return this->alphaTable[n];
    }

    Real p_int_r_table(const Real r, const RealVector& Y0_aAnTable,
                                     const RealVector& J0_aAnTable,
                                     const RealVector& Y0J1J0Y1Table) const;

    Real ip_theta_table(const Real theta, const RealVector& p_nTable) const;

    Real p_survival_i_exp_table(const unsigned int i, const Real t,
                                const RealVector& table ) const;

    Real leavea_i_exp(const unsigned int i, const Real alpha) const;
    Real leaves_i_exp(const unsigned int i, const Real alpha) const;

    Real ip_theta_n(const unsigned int m, const Real theta,
                    const RealVector& p_nTable) const;


    Real p_int_r_i_exp_table(const unsigned int i, const Real r,
                             const RealVector& Y0_aAnTable,
                             const RealVector& J0_aAnTable,
                             const RealVector& Y0J1J0Y1Table) const;

    void createPsurvTable(RealVector& table) const;

    void createY0J0Tables(RealVector& Y0_Table, RealVector& J0_Table,
                          RealVector& Y0J1J0Y1_Table, const Real t) const;

    void makep_mTable(RealVector& p_mTable, const Real r, const Real t) const;

    void makedp_m_at_aTable(RealVector& p_mTable, const Real t) const;

    unsigned int guess_maxi(const Real t) const;

    struct f_alpha0_aux_params
    {
        const GreensFunction2DRadAbs* const gf;
        const Real value;
    };

    static Real f_alpha0_aux_F(const Real alpha,
                               const f_alpha0_aux_params* const params)
    {
        return params->gf->f_alpha0(alpha);
    }

    struct f_alpha_aux_params
    {
        const GreensFunction2DRadAbs* const gf;
        const Integer n;
        Real value;
    };

    static Real f_alpha_aux_F(const Real alpha,
                              const f_alpha_aux_params* const params)
    {
        return params->gf->f_alpha(alpha, params->n);
    }

    struct p_survival_table_params
    {
        const GreensFunction2DRadAbs* const gf;
        // const Real r0;
        RealVector& table;
        const Real rnd;
    };

    static Real p_survival_table_F(const Real t,
                                   const p_survival_table_params* const params)
    {
        return params->rnd - params->gf->p_survival_table(t, params->table);
    }

    struct p_int_r_params
    {
        const GreensFunction2DRadAbs* const gf;
        const Real t;
        // const Real r0;
        const RealVector& Y0_aAnTable;
        const RealVector& J0_aAnTable;
        const RealVector& Y0J1J0Y1Table;
        const Real rnd;
    };

    static Real p_int_r_F(const Real r, const p_int_r_params* const params)
    {
        return params->gf->p_int_r_table(r, params->Y0_aAnTable,
                params->J0_aAnTable, params->Y0J1J0Y1Table) - params->rnd;
    }

    struct ip_theta_params
    {
        const GreensFunction2DRadAbs* const gf;
        const Real r;
        // const Real r0;
        const Real t;
        const RealVector& p_nTable;
        const Real value;
    };

    static Real ip_theta_F(const Real theta,
                           const ip_theta_params* const params)
    {
        const Real one_div_two_pi = boost::math::constants::one_div_two_pi<Real>();
        return theta * one_div_two_pi + params->gf->ip_theta_table(
                theta, params->p_nTable) * one_div_two_pi * 2 - params->value;
    }

  private:
    // Error tolerance used by default.
    static const Real TOLERANCE;

    // SphericalBesselGenerator's accuracy, used by some
    // theta-related calculations.
    static const Real MIN_T_FACTOR;

    static const Real L_TYPICAL; // typical length scale
    static const Real T_TYPICAL; // typical time scale
    static const Real EPSILON; // relative numeric error
    // TESTING temporarily increased; was 1e-12

    static const unsigned int MAX_ORDER = 30; // The maximum number of m terms
    static const unsigned int MAX_ALPHA_SEQ = 500; // The maximum number of n terms

    // Parameters for alpha-root finding
    // ======
    // See getAlpha() in cpp file for more information.
    //
    // Parameters for scanning method
    // Left boundary of 1st search interval 1st root
    static const Real SCAN_START;
    // Length of the scanning interval relative to estimated interval
    static const Real FRACTION_SCAN_INTERVAL; // TODO CHANGED THIS FROM .5 to .2

    // Other paramters
    // After CONVERGENCE_ASSUMED subsequent roots that lay within +/-
    // INTERVAL_MARGIN from the distance to which the distance is known to
    // converge, it is assumed all following roots have a distances inbetween
    // that don't deviate for more than INTERVAL_MARGIN from the distance to
    // which the roots are known to converge (Pi/(a-sigma)).
    static const Real CONVERGENCE_ASSUMED;
    static const Real INTERVAL_MARGIN;

private:

    const Real a;
    const Real h;

    // Tables that hold calculated roots (y=0) of "alpha" function for each
    // order n.
    mutable boost::array<RealVector, MAX_ORDER+1> alphaTable;

    // Constants used in the roots of f_alpha() finding algorithm.
    // ====
    //
    // This constant will simply be M_PI/(a-Sigma), the value to which the
    // distance between roots of f_alpha() should converge.
    const Real estimated_alpha_root_distance_;

    // Table which tells us at which x we're left with scanning the alpha
    // function for a sign change, for a given order n. (A sign change would
    // indicate a root (y=0) lies between the boundaries of the "scanned"
    // interval.)
    //      If x_scan[n] < 0, this indicates scanning is no longer required
    // because the distance between the roots is converging and within
    // boundaries that allow the direct use of the estimate interval width
    // pi/(sigma-a).
    //      Initial values are set by constructor.
    mutable boost::array<Real, MAX_ORDER+1> alpha_x_scan_table_;
    //
    // Table that keeps track of the number of previous subsequent roots that
    // we're within margin of the distance to which they're expected to
    // converge.
    mutable boost::array<int, MAX_ORDER+1> alpha_correctly_estimated_;
};

} // greens_functions
#endif // GREENS_FUNCTION_2D_RAD_ABS_HPP
