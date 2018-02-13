#include "GreensFunction2DAbsSym.hpp"
#include <gsl/gsl_sf_bessel.h>
// #include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/cstdint.hpp>
#include <exception>
#include <stdexcept>
#include <cmath>

namespace greens_functions
{

const Real GreensFunction2DAbsSym::CUTOFF = 1e-10;
const Real GreensFunction2DAbsSym::CUTOFF_H = 6.0;

// an alternative form, which is not very convergent.
Real GreensFunction2DAbsSym::p_survival(const Real t) const
{
    const Integer N_MAX = 100;
    const Real    Dt    = this->D_ * t;

    Real sum(0.0);
    for(Integer n=1; n<=N_MAX; ++n)
    {
        const Real aAn    = gsl_sf_bessel_zero_J0(n);
        const Real An     = aAn / this->a_;
        const Real J1_aAn = gsl_sf_bessel_J1(aAn);
        const Real term   = std::exp(-Dt * An * An) / (An * J1_aAn);
        sum += term;

        if(std::abs(term / sum) < CUTOFF)
        {
            break;
        }
    }
    return (2.0 / a_) * sum;
}

Real GreensFunction2DAbsSym::p_int_r_free(const Real r, const Real t) const
{
    const Real Dt     = this->D_ * t;
    const Real sqrtDt = std::sqrt(Dt);
    const Real sqrtPI = boost::math::constants::root_pi<Real>();

    return boost::math::erf(r / (sqrtDt + sqrtDt)) -
           r * std::exp(-r * r / (4.0 * Dt)) / (sqrtPI * sqrtDt);
}

Real GreensFunction2DAbsSym::p_int_r(const Real r, const Real t) const
{
    const Real Dt = this->D_ * t;
    const Integer N_MAX = 10000;
//    const Real maxn( ( a / M_PI ) * sqrt( log( exp( DtPIsq_asq ) / CUTOFF ) /
//                                          ( D * t ) ) );
    Real sum = 0.0;
    for(Integer n=1; n<=N_MAX; ++n)
    {
        const Real aAn    = gsl_sf_bessel_zero_J0(n);
        const Real An     = aAn / a_;
        const Real rAn    = r * An;
        const Real J1_aAn = gsl_sf_bessel_J1(aAn);
        const Real J1_rAn = gsl_sf_bessel_J1(rAn);
        const Real term   = std::exp(-Dt * An * An) * r * J1_rAn /
                            (An * J1_aAn * J1_aAn);
        sum += term;

        if(std::abs(term / sum) < CUTOFF)
        {
            break;
        }
    }
    return (2.0 / (a_ * a_)) * sum;
}

Real GreensFunction2DAbsSym::drawTime(const Real rnd) const
{
    if(!(0.0 <= rnd && rnd < 1.0))
    {
        throw std::invalid_argument(boost::str(boost::format(
            "GreensFunction2DAbsSym: rnd(%1%) must be in [0, 1)") % rnd));
    }
    if(this->D_ == 0.0 || this->a_ == std::numeric_limits<Real>::infinity())
    {
        return std::numeric_limits<Real>::infinity();
    }
    if(this->a_ == 0.0)
    {
        return 0.0;
    }

    p_survival_equation_t p_survival_eq(*this, rnd);

    // Find a good interval to determine the first passage time in
    // a guess: msd = sqrt(2*d*D*t)
    const Real t_guess = a_ * a_ / (4.0 * D_);
    const Real value   = p_survival_eq(t_guess);

    Real low        = t_guess;
    Real high       = t_guess;
    Real low_value  = value;
    Real high_value = value;

    // scale the interval around the guess such that the function straddles
    if(value < 0.0)
    {
        while(high_value < 0.0)
        {
            high *= 10.0;
            high_value = p_survival_eq(high);
            if(std::abs(high) >= t_guess * 1e6)
            {
                throw std::runtime_error(
                        "GreensFunction2DAbsSym: couldn't adjust high.");
            }
        }
    }
    else
    {
        Real low_value_prev = low_value;
        while(low_value >= 0.0)
        {
            low *= 0.1;
            low_value = p_survival_eq(low);

            if(std::abs(low) <= t_guess * 1e-6 ||
               std::abs(low_value - low_value_prev) < CUTOFF)
            {
                return low;
            }
            low_value_prev = low_value;
        }
    }

    tolerance_t tol(/*absolute = */1e-18, /*relative = */1e-12);
    boost::uintmax_t iter = 100;
    const std::pair<Real, Real> t = boost::math::tools::toms748_solve(
            p_survival_eq, low, high, low_value, high_value, tol, iter);
    if(iter == 100)
    {
        throw std::runtime_error(boost::str(boost::format(
            "GreensFunction2DAbsSym::drawTime: failed to find a root:"
            "rnd = %1%, high = %2%, low = %3%, iter = %4%") %
            rnd % high % low % iter));
    }
    return t.first;
}

Real GreensFunction2DAbsSym::drawR(const Real rnd, const Real t) const
{
    if(!(0.0 <= rnd && rnd < 1.0))
    {
        throw std::invalid_argument(boost::str(boost::format(
            "GreensFunction2DAbsSym::drawR: rnd(%1%) must be in [0,1)") % rnd));
    }
    if(t < 0.0)
    {
        throw std::invalid_argument(boost::str(boost::format(
            "GreensFunction2DAbsSym::drawR: t(%1%) must be positive") % t));
    }
    if(a_ == 0.0 || t == 0.0 || D_ == 0.0)
    {
        return 0.0;
    }
    const Real psurv = this->p_survival(t);
    const p_int_r_equation_t p_int_r_eq(*this, t, psurv * rnd);

    assert(psurv > 0.0);

    const Real low  = 0.0;
    const Real high = this->a_;
    const tolerance_t tol(/*absolute = */1e-18, /*relative = */1e-12);
    boost::uintmax_t  iter = 100;

    const std::pair<Real, Real> r = boost::math::tools::toms748_solve(
            p_int_r_eq, low, high, tol, iter);
    if(iter == 100)
    {
        throw std::runtime_error(boost::str(boost::format(
            "GreensFunction2DAbsSym::drawR: failed to find a root:"
            "rnd = %1%, t = %2%, p_survival = %3%, iter = %4%") %
            rnd % t % psurv % iter));
    }
    return r.first;

// //  if( a <= thresholdDistance )	// if the domain is not so big, the boundaries are felt
// //  {
//         psurv = p_survival( t );
//         //psurv = p_int_r( a, t );
//         //printf("dr %g %g\n",psurv, p_survival( t ));
//         //assert( fabs(psurv - p_int_r( a, t )) < psurv * 1e-8 );
//
//         assert( psurv > 0.0 );
//         F.function = reinterpret_cast<double (*)(double, void*)>( &p_r_F );
// /*  }
//     else				// if the domain is very big, just use the free solution
//     {
//         // p_int_r < p_int_r_free
//         if( p_int_r_free( a, t ) < rnd )	// if the particle is outside the domain?
//         {
//             std::cerr << "p_int_r_free( a, t ) < rnd, returning a." 
//                       << std::endl;
//             return a;
//         }
//
//         psurv = 1.0;
//         F.function = reinterpret_cast<double (*)(double, void*)>( &p_r_free_F );
//     }
// */
//
//     //const Real high( std::min( thresholdDistance, a ) );
//
//     const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
//     gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
//
//     const Real r( findRoot( F, solver, low, high, 1e-18, 1e-12,
//                             "GreensFunction2DAbsSym::drawR" ) );
//     gsl_root_fsolver_free( solver );
//     return r;
}
} // greens_functions
