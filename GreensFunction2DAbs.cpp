#include <iomanip>
#include <cmath>
#include <boost/format.hpp>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_roots.h>
#include "GreensFunction2DAbs.hpp"
#include "findRoot.hpp"

namespace greens_functions
{
const Real GreensFunction2DAbs::CUTOFF = 1e-10;

GreensFunction2DAbs::GreensFunction2DAbs(const Real D,
                                         const Real r0,
                                         const Real a)
    : D(D), a(a), r0(r0)
{
    ;
}

GreensFunction2DAbs::~GreensFunction2DAbs()
{
    ;
}

const Real GreensFunction2DAbs::p_survival(const Real t) const
{
    // when t == 0.0, return value become eventually 1.0,
    // but the speed of convergence is too slow.
    if(t == 0.0) return 1.0;

    const Real r_0(this->getr0());
    const Real a(this->geta());
    const Real Dt(this->getD() * t);
    const Integer num_term_use(100);
    const Real threshold(CUTOFF);

    Real sum(0e0);
    Real term1(0e0);
    Real term2(0e0);
    Real term(0e0);

    Real a_alpha_n(0e0);
    Real alpha_n(0e0);
    Real J0_r0_alpha_n(0e0);
    Real J1_a_alpha_n(0e0);

    Integer n(1);
    for(; n < num_term_use; ++n)
    {
        a_alpha_n = gsl_sf_bessel_zero_J0(n);
        alpha_n = a_alpha_n / a;
        J0_r0_alpha_n = gsl_sf_bessel_J0(r_0 * alpha_n);
        J1_a_alpha_n  = gsl_sf_bessel_J1(a_alpha_n);

        term1 = std::exp(-1e0 * alpha_n * alpha_n * Dt) * J0_r0_alpha_n;
        term2 = alpha_n * J1_a_alpha_n;
        term = term1 / term2;
        sum += term;

//             std::cout << "sum " << sum << ", term" << term << std::endl;

        if(fabs(term/sum) < threshold)
        {
//                 std::cout << "normal exit. " << n << std::endl;
            break;
        }
    }
    if(n == num_term_use)
        std::cout << "warning: use term over num_term_use" << std::endl;
    return (2e0 * sum / a);
}

const Real GreensFunction2DAbs::p_int_r(const Real r, const Real t) const
{
    //speed of convergence is too slow
    if(r == 0e0) return 0e0;

    const Real r_0(this->getr0());
    const Real a(this->geta());
    const Real Dt(this->getD() * t);
    const Integer num_term_use(100);
    const Real threshold(CUTOFF);

    Real sum(0e0);
    Real term(0e0);
    Real term1(0e0);
    Real term2(0e0);
    Real term3(0e0);

    Real a_alpha_n(0e0);
    Real alpha_n(0e0);
    Real J0_r0_alpha_n(0e0);
    Real J1_r_alpha_n(0e0);
    Real J1_a_alpha_n(0e0);

    Integer n(1);
    for(; n < num_term_use; ++n)
    {
        a_alpha_n = gsl_sf_bessel_zero_J0(n);
        alpha_n = a_alpha_n / a;
        J0_r0_alpha_n = gsl_sf_bessel_J0(r_0 * alpha_n);
        J1_r_alpha_n  = gsl_sf_bessel_J1(r * alpha_n);
        J1_a_alpha_n  = gsl_sf_bessel_J1(a_alpha_n);

        term1 = std::exp(-1e0 * alpha_n * alpha_n * Dt);
        term2 = r * J1_r_alpha_n * J0_r0_alpha_n;
        term3 = (alpha_n * J1_a_alpha_n * J1_a_alpha_n);

        term = term1 * term2 / term3;
        sum += term;

//             std::cout << "sum " << sum << ", term" << term << std::endl;

        if(fabs(term/sum) < threshold)
        {
//                 std::cout << "normal exit. " << n << std::endl;
            break;
        }
    }
    if(n == num_term_use)
        std::cout << "warning: use term over num_term_use" << std::endl;

    return (2e0 * sum / (a*a));
}

const Real GreensFunction2DAbs::p_int_theta(const Real r,
                                            const Real theta,
                                            const Real t) const
{
    if(fabs(r) < CUTOFF)
    {
        return theta * 0.5 / M_PI;
    }

    if(fabs(r-a) < CUTOFF)
    {
        return 0e0;
    }

    if(theta == 0e0) return 0e0;
//         if(fabs(theta - 2*M_PI) < CUTOFF) return 1e0;

    const Real first_term(p_int_theta_first(r, theta, t));
    const Real second_term(p_int_theta_second(r, theta, t));
//         const Real denominator(p_int_2pi(r, t));

    return (first_term + second_term);
}

const Real GreensFunction2DAbs::p_int_theta_first(const Real r,
                                                  const Real theta,
                                                  const Real t) const
{
    const Real r_0(this->getr0());
    const Real a(this->geta());
    const Real minusDt(-1e0 * this->getD() * t);

    const Integer num_term_use(100);
    const Real threshold(CUTOFF);

    Real sum(0e0);
    Real term(0e0);
    Real term1(0e0);
    Real term2(0e0);
    Real term3(0e0);

    Real a_alpha_n(0e0);
    Real alpha_n(0e0);
    Real J0_r_alpha_n(0e0);
    Real J0_r0_alpha_n(0e0);
    Real J1_a_alpha_n(0e0);

    Integer n(1);
    for(; n < num_term_use; ++n)
    {
        a_alpha_n = gsl_sf_bessel_zero_J0(n);
        alpha_n = a_alpha_n / a;
        J0_r_alpha_n  = gsl_sf_bessel_J0(r * alpha_n);
        J0_r0_alpha_n = gsl_sf_bessel_J0(r_0 * alpha_n);
        J1_a_alpha_n  = gsl_sf_bessel_J1(a_alpha_n);

        term1 = std::exp(alpha_n * alpha_n * minusDt);
        term2 = J0_r_alpha_n * J0_r0_alpha_n;
        term3 = J1_a_alpha_n * J1_a_alpha_n;

        term = term1 * term2 / term3;
        sum += term;

//             std::cout << "sum " << sum << ", term" << term << std::endl;

        if(fabs(term/sum) < threshold)
        {
//                 std::cout << "normal exit. n = " << n << " first term" << std::endl;
            break;
        }
    }
    if(n == num_term_use)
        std::cout << "warning: use term over num_term_use" << std::endl;

//         return (sum / (M_PI * a * a));
    return (theta * sum / (M_PI * a * a));
}

const Real GreensFunction2DAbs::p_int_theta_second(const Real r,
                                                   const Real theta,
                                                   const Real t) const
{
    const Real r_0(this->getr0());
    const Real a(this->geta());
    const Real minusDt(-1e0 * this->getD() * t);

    const Integer num_in_term_use(100);
    const Integer num_out_term_use(100);
    const Real threshold(CUTOFF);

    Real sum(0e0);
    Real term(0e0);
    Integer n(1);
    for(; n < num_out_term_use; ++n)
    {
        Real in_sum(0e0);
        Real in_term(0e0);
        Real in_term1(0e0);
        Real in_term2(0e0);
        Real in_term3(0e0);

        Real a_alpha_mn(0e0);
        Real alpha_mn(0e0);
        Real Jn_r_alpha_mn(0e0);
        Real Jn_r0_alpha_mn(0e0);
        Real Jn_d_1_a_alpha_mn(0e0);// J_n-1(a alpha_mn)
        Real Jn_p_1_a_alpha_mn(0e0);// J_n+1(a alpha_mn)

        Real n_real(static_cast<double>(n));
        int n_int(static_cast<int>(n));
        Integer m(1);

        for(; m < num_in_term_use; ++m)
        {
            a_alpha_mn = gsl_sf_bessel_zero_Jnu(n_real, m);
            alpha_mn = a_alpha_mn / a;
            Jn_r_alpha_mn     = gsl_sf_bessel_Jn(n_int, r * alpha_mn);
            Jn_r0_alpha_mn    = gsl_sf_bessel_Jn(n_int, r_0 * alpha_mn);
            Jn_d_1_a_alpha_mn = gsl_sf_bessel_Jn(n_int - 1, a_alpha_mn);
            Jn_p_1_a_alpha_mn = gsl_sf_bessel_Jn(n_int + 1, a_alpha_mn);

            in_term1 = std::exp(alpha_mn * alpha_mn * minusDt);
            in_term2 = Jn_r_alpha_mn * Jn_r0_alpha_mn;
            in_term3 = Jn_d_1_a_alpha_mn - Jn_p_1_a_alpha_mn;

            in_term = in_term1 * in_term2 / (in_term3 * in_term3);
            in_sum += in_term;

//                 std::cout << "inner sum " << in_sum << ", term" << in_term << std::endl;

            if(fabs(in_term/in_sum) < threshold)
            {
//                     std::cout << "normal exit. m = " << m << " second term" << std::endl;
                break;
            }
        }
        if(m == num_in_term_use)
            std::cout << "warning: use term over num_in_term_use" << std::endl;

//             term = in_sum * std::cos(n_real * theta);
        term = in_sum * std::sin(n_real * theta) / n_real;
        sum += term;

//             std::cout << "outer sum " << sum << ", term" << term << std::endl;

        if(fabs(in_sum / (n_real * sum)) < threshold)
        {
            /* if n * theta is a multiple of \pi, the term may be zero and *
             * term/sum become also zero. this is a problem. sin is in a   *
             * regeon [-1, 1], so the order of term does not depend on     *
             * value of sin, so this considers only (in_sum / n_real).     */

//                 std::cout << "normal exit. n = " << n << " second term" << std::endl;
            break;
        }
    }
    if(n == num_out_term_use)
        std::cout << "warning: use term over num_out_term_use" << std::endl;

    return (8e0 * sum / (M_PI * a * a));
}

const Real GreensFunction2DAbs::p_int_2pi(const Real r, const Real t) const
{
    const Real r_0(this->getr0());
    const Real a(this->geta());
    const Real minusDt(-1e0 * this->getD() * t);
    const Integer num_term_use(100);
    const Real threshold(CUTOFF);

    Real sum(0e0);
    Real term(0e0);
    Real term1(0e0);
    Real term2(0e0);
    Real term3(0e0);

    Real a_alpha_n(0e0);
    Real alpha_n(0e0);
    Real J0_r0_alpha_n(0e0);
    Real J0_r_alpha_n(0e0);
    Real J1_a_alpha_n(0e0);

    Integer n(1);
    for(; n < num_term_use; ++n)
    {
        a_alpha_n = gsl_sf_bessel_zero_J0(n);
        alpha_n = a_alpha_n / a;
        J0_r0_alpha_n = gsl_sf_bessel_J0(r_0 * alpha_n);
        J0_r_alpha_n  = gsl_sf_bessel_J0(r * alpha_n);
        J1_a_alpha_n  = gsl_sf_bessel_J1(a_alpha_n);

        term1 = std::exp(alpha_n * alpha_n * minusDt);
        term2 = J0_r_alpha_n * J0_r0_alpha_n;
        term3 = J1_a_alpha_n * J1_a_alpha_n;
        term  = term1 * term2 / term3;
        sum  += term;

//             std::cout << "sum " << sum << ", term" << term << std::endl;

        if(fabs(term/sum) < threshold)
        {
//                 std::cout << "normal exit. n = " << n << " denominator" << std::endl;
            break;
        }
    }
    if(n == num_term_use)
        std::cout << "warning: use term over num_term_use" << std::endl;

    return (2e0 * sum / (a * a));
}

const std::string GreensFunction2DAbs::dump() const
{
    std::ostringstream ss;
    ss << "D = " << this->getD() << ", a = " << this->geta()
       << ", = " << this->getr0() << std::endl;
    return ss.str();
}

//******************************* drawTime ***********************************//
const Real
    GreensFunction2DAbs::p_survival_F(const Real t,
                                      const p_survival_params* params)
{
    const GreensFunction2DAbs* const gf(params->gf);
    const Real rnd(params->rnd);

    //seek certain t that satisfies 1 - p_survival(t) = rnd.
    return 1e0 - gf->p_survival(t) - rnd;
}

const Real GreensFunction2DAbs::drawTime(const Real rnd) const
{
    THROW_UNLESS(std::invalid_argument, 0.0<=rnd && rnd <= 1.0);
    if(D == 0e0 || a == INFINITY || rnd == 1e0)
        return INFINITY;
    if(a == r0 || rnd == 0e0)
        return 0e0;

    p_survival_params params = {this, rnd};

    gsl_function F = 
    {
        reinterpret_cast<typeof(F.function)>(&p_survival_F), &params
    };

    // this is not so accurate because
    // initial position is not the center of this system.
    const Real t_guess(a * a * 0.25 / D);
    Real value(GSL_FN_EVAL(&F, t_guess));
    Real low(t_guess);
    Real high(t_guess);

    // to determine high and low border
    if(value < 0.0)
    {
        do
        {
            high *= 1e1;
            value = GSL_FN_EVAL(&F, high);
            if(fabs(high) > t_guess * 1e6)
                throw std::invalid_argument("could not adjust higher border");
        }
        while(value <= 0e0);
    }
    else
    {
        Real value_prev = value;
        do
        {
            low *= 1e-1;
            value = GSL_FN_EVAL(&F, low);
            if(fabs(low) <= t_guess * 1e-6 || fabs(value - value_prev) < CUTOFF)
                throw std::invalid_argument("could not adjust lower border");
            value_prev = value;
        }
        while(value >= 0e0);
    }

    //find the root
    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));

    const Real t(findRoot(F, solver, low, high, 1e-18, 1e-12,
                          "GreensFunction2DAbs::drawTime"));

    gsl_root_fsolver_free(solver);

    return t;
}

//********************************* drawR ************************************//
const Real
    GreensFunction2DAbs::p_r_F(const Real r, const p_r_params* params)
{
    const GreensFunction2DAbs* const gf(params->gf);
    const Real t(params->t);
    const Real target(params->target);

    return gf->p_int_r(r, t) - target;
}

const Real GreensFunction2DAbs::drawR(const Real rnd, const Real t) const
{
    THROW_UNLESS(std::invalid_argument, 0.0<=rnd && rnd <= 1.0);

    if(a == r0)
        throw std::invalid_argument("a equal r0");

    if(t == 0e0 || D == 0e0)
        return r0;

    if(rnd == 1e0)
        return a;//!?

    Real p_surv(p_survival(t));
    assert(p_surv > 0e0);

    p_r_params params = {this, t, rnd * p_surv};

    gsl_function F =
    {
        reinterpret_cast<typeof(F.function)>(&p_r_F), &params
    };

    const Real low(0e0);
    const Real high(a);

    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));

    const Real r(findRoot(F, solver, low, high, 1e-18, 1e-12,
                          "GreensFunction2DAbsSym::drawR"));

    gsl_root_fsolver_free(solver);

    return r;
}


//********************************* drawTheta ********************************//
const Real GreensFunction2DAbs::p_theta_F(const Real theta,
                                          const p_theta_params* params)
{
    const GreensFunction2DAbs* const gf(params->gf);
    const Real t(params->t);
    const Real r(params->r);
    const Real rnd(params->rnd);

    return gf->p_int_theta(r, theta, t) - rnd;
}

const Real GreensFunction2DAbs::drawTheta(const Real rnd,
                                          const Real r,
                                          const Real t) const
{
    THROW_UNLESS(std::invalid_argument, 0.0<=rnd && rnd <= 1.0);

    if(rnd == 1e0)
        return 2e0 * M_PI;

    if(fabs(r) < CUTOFF)// r == 0e0 ?
    {
        throw std::invalid_argument(
                (boost::format("2DAbs::drawTheta r is too small: r=%f10") % r).str());
    }

    if(fabs(r-a) < CUTOFF)// r == a ?
    {
        //when R equal a, p_int_theta is zero at any theta
        throw std::invalid_argument(
                (boost::format("2DAbs::drawTheta r is nealy a: r=%f10, a=%f10") % r % a).str());
    }

    if(t == 0e0 || D == 0e0)
        return 0e0;

    Real int_2pi = p_int_2pi(r, t);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
     * When t is too large, int_2pi become zero and drawR returns 2pi    *
     * at any value of rnd. To avoid this,  return rnd * theta / 2pi     *
     * because when t -> \infty the probability density function of theta*
     * become uniform distribution                                       *
     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    if(int_2pi == 0e0)
    {
        std::cout << dump();
        std::cout << "Warning: t is too large. t = " << t << std::endl;
    }

    p_theta_params params = {this, t, r, rnd * int_2pi};

    gsl_function F =
    {
        reinterpret_cast<typeof(F.function)>(&p_theta_F), &params
    };

    const Real low(0e0);
    const Real high(2e0 * M_PI);

    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));

    const Real theta(findRoot(F, solver, low, high, 1e-18, 1e-12,
                              "GreensFunction2DAbsSym::drawTheta"));

    gsl_root_fsolver_free(solver);

    return theta;
}
}
