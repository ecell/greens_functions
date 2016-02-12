#include <iomanip>
#include <cmath>
#include <boost/format.hpp>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_roots.h>
#include "GreensFunction2DRefWedgeAbs.hpp"
#include "findRoot.hpp"

namespace greens_functions
{
    const Real GreensFunction2DRefWedgeAbs::CUTOFF = 1e-10;

    GreensFunction2DRefWedgeAbs::GreensFunction2DRefWedgeAbs(const Real D_,
                                                             const Real r0_,
                                                             const Real a_,
                                                             const Real phi_)
        : D(D_), a(a_), r0(r0_), phi(phi_)
    {
        ;
    }

    GreensFunction2DRefWedgeAbs::~GreensFunction2DRefWedgeAbs()
    {
        ;
    }

    const Real GreensFunction2DRefWedgeAbs::p_survival(const Real t) const
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

        // calculating
        //   2 / a * sum_m[ J0(r0 * alpha_m0) / {alpha_m0 * J1(a * alpha_m0)} *
        //                  exp(-alpha_m0^2 * Dt)]
        // same as 2DAbs.
        for(/*Integer n = 1*/; n < num_term_use; ++n)
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

    const Real GreensFunction2DRefWedgeAbs::p_int_r(const Real r, const Real t) const
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
        // calculating
        //   2r / a^2 *
        //   sum_m[ J0(r0 * alpha_m0)J1(r * alpha_m0) /  <= term 2
        //          {alpha_m0 * {J1(a * alpha_m0)}^2} *  <= term 3
        //          exp(-alpha_m0^2 * Dt)]               <= term 1
        // same as 2DAbs.
        for(/*Integer n = 1*/; n < num_term_use; ++n)
        {
            a_alpha_n = gsl_sf_bessel_zero_J0(n);// n-th root (a * alpha_n)
            alpha_n = a_alpha_n / a;
            J0_r0_alpha_n = gsl_sf_bessel_J0(r_0 * alpha_n);
            J1_r_alpha_n  = gsl_sf_bessel_J1(r * alpha_n);
            J1_a_alpha_n  = gsl_sf_bessel_J1(a_alpha_n);

            // exponential term
            term1 = std::exp(-1e0 * alpha_n * alpha_n * Dt);
            // numerator. product of bessel functions
            term2 = J1_r_alpha_n * J0_r0_alpha_n;
            // denominator. product of alpha and square of bessel functions
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

        return (2e0 * r * sum / (a*a));
    }

    const Real GreensFunction2DRefWedgeAbs::p_int_theta(const Real r,
                                                        const Real theta,
                                                        const Real t) const
    {
        /* XXX: Theta definition isn't same as other GF. The initial position *
         *      is not zero but phi/2. Compensation will be done by drawTheta.*
         *      This function asserts the accepting value theta is in the     *
         *      range [0, phi]                                                */

        // where the r is zero, theta is not defined.
        // return linearly incleasing function(integrated uniform distribution).
        if(fabs(r) < CUTOFF)
        {
            return theta / this->phi;
        }

        // when the initial r equals the abs boundary,
        // particle cannot go anywhere.
        if(fabs(r-a) < CUTOFF)
        {
            return 0e0;
        }

        // this returns accumurate probability distribution,
//         if(theta == phi) return 0e0; // theta == 0e0?
//         if(fabs(theta - 2*M_PI) < CUTOFF) return 1e0;

        const Real first_term(p_int_theta_first(r, theta, t));
        const Real second_term(p_int_theta_second(r, theta, t));

        return (first_term + second_term);
    }

    const Real GreensFunction2DRefWedgeAbs::p_int_theta_first(const Real r,
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

        // calculating
        // 2 * theta / (phi * a^2) *
        // sum_m[ exp(-alpha_m0^2 * Dt) *                <= term 1
        //        J0(r0 * alpha_m0) * J0(r * alpha_m0) / <= term 2
        //        {J1(a * alpha_m0)}^2                   <= term 3

        Integer n(1);
        for(/*Integer n = 1*/; n < num_term_use; ++n)
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

        return (2e0 * theta * sum / (this->phi * a * a));
    }

    const Real
        GreensFunction2DRefWedgeAbs::p_int_theta_second(const Real r,
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

        // prepair (2pi / phi) and (2pi theta / phi). using tau = 2pi.
        const Real tau_phi = 2e0 * M_PI / this->phi;
//         const Real tau_theta_phi = tau_phi * theta;

        // for (-1)^n
        Real sgn = -1e0;

        // calculating 
        // 8 / (pi * a^2) * 
        // sum_n^inf [
        //     ((-1)^n / n) * sin(n * 2 * pi * theta / phi) * 
        //     sum_m [
        //         exp(-alpha_mn^2 * Dt) *
        //         J_{n * 2pi / phi}(r0 * alpha_mn) *
        //         J_{n * 2pi / phi}(r * alpha_mn) / 
        //        (J_{n * 2pi / phi - 1}(a * alpha_mn) -
        //         J_{n * 2pi / phi + 1}(a * alpha_mn))^2 
        //     ]
        // ]

        Integer n(1);
        for(/*unsigned int n = 1*/; n < num_out_term_use; ++n)
        {
            Real in_sum(0e0);
            Real in_term(0e0);
            Real in_term1(0e0);
            Real in_term2(0e0);
            Real in_term3(0e0);

            Real a_alpha_mn(0e0);
            Real alpha_mn(0e0);
            Real Jnpp_r_alpha_mn(0e0);
            Real Jnpp_r0_alpha_mn(0e0);
            Real Jnpp_d_1_a_alpha_mn(0e0);// J_n-1(a alpha_mn)
            Real Jnpp_p_1_a_alpha_mn(0e0);// J_n+1(a alpha_mn)

            Real bessel_order(n * tau_phi);

            Integer m(1);
            for(/*unsigned int m = 1*/; m < num_in_term_use; ++m)
            {
                a_alpha_mn = gsl_sf_bessel_zero_Jnu(bessel_order, m);
                alpha_mn = a_alpha_mn / a;
                Jnpp_r_alpha_mn
                    = gsl_sf_bessel_Jnu(bessel_order, r * alpha_mn);
                Jnpp_r0_alpha_mn
                    = gsl_sf_bessel_Jnu(bessel_order, r_0 * alpha_mn);
                Jnpp_d_1_a_alpha_mn
                    = gsl_sf_bessel_Jnu(bessel_order - 1, a_alpha_mn);
                Jnpp_p_1_a_alpha_mn
                    = gsl_sf_bessel_Jnu(bessel_order + 1, a_alpha_mn);

                in_term1 = std::exp(alpha_mn * alpha_mn * minusDt);
                in_term2 = Jnpp_r_alpha_mn * Jnpp_r0_alpha_mn;
                in_term3 = Jnpp_d_1_a_alpha_mn - Jnpp_p_1_a_alpha_mn;

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
                std::cout << "warning: too slow convergence in p_int_theta_2nd m."
                          << std::endl;

            term = sgn * in_sum * sin(bessel_order * theta) / n;
            sum += term;

            sgn *= -1e0;

//             std::cout << "outer sum " << sum << ", term" << term << std::endl;

            if(fabs(in_sum / (n * sum)) < threshold)
            {
//                 std::cout << "normal exit. n = " << n << " second term" << std::endl;
                break;
                /* if bessel_order * theta = product of natural number and pi,*
                 * term become zero and this code breaks loop too early.      *
                 * Therefore, consider the difference between sum and         *
                 * in_sum / n only. sin is always in range(-1,1). the effect  *
                 * is small.                                                  */
            }
        }
        if(n == num_out_term_use)
            std::cout << "warning: too slow convergence in p_int_theta_2nd n."
                      << std::endl;

        return (8e0 * sum / (M_PI * a * a));
    }

    const Real GreensFunction2DRefWedgeAbs::p_int_phi(const Real r, const Real t) const
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

        // calculating
        // 2 / a^2 sum_m[ exp(-alpha_m0^2 * Dt) *                <= term 1
        //                J0(r0 * alpha_m0) * J0(r * alpha_m0) / <= term 2
        //                {J1(a * alpha_m0)}^2                   <= term 3
        // same as 2DAbs::p_int_2pi.
        Integer n(1);
        for(/*Integer n = 1*/; n < num_term_use; ++n)
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
        if(n == num_term_use)std::cout << "warning: use term over num_term_use" << std::endl;
            std::cout << "warning: use term over num_term_use" << std::endl;

        return (2e0 * sum / (a * a));
    }

    const std::string GreensFunction2DRefWedgeAbs::dump() const
    {
        std::ostringstream ss;
        ss << "D = " << this->getD() << ", a = " << this->geta()
           << ", r_0 = " << this->getr0() << ", phi = " << this->getphi()
           << std::endl;
        return ss.str();
    }

//******************************* drawTime ***********************************//
    const Real
        GreensFunction2DRefWedgeAbs::p_survival_F(const Real t,
                                          const p_survival_params* params)
    {
        const GreensFunction2DRefWedgeAbs* const gf(params->gf);
        const Real rnd(params->rnd);

        // to search the t that satisfies 1 - p_survival(t) = rnd.
        return 1e0 - gf->p_survival(t) - rnd;
    }

    const Real GreensFunction2DRefWedgeAbs::drawTime(const Real rnd) const
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
                              "GreensFunction2DRefWedgeAbs::drawTime"));

        gsl_root_fsolver_free(solver);

        return t;
    }

//********************************* drawR ************************************//
    const Real
        GreensFunction2DRefWedgeAbs::p_r_F(const Real r, const p_r_params* params)
    {
        const GreensFunction2DRefWedgeAbs* const gf(params->gf);
        const Real t(params->t);
        const Real target(params->target);

        return gf->p_int_r(r, t) - target;
    }

    const Real GreensFunction2DRefWedgeAbs::drawR(const Real rnd, const Real t) const
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
                              "GreensFunction2DRefWedgeAbsSym::drawR"));

        gsl_root_fsolver_free(solver);

        return r;
    }


//********************************* drawTheta ********************************//
    const Real GreensFunction2DRefWedgeAbs::p_theta_F(const Real theta,
                                              const p_theta_params* params)
    {
        const GreensFunction2DRefWedgeAbs* const gf(params->gf);
        const Real t(params->t);
        const Real r(params->r);
        const Real rnd(params->rnd);

        return gf->p_int_theta(r, theta, t) - rnd;
    }

    const Real GreensFunction2DRefWedgeAbs::drawTheta(const Real rnd,
                                              const Real r,
                                              const Real t) const
    {
        THROW_UNLESS(std::invalid_argument, 0.0<=rnd && rnd <= 1.0);

        if(rnd == 1e0)
            return phi;

        if(fabs(r) < CUTOFF)// r == 0e0 ?
        {
            throw std::invalid_argument(
                    (boost::format(
                         "2DRefWedgeAbs::drawTheta r is too small: r=%f10"
                         ) % r).str());
        }

        if(fabs(r-a) < CUTOFF)// r == a ?
        {
            //when R equal a, p_int_theta is zero at any theta
            throw std::invalid_argument(
                    (boost::format(
                         "2DAbs::drawTheta r is nealy a: r=%f10, a=%f10"
                         ) % r % a).str());
        }

        if(t == 0e0 || D == 0e0)
            return 0e0;

        Real int_2pi = p_int_phi(r, t);

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
                                  "GreensFunction2DRefWedgeAbsSym::drawTheta"));

        gsl_root_fsolver_free(solver);

        return theta;
    }
}
