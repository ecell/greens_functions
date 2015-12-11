#include "GreensFunction2DAbs.hpp"
#include <cmath>
#include <iomanip>
#include <gsl/gsl_sf_bessel.h>

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
        if(fabs(r) < 1e-12)
        {
            // when r=0, theta does not have any meaning.
            throw std::invalid_argument("too small r in p_int_theta");
        }

        if(fabs(r-a) < 1e-12)
        {
            // when r=a, value of this function become 0 in everywhere.
            throw std::invalid_argument("too big r(nearly equal a) in p_int_theta");
        }

        const Real first_term(p_int_theta_first(r, theta, t));
        const Real second_term(p_int_theta_second(r, theta, t));
        const Real denominator(p_int_2pi(r, t));
        return (first_term + second_term) / denominator;
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
                 * value of sin, so this considers only in_sum / n_real.       */

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
}
