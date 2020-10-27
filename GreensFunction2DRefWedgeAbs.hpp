#ifndef GREENS_FUNCTION_2D_REF_WEDGE_ABS_HPP
#define GREENS_FUNCTION_2D_REF_WEDGE_ABS_HPP
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "Defs.hpp"

namespace greens_functions
{
    class GreensFunction2DRefWedgeAbs
    {
    public:

        GreensFunction2DRefWedgeAbs(const Real D, const Real r0, const Real a,
                                    const Real phi);
        virtual ~GreensFunction2DRefWedgeAbs();

        const Real drawTime (const Real rnd) const;
        const Real drawR    (const Real rnd, const Real t) const;
        const Real drawTheta(const Real rnd, const Real r, const Real t) const;

        const Real p_survival  (const Real t) const;
        const Real p_int_r     (const Real r, const Real t) const;
        const Real p_int_phi   (const Real r, const Real t) const;
        const Real dp_int_phi  (const Real t) const;
        const Real p_int_theta (const Real r, const Real theta, const Real t) const;
        const Real dp_int_theta(const Real theta, const Real t) const;

        const Real getD()   const {return this->D_;  }
        const Real geta()   const {return this->a_;  }
        const Real getr0()  const {return this->r0_; }
        const Real getphi() const {return this->phi_;}
        const std::string dump() const;
        const char* getName() const
        {
            return "GreensFunction2DRefWedgeAbs";
        }

    private:

        const Real p_int_theta_first(const Real r,
                                     const Real theta,
                                     const Real t) const;

        const Real p_int_theta_second(const Real r,
                                      const Real theta,
                                      const Real t) const;

        const Real dp_int_theta_first(const Real theta,
                                      const Real t) const;

        const Real dp_int_theta_second(const Real theta,
                                       const Real t) const;


    private:

        struct p_survival_params
        {
            const GreensFunction2DRefWedgeAbs* const gf;
            const Real rnd;
        };

        struct p_r_params
        {
            const GreensFunction2DRefWedgeAbs* const gf;
            const Real t;
            const Real target;
            /* when time = t, probability of existence is p_surv(t).         *
             * to seek the r, the random value that satisfy the equation     *
             * must be in [0, p_surv). the value, rnd * p_surv(t) is target. */
        };

        struct p_theta_params
        {
            const GreensFunction2DRefWedgeAbs* const gf;
            const Real t;
            const Real r;
            const Real rnd;
        };

        struct dp_theta_params
        {
            const GreensFunction2DRefWedgeAbs* const gf;
            const Real t;
            const Real rnd;
        };

        static const Real p_survival_F(const Real t,
                                       const p_survival_params* params);

        static const Real p_r_F(const Real r,
                                const p_r_params* params);

        static const Real p_theta_F(const Real theta,
                                    const p_theta_params* params);

        static const Real dp_theta_F(const Real theta,
                                     const dp_theta_params* params);


    private:

        static const Real CUTOFF;//1e-10, provisionally.
        static const Real maximum_alpha2_Dt;

        const Real D_;
        const Real a_;
        const Real r0_;
        const Real phi_;// the angle of apical. initial theta = \phi / 2
    };
}
#endif//GREENS_FUNCTION_2D_REF_WEDGE_ABS
