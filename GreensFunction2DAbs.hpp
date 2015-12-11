#ifndef GREENS_FUNCTION_2D_ABS
#define GREENS_FUNCTION_2D_ABS
#include <vector>
#include <string>
#include <iostream>
#include <sstream>

#include "Defs.hpp"

namespace greens_functions
{
    class GreensFunction2DAbs
    {
    public:

        GreensFunction2DAbs(const Real D, const Real r0, const Real a);
        virtual ~GreensFunction2DAbs();

//         const Real drawTime(const Real rnd) const;
//         const Real drawR(const Real rnd, const Real t) const;
//         const Real drawTheta(const Real rnd, const Real r, const Real t) const;

        const Real p_survival(const Real t) const;
        const Real p_int_r(const Real r, const Real t) const;
        const Real p_int_theta(const Real r, const Real theta, const Real t) const;
        const Real p_int_2pi(const Real r, const Real t) const;

        const Real getD()  const {return this->D;}
        const Real geta()  const {return this->a;}
        const Real getr0() const {return this->r0;}
        const std::string dump() const;
        const char* getName() const
        {
            return "GreensFunction2DAbs";
        }

//     private:

        //calculate the first term of p_int_theta
        const Real p_int_theta_first(const Real r,
                                     const Real theta,
                                     const Real t) const;

        //calculate the second term of p_int_theta
        const Real p_int_theta_second(const Real r,
                                      const Real theta,
                                      const Real t) const;

    private:

        static const Real CUTOFF;

        const Real D;
        const Real a;
        const Real r0;
    };

}
#endif//GREENS_FUNCTION_2D_ABS
