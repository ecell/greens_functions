#ifndef GREENS_FUNCTION_2D_ABS_SYM_HPP
#define GREENS_FUNCTION_2D_ABS_SYM_HPP
#include "Defs.hpp"
#include <boost/format.hpp>
#include <string>

namespace greens_functions
{

class GreensFunction2DAbsSym
{
  public:
    GreensFunction2DAbsSym(const Real D, const Real a): D_(D), a_(a) {}
    virtual ~GreensFunction2DAbsSym() {}

    Real D() const throw() {return this->D_;}
    Real a() const throw() {return this->a_;}

    Real drawTime(const Real rnd) const;
    Real drawR   (const Real rnd, const Real t) const;

    Real p_survival  (const Real t) const;
    Real p_int_r     (const Real r, const Real t) const;
    Real p_int_r_free(const Real r, const Real t) const;

    const char* getName() const throw() {return "GreensFunction2DAbsSym";}
    std::string dump() const
    {
        return boost::str(boost::format(
            "GreensFunction2DAbsSym: D = %1%, a = %2%") % this->D_ % this->a_);
    }

  private:

    struct tolerance_t
    {
        const Real absolute_;
        const Real relative_;

        tolerance(const Real abs, const Real rel)
            : absolute_(abs), relative_(rel)
        {}

        bool operator()(const Real lhs, const Real rhs) const throw()
        {
            return std::abs(lhs - rhs)       <= absolute_ ||
                   std::abs(lhs / rhs - 1.0) <= relative_;
        }
    };

    struct p_survival_equation_t
    {
        const GreensFunction2DAbsSym& gf_;
        const Real target_;

        p_survival_equation_t(const GreensFunction2DAbsSym& gf, const Real rnd)
            : gf(gf), target_(1.0 - rnd)
        {}

        Real operator()(const Real t) const
        {
            return target_ - gf_.p_survival(t);
        }
    };

    struct p_int_r_equation_t
    {
        const GreensFunction2DAbsSym& gf_;
        const Real t_;
        const Real target_;

        p_int_r_equation_t(const GreensFunction2DAbsSym& gf,
                           const Real target, const Real t)
            : gf(gf), t_(t) target_(target)
        {}

        Real operator()(const Real r) const
        {
            return gf_.p_int_r(r, t_) - target_;
        }
    };

    struct p_int_r_free_equation_t
    {
        const GreensFunction2DAbsSym& gf_;
        const Real t_;
        const Real target_;

        p_int_r_free_equation_t(const GreensFunction2DAbsSym& gf,
                                const Real target, const Real t)
            : gf(gf), t_(t) target_(target)
        {}

        Real operator()(const Real r) const
        {
            return gf_.p_int_r_free(r, t_) - target_;
        }
    };

  private:
    // H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7,
    // 5.6: ~1e-8, 6.0: ~1e-9
    static const Real CUTOFF_H;
    static const Real CUTOFF;

    const Real D_;
    const Real a_;
};

} // greens_functions 
#endif //GREENS_FUNCTION_2D_ABS_SYM_HPP
