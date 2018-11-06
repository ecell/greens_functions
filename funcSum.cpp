#include <vector>
#include <cmath>
#include <boost/bind.hpp>
#include <gsl/gsl_sum.h>

//#include "Logger.hpp"
#include "funcSum.hpp"

namespace greens_functions
{

//static Logger& _log(Logger::get_logger("funcSum"));

Real
funcSum_all(boost::function<Real(unsigned int i)> f, std::size_t max_i)
{
    const Real p_0(f(0));
    if (p_0 == 0.0)
    {
        return 0.0;
    }

    Real sum = p_0;
    for(std::size_t i=1; i < max_i; ++i)
    {
        sum += f(i);
    }
    return sum;
}


Real
funcSum_all_accel(boost::function<Real(unsigned int i)> f,
                  std::size_t max_i, Real tolerance)
{
    const Real p_0(f(0));
    if (p_0 == 0.0)
    {
        return 0.0;
    }

    std::vector<Real> pTable(max_i);
    pTable[0] = p_0;
    for(std::size_t i=1; i < max_i; ++i)
    {
        pTable[i] = f(i);
    }

    Real sum;
    Real error;
    gsl_sum_levin_utrunc_workspace* workspace(gsl_sum_levin_utrunc_alloc(max_i));
    gsl_sum_levin_utrunc_accel(
            pTable.data(), pTable.size(), workspace, &sum, &error);
    if (std::abs(error) >= std::abs(sum * tolerance))
    {
/*        _log.error("series acceleration error: %.16g"
                  " (rel error: %.16g), terms_used = %d (%d given)",
                  std::abs(error), std::abs(error / sum),
                  workspace->terms_used, pTable.size());
        // TODO look into this crashing behaviour */
    }
    gsl_sum_levin_utrunc_free(workspace);

    return sum;
}


Real
funcSum(boost::function<Real(unsigned int i)> f, std::size_t max_i, Real tolerance)
// funcSum
// ==
// Will simply calculate the sum over a certain function f, until it converges
// (i.e. the sum > tolerance*current_term for a CONVERGENCE_CHECK number of
// terms), or a maximum number of terms is summed (usually 2000).
//
// Input:
// - f: A       function object
// - max_i:     maximum number of terms it will evaluate
// - tolerance: convergence condition, default value 1-e8 (see .hpp)
//
// About Boost::function
// ===
// Boost::function doesn't do any type checking: It will take any object
// and any signature you provide in its template parameter, and create an
// object that's callable according to your signature and calls the object.
// If that's impossible, it's a compile error.
// (From: http://stackoverflow.com/questions/527413/how-boostfunction-and-boostbind-work)
{
    // DEFAULT = 4
    const unsigned int CONVERGENCE_CHECK(4);

    const Real p_0(f(0));
    if (p_0 == 0.0)
    {
        return 0.0;
    }

    Real sum(p_0);
    std::vector<Real> pTable;
    pTable.push_back(p_0);

    bool extrapolationNeeded(true);

    unsigned int convergenceCounter(0);

    for(std::size_t i=1; i < max_i; ++i)
    {
        const Real p_i(f(i));
        pTable.push_back(p_i);
        sum += p_i;

        if (std::abs(sum) * tolerance >= std::abs(p_i)) // '=' is important
        {
            ++convergenceCounter;
        }
        // this screws it up; why?
        else
        {
            convergenceCounter = 0;
        }

        if (convergenceCounter >= CONVERGENCE_CHECK)
        {
            extrapolationNeeded = false;
            break;
        }
    }

    if (extrapolationNeeded)
    {
        Real error;
        gsl_sum_levin_utrunc_workspace*
            workspace(gsl_sum_levin_utrunc_alloc(max_i));
        gsl_sum_levin_utrunc_accel(
            pTable.data(), pTable.size(), workspace, &sum, &error);
        if (std::abs(error) >= std::abs(sum * tolerance * 10))
        {
/*            _log.error("series acceleration error: %.16g"
                      " (rel error: %.16g), terms_used = %d (%d given)",
                      std::abs(error), std::abs(error / sum),
                      workspace->terms_used, pTable.size()); */
        }

        gsl_sum_levin_utrunc_free(workspace);
    }

    return sum;
}

}
