#ifndef ECELL_GREENS_FUNCTIONS_COMPAT_HPP
#define ECELL_GREENS_FUNCTIONS_COMPAT_HPP

#include <greens_functions/config.h>
#include <cmath>
#include <limits>

#ifndef ECELL_GREENS_FUNCTIONS_HAVE_SINCOS
namespace greens_functions
{
inline void sincos( double x, double* s, double* c )
{
    *s = std::sin(x);
    *c = std::cos(x);
}
} // greens_functions
#endif /* HAVE_SINCOS */

#if __cplusplus >= 201103L

// after C++11, std::isfinite is defined in <cmath> header.
namespace greens_functions
{

inline bool isfinite(float x)
{
    return std::isfinite(x);
}
inline bool isfinite(double x)
{
    return std::isfinite(x);
}
inline bool isfinite(long double x)
{
    return std::isfinite(x);
}

} // greens_functions
#else

// it means that the compiler still does not compatible with C++11 or
// the user compiles this without -std=c++(11|14|17) flag.
namespace greens_functions
{

namespace detail
{
template<typename T>
inline bool isfinite_impl(T x)
{
    // return (not nan && not +inf && not -inf);
    // if x was nan, all the comparison including `x == x` would be false!
    return (x == x) && (x !=  std::numeric_limits<T>::infinity()) &&
                       (x != -std::numeric_limits<T>::infinity());
}
} // detail

inline bool isfinite(float x)
{
    return detail::isfinite_impl(x);
}
inline bool isfinite(double x)
{
    return detail::isfinite_impl(x);
}
inline bool isfinite(long double x)
{
    return detail::isfinite_impl(x);
}

} // greens_functions
#endif // ISFINITE
#endif // ECELL_GREENS_FUNCTIONS_COMPAT_HPP
