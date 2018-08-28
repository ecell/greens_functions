#if !defined( COMPAT_HPP )
#define COMPAT_HPP

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include <cmath>
#include <limits>

#if !HAVE_DECL_INFINITY
#  ifndef INFINITY
#    define INFINITY ( std::numeric_limits< double >::infinity() )
#  endif
#endif /* HAVE_DECL_INFINITY */

#if !defined( HAVE_SINCOS )
namespace greens_functions
{
inline void sincos( double x, double* s, double* c )
{
    *s = sin( x );
    *c = cos( x );
}
} // greens_functions
#endif /* !HAVE_SINCOS */

#if !defined( HAVE_ISFINITE )
namespace greens_functions
{
inline int isfinite( double x )
{
	return x == x && x != INFINITY && -x != INFINITY;
}
} // greens_functions
#else // std::isfinite found (compiler supports c++11!)
namespace greens_functions
{
inline int isfinite( double x )
{
	return static_cast<int>(std::isfinite(x));
}
} // greens_functions
#endif

#endif // COMPAT_HPP
