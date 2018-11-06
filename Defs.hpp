#if !defined( __DEFS_HPP )
#define __DEFS_HPP

#include <cstddef>
#include <limits>
#include <cmath>

namespace greens_functions{

typedef double Real;
typedef long int Integer;
typedef unsigned long int UnsignedInteger;

#define STR( S ) #S
#define THROW_UNLESS( CLASS, EXPRESSION )\
    if(!(EXPRESSION))\
      throw CLASS("Check ["+std::string(STR(EXPRESSION)) + "] failed.");

const Real SEPARATION_TOLERANCE(1e-07);
const Real MINIMAL_SEPARATION_FACTOR(1.0 + SEPARATION_TOLERANCE);

}
#endif // __DEFS_HPP
