#if !defined( _DEFS_HPP )
#define _DEFS_HPP

#include <cstddef>
#include <boost/math/quaternion.hpp>
#include <cmath>
#include <stdexcept>
#include <cassert>
#include <utility>
#include "Vector3.hpp"

typedef double Real;
typedef long int Integer;
typedef unsigned long int UnsignedInteger;
typedef Vector3<double> Realvec;
typedef boost::math::quaternion<double> Realquat;

const Real GLOBAL_TOLERANCE(1e-7);
const int RENEWLOOP_UPPER_LIMIT(100);

#define STR( S ) #S
#define THROW_UNLESS( CLASS, EXPRESSION )\
    if(!(EXPRESSION))\
      throw CLASS("Check ["+std::string(STR(EXPRESSION)) + "] failed.");

// const Real SEPARATION_TOLERANCE( 1e-07  );
// const Real MINIMAL_SEPARATION_FACTOR( 1.0 + SEPARATION_TOLERANCE );

Realvec rotation(const double angle, const Realvec& axis, const Realvec& target)
{
  const double cos_t( cos( angle / 2 ) );
  const double sin_t( sin( angle / 2 ) );
  const double sin_normalize( sin_t / length(axis) );
 
  Realquat Q(cos_t, axis[0] * sin_normalize, 
                    axis[1] * sin_normalize,
                    axis[2] * sin_normalize );
  Realquat P(0e0, target[0], target[1], target[2] );
  Realquat S( Q * P * conj(Q) );
 
  Realvec temp;
  temp[0] = S.R_component_2();
  temp[1] = S.R_component_3();
  temp[2] = S.R_component_4();
  return temp;
}

//**************************** std::pair ****************************************
std::pair<Real, Real>
sum( const std::pair<Real, Real>& lhs, const std::pair<Real, Real>& rhs )
{
  std::pair<Real, Real> retpair( lhs.first + rhs.first, lhs.second + rhs.second );
  return retpair;
}

std::pair<Real, Real>
subtract(const std::pair<Real, Real>& lhs, const std::pair<Real, Real>& rhs)
{
  std::pair<Real, Real> retpair( lhs.first - rhs.first, lhs.second - rhs.second );
  return retpair;
}

std::pair<Real, Real>
multiple( const Real lhs, const std::pair<Real, Real>& rhs )
{
  std::pair<Real, Real> retpair( lhs * rhs.first, lhs * rhs.second );
  return retpair;
}

std::pair<Real, Real>
multiple( const std::pair<Real, Real>& lhs, const Real rhs )
{
  std::pair<Real, Real> retpair( rhs * lhs.first, rhs * lhs.second );
  return retpair;
}


#endif // _DEFS_HPP
