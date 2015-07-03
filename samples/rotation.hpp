#ifndef ROTATION_HPP
#define ROTATION_HPP
#include <boost/math/quaternion.hpp>
#include <cmath>
#include "Vector3.hpp"
typedef boost::math::quaternion<double> Realquat;
// typedef Vector3<double> Realvec;

Realvec rotation(const double& angle, const Realvec& axis, const Realvec& target)
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
};

#endif /*ROTATION_HPP*/
