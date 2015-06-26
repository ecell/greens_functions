#include "face_base.hpp"
#include <boost/math/quaternion.hpp>

typedef boost::shared_ptr<face_base> face_sptr;
typedef boost::math::quaternion<double> Realquat;
Realvec rotation(const double& angle, const Realvec& axis, const Realvec& target);

class particle 
{
  Realvec position;
  int particle_id;
  boost::shared_ptr<face_base> face_ptr;
public:
  particle( int id, const Realvec& pos,
	    boost::shared_ptr<face_base> ptr )
  {
    particle_id = id;
    position = pos;
    face_ptr = ptr;
  };

  void move( const Real& r, const Real& theta );

  int get_id(){ return particle_id; };

  friend std::ostream& operator<<(std::ostream& os, const particle& part);
};

void particle::move( const Real& r, const Real& theta )
{
  Realvec axis( face_ptr->get_normal_vector() );
  Realvec target( face_ptr->get_represent_vector() );
  Realvec displacement(3);

  displacement = r * rotation(theta, axis, target);

  position = face_ptr->move(position, displacement, face_ptr);
  return;
}

std::ostream& operator<<( std::ostream& os, const particle& part)
{
  os << part.position[0] << " ";
  os << part.position[1] << " ";
  os << part.position[2] << " ";
  return os;
}

Realvec rotation(const double& angle, const Realvec& axis, const Realvec& target)
{
  const double cos_t( cos( angle / 2 ) );
  const double sin_t( sin( angle / 2 ) );
  const double normalize( sqrt( axis[0] * axis[0] +
                                axis[1] * axis[1] +
                                axis[2] * axis[2] ) );
  const double sin_normalize( sin_t / normalize );
 
  Realquat Q(cos_t, axis[0] * sin_normalize, 
                    axis[1] * sin_normalize,
                    axis[2] * sin_normalize );
  Realquat P(0e0, target[0], target[1], target[2] );
  Realquat S( Q * P * conj(Q) );
 
  Realvec temp(3);
  temp[0] = S.R_component_2();
  temp[1] = S.R_component_3();
  temp[2] = S.R_component_4();
  return temp;
}
