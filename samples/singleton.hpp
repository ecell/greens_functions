#include "face_base.hpp"
#include <boost/math/quaternion.hpp>

typedef boost::shared_ptr<face_base> face_sptr;
typedef boost::math::quaternion<double> Realquat;
Realvec rotation(const double& angle, const Realvec& axis, const Realvec& target);

class particle 
{
  Realvec position;
  int particle_id;
  face_sptr face_ptr;
public:
  particle( int id, const Realvec& pos, face_sptr ptr )
  {
    //caution: using epsilon for decision of (double) ~= 0
    Realvec vertex( ptr->get_vertex() );
    Real epsilon(1e-10);
    THROW_UNLESS( std::invalid_argument,
		  fabs( dot_product(pos-vertex, ptr->get_normal_vector() ) ) < epsilon );
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

  Realvec temp;
  temp = rotation(theta, axis, target);

  Realvec displacement( temp * r );

  Realvec newpos( face_ptr->move(position, displacement, face_ptr) );
  position = newpos;
  return;
}

std::ostream& operator<<( std::ostream& os, const particle& part)
{
  os << part.position[0] << " ";
  os << part.position[1] << " ";
  os << part.position[2];
  return os;
}

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
}
