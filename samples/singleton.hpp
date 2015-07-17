#include "FaceBase.hpp"
#include "rotation.hpp"
typedef boost::shared_ptr<FaceBase> face_sptr;

class particle 
{
  Realvec position;
  int particle_id;
  boost::shared_ptr<FaceBase> face_ptr;
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

  Real get_max_a();

  int get_id(){ return particle_id; };

  boost::shared_ptr<FaceBase> get_face_sptr();

  friend std::ostream& operator<<(std::ostream& os, const particle& part);
};

//TODO make it possible that particle move in 2d and 3d using same function
void particle::move( const Real& r, const Real& theta )
{
  Realvec axis( face_ptr->get_normal_vector() );
  Realvec target( face_ptr->get_represent_vector() );

  Realvec temp;
  temp = rotation(theta, axis, target);

  Realvec displacement( temp * r );

//   std::cout << "call face that have this id: " << face_ptr->get_id() << std::endl;

  Realvec newpos( face_ptr->renew_position(position, displacement, face_ptr) );
  position = newpos;
  return;
}

Real particle::get_max_a()
{
  Real retval( face_ptr->get_max_a(position) );
  return retval;
}

boost::shared_ptr<FaceBase> particle::get_face_sptr(){ return face_ptr; };

std::ostream& operator<<( std::ostream& os, const particle& part)
{
  os << part.position[0] << " ";
  os << part.position[1] << " ";
  os << part.position[2];
  return os;
}

