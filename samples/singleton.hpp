#ifndef SINGLETON_HPP
#define SINGLETON_HPP

#include "FaceBase.hpp"
#include "rotation.hpp"

class particle 
{
  Realvec position;
  int particle_id;
  boost::shared_ptr<FaceBase> face_ptr;
  bool shell_involves_vertex;
public:
  particle( int id, const Realvec& pos, boost::shared_ptr<FaceBase> ptr )
  : shell_involves_vertex(false)
  {
    if( fabs( dot_product( pos - ptr->get_vertex(), ptr->get_normal_vector() ) ) > 1e-12 )
    {
      std::cout << "dot_product: " << dot_product( pos - ptr->get_vertex(), ptr->get_normal_vector() ) << std::endl;
      throw std::invalid_argument("particle is not on the plane");
    //TODO if !in_this_face( pos, face_ptr ) then throw invalid_argument
    }

    particle_id = id;
    position = pos;
    face_ptr = ptr;
  };

  void move( const Real& r, const Real& theta );

  Real get_max_a();

  int get_id(){ return particle_id; };

  int get_face_id();

  bool involve_vertex(){ return shell_involves_vertex; };

  boost::shared_ptr<FaceBase> get_face_sptr();

  friend std::ostream& operator<<(std::ostream& os, const particle& part);
};

//TODO make it possible that particle move in 2d and 3d using same function
void particle::move( const Real& r, const Real& theta )
{

  if(shell_involves_vertex)
  {
    //TODO
    //theta?
//     Realvec axis( face_ptr->get_normal_vector() );
//     Realvec target( face_ptr->get_represent_vector() );
//
//     Realvec temp;
//     temp = rotation(theta, axis, target);
//
//     Realvec displacement( temp * r );
//
//     Realvec newpos( face_ptr->renew_position(position, displacement, face_ptr) );
//     position = newpos;
//     return;

  }else{

    Realvec axis( face_ptr->get_normal_vector() );
    Realvec target( face_ptr->get_represent_vector() );
    Realvec disp_direction( rotation(theta, axis, target) );
    Realvec displacement( disp_direction * r );

    if( fabs(length(displacement) - r) > 1e-12 )
    {
      std::cerr << "r: " << std::setprecision(16) << r << std::endl;
      std::cerr << "disp_direction length: " << std::setprecision(16)
		<< length(disp_direction) << std::endl;
      std::cerr << "displacement length: " << std::setprecision(16)
		<< length(displacement) << std::endl;
      std::cerr << "fabs(length(displacement) - r)" << std::setprecision(16) 
		<< fabs(length(displacement) - r) << std::endl;
      throw std::invalid_argument("displacement length is not equal to shell size.");
    }

    std::pair<Realvec, FaceBase_sptr> pos_ptr;
    pos_ptr = face_ptr->apply_displacement(position, displacement, face_ptr);

    position = pos_ptr.first;
    face_ptr = pos_ptr.second;

    return;
  }

}

Real particle::get_max_a()
{
  Real retval( face_ptr->get_max_a(position, shell_involves_vertex) );
  return retval;
}

boost::shared_ptr<FaceBase> particle::get_face_sptr()
{
  return face_ptr;
};

int particle::get_face_id()
{
  int ret_id(face_ptr->get_id());
  return ret_id;
}

std::ostream& operator<<( std::ostream& os, const particle& part)
{
  os << std::setprecision(16) << std::setw(25) << part.position[0] << " ";
  os << std::setprecision(16) << std::setw(25) << part.position[1] << " ";
  os << std::setprecision(16) << std::setw(25) << part.position[2];
  return os;
}

#endif /*SINGLETON_HPP*/
