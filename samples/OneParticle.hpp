#ifndef SINGLETON_HPP
#define SINGLETON_HPP
#include <iostream>
#include <stdexcept>

#include "FaceBase.hpp"
#include "rotation.hpp"

class particle 
{
  bool shell_involves_vertex;
  int particle_id;
  Realvec position;
  FaceBase_sptr face_ptr;

public:
  particle( int id, const Realvec& pos, FaceBase_sptr ptr )
  : shell_involves_vertex(false)
  {
    if( fabs( dot_product( pos - ptr->get_para_origin(), ptr->get_normal_vector() ) ) > 1e-12 )
    {
      std::cout << "dot_product: " << dot_product( pos - ptr->get_para_origin(), ptr->get_normal_vector() ) << std::endl;
      throw std::invalid_argument("particle is not on the plane");
    //TODO if !in_this_face( pos, face_ptr ) then throw invalid_argument
    }

    particle_id = id;
    position = pos;
    face_ptr = ptr;
  };

  void move( const Real r, const Real theta );

  Real get_max_a();

  int get_id(){ return particle_id; };

  int get_face_id();

  bool involve_vertex(){ return shell_involves_vertex; };

  FaceBase_sptr get_face_sptr();

  friend std::ostream& operator<<(std::ostream& os, const particle& part);
};

//TODO make it possible that particle move in 2d and 3d using same function
void particle::move( const Real r, const Real theta )
{

  if(shell_involves_vertex)
  {
    //TODO
    //theta?
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

FaceBase_sptr particle::get_face_sptr()
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
