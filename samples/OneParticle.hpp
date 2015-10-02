#ifndef SINGLETON_HPP
#define SINGLETON_HPP
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <utility>
#include "FaceBase.hpp"
#include "Defs.hpp"

class OneParticle 
{
  bool shell_includes_vertex;
  int particle_id;
  Realvec position;
  FaceBase_sptr face_ptr;
  int times_shell_include_vertex;

public:
  OneParticle( int id, const Realvec& pos, FaceBase_sptr ptr )
  : shell_includes_vertex(false), times_shell_include_vertex(0)
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

  ~OneParticle()
  {
    std::cout << "particle id: " << particle_id << " : number of times shell includes vertex: " << times_shell_include_vertex << std::endl;
  }

  void move( const Real r, const Real theta );

  Real get_max_a();

  int get_id(){ return particle_id; };

  int get_face_id();

  bool include_vertex(){ return shell_includes_vertex; };

  FaceBase_sptr get_face_sptr();

  friend std::ostream& operator<<(std::ostream& os, const OneParticle& part);
};

//TODO make it possible that particle move in 2d and 3d using same function
void OneParticle::move( const Real r, const Real theta )
{

  if(shell_includes_vertex)
  {
    //TODO
    ++times_shell_include_vertex;
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

// std::pair<Realvec, FaceBase_sptr>

Real OneParticle::get_max_a()
{
  Real retval( face_ptr->get_max_a(position, shell_includes_vertex) );
  return retval;
}

FaceBase_sptr OneParticle::get_face_sptr()
{
  return face_ptr;
};

int OneParticle::get_face_id()
{
  int ret_id(face_ptr->get_id());
  return ret_id;
}

std::ostream& operator<<( std::ostream& os, const OneParticle& part)
{
  os << std::setprecision(16) << std::setw(25) << part.position[0] << " ";
  os << std::setprecision(16) << std::setw(25) << part.position[1] << " ";
  os << std::setprecision(16) << std::setw(25) << part.position[2];
  return os;
}

#endif /*SINGLETON_HPP*/
