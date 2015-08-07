#ifndef FACE_INF_HPP
#define FACE_INF_HPP
#include "FaceBase.hpp"

using namespace greens_functions;

// this infinite plane must pass through origin of coordinate
class FaceInf : public FaceBase
{
public:
  FaceInf( const int& id, const Realvec& norm, const Realvec& rep )
  : FaceBase( id, norm, rep )
  {
    ;
  };

  ~FaceInf(){ };

  virtual Realvec renew_position( const Realvec& position, const Realvec& displacement, boost::shared_ptr<FaceBase>& p);

  virtual bool still_in_the_face( const Realvec& position, const Realvec& displacement );

  virtual Realvec get_vertex();
  
  virtual void print_class_name();
};

Realvec FaceInf::renew_position( const Realvec& position, const Realvec& displacement, boost::shared_ptr<FaceBase>& p)
{
    bool in_the_infty_plane;
    in_the_infty_plane = still_in_the_face(position, displacement);
    THROW_UNLESS(std::invalid_argument, in_the_infty_plane);

    Realvec temp(position + displacement);
    return temp;
}

bool FaceInf::still_in_the_face( const Realvec& position, const Realvec& displacement )
{
  Real epsilon(1e-12);
  Realvec v1( get_vertex() );
  Realvec v2( position + displacement );
  Realvec v3( v2 - v1 );
  Realvec norm( get_normal_vector() );

  Real orth( dot_product(v3, norm) );
  
  return ( orth < epsilon );
};

Realvec FaceInf::get_vertex()
{
  Realvec zero(0e0, 0e0, 0e0);
  return zero;
};

void FaceInf::print_class_name()
{
  std::cout << "class: FaceInf" << std::endl;
}

#endif /*FACE_INF_HPP*/
