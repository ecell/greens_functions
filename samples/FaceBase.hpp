#ifndef FACE_BASE_HPP
#define FACE_BASE_HPP
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include "../Defs.hpp"
#include "Vector3.hpp"

using namespace greens_functions;
class Polygon;
typedef Vector3<Real> Realvec;

class FaceBase
{
protected:
  int face_id;
  Realvec normal;
  Realvec represent;

public:
  FaceBase( const int& id, const Realvec& norm, const Realvec& rep )
  : face_id(id), normal( norm/length(norm) ), represent( rep/length(rep) )
  {
    THROW_UNLESS( std::invalid_argument, dot_product(norm, rep) == 0 );
    THROW_UNLESS( std::invalid_argument, length( norm ) != 0 );
    THROW_UNLESS( std::invalid_argument, length( rep ) != 0 );
  };

  virtual Realvec renew_position( const Realvec& position, const Realvec& displacement, boost::shared_ptr<FaceBase>& p) = 0;

  virtual bool still_in_the_face( const Realvec& position, const Realvec& displacement ) = 0;

  virtual Realvec get_vertex() = 0;

  virtual void set_belonging_polygon( boost::shared_ptr<Polygon> p_sptr ){};//do nothing

  virtual void set_near_vertexs(){};//do_nothing

  virtual Realvec get_another_vertex( const Realvec& edge )
  {
    throw std::invalid_argument( "this face must have no neighbor" );
    Realvec zero;
    return zero;
  };

  virtual Real get_max_a(const Realvec& position)
  {
    throw std::invalid_argument( "this face has no vertex" );
    return 0e0;
  };

  int get_id(){ return face_id; };
  Realvec get_normal_vector(){ return normal; };
  Realvec get_represent_vector(){ return represent; };
};

#endif /*FACE_BASE_HPP*/
