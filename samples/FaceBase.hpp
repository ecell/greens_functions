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
  virtual Realvec renew_position( const Real& pos_alpha, const Real& pos_beta,
				  const Real& dis_alpha, const Real& dis_beta, boost::shared_ptr<FaceBase>& p)
  {
    print_class_name();
    throw std::invalid_argument("this class doesnt have this function or still not overloaded");
    Realvec zero;
    return zero;
   
  }

  virtual bool still_in_the_face( const Realvec& position, const Realvec& displacement ) = 0;

  virtual Realvec get_vertex() = 0;

  virtual void set_belonging_polygon( boost::shared_ptr<Polygon> p_sptr ){};//do nothing

  virtual void set_near_vertexs()
  {
    print_class_name();
    throw std::invalid_argument( "this class has no neighbor or still not overloaded" );
  }

//   virtual void set_neighbors_edge(){};
//   {
//     print_class_name();
//     throw std::invalid_argument( "this class has no neighbor or still not overloaded" )
//   }

  virtual Realvec get_another_vertex( const Realvec& edge )
  {
    print_class_name();
    throw std::invalid_argument( "this class has no neighbor or still not overloaded" );
    Realvec zero;
    return zero;
  };

  virtual Real get_max_a(const Realvec& position, bool& vertex_involve_flag)
  {
    print_class_name();
    throw std::invalid_argument( "this class has no neighbor or still not overloaded" );
    return 0e0;
  };

  virtual Real get_minimum_height(const Realvec& neighbors_edge)
  {
    print_class_name();
    throw std::invalid_argument( "this class has no neighbor or still not overloaded" );
    return 0e0;
  };

  virtual Real get_left_angle( const Realvec& neighbors_edge )
  {
    print_class_name();
    throw std::invalid_argument( "this class has no angle or still not overloaded" );
    return 0e0;
  }

  virtual Real get_right_angle( const Realvec& neighbors_edge )
  {
    print_class_name();
    throw std::invalid_argument( "this class has no angle or still not overloaded" );
    return 0e0;
  }

  virtual Realvec get_para_a()
  {
    print_class_name();
    throw std::invalid_argument( "this class has no neighbor or still not overloaded" );
    Realvec zero;
    return zero;
  };

  virtual Realvec get_para_b()
  {
    print_class_name();
    throw std::invalid_argument( "this class has no neighbor or still not overloaded" );
    Realvec zero;
    return zero;
  };

  virtual void print_class_name()
  {
    std::cout << "class: FaceBase" << std::endl;
    return;
  }

  int get_id(){ return face_id; };
  Realvec get_normal_vector(){ return normal; };
  Realvec get_represent_vector(){ return represent; };

  Real smaller_angle(const Realvec& v1, const Realvec& v2);
};

Real FaceBase::smaller_angle(const Realvec& v1, const Realvec& v2)
{
  Real len1( length(v1) );
  Real len2( length(v2) );
  Real inner( dot_product( v1, v2 ) );

  Real angle( acos( inner / len1 / len2 ) );
  return angle;
}

#endif /*FACE_BASE_HPP*/
