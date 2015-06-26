#ifndef FACE_BASE_HPP
#define FACE_BASE_HPP
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include "../Defs.hpp"
#include "Vector3.hpp"

using namespace greens_functions;

typedef Vector3<Real> Realvec;

class face_base
{
protected:
  int face_id;
  Realvec normal;
  Realvec represent;

public:
  face_base( const int& id, const Realvec& norm,
	     const Realvec& rep )
  {
    THROW_UNLESS( std::invalid_argument,
		  dot_product(norm, rep) == 0 );
    THROW_UNLESS( std::invalid_argument,
		  length( norm ) != 0);
    THROW_UNLESS( std::invalid_argument,
		  length( rep ) != 0);
    
    face_id = id;
    normal = norm / length( norm );
    represent = rep / length( rep );
  };

  virtual Realvec move(Realvec& position,
		       Realvec& displacement,
		       boost::shared_ptr<face_base>& p)
  {
//     std::cout << "face_base->move" << std::endl;
    Realvec zero;
    return zero;
  };

  virtual bool still_in_the_face( const Realvec& position,
				  const Realvec& displacement )
  { return false; };

  int get_id(){ return face_id; };
  Realvec get_normal_vector(){ return normal; };
  Realvec get_represent_vector(){ return represent; };
  virtual Realvec get_vertex()
  {
    Realvec zero;
    return zero;
  };
};

#endif /*FACE_BASE_HPP*/
