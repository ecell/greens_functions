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

using namespace greens_functions;

typedef boost::numeric::ublas::vector<Real> Realvec;

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
		  boost::numeric::ublas::inner_prod(norm, rep) == 0 );
    face_id = id;
    normal = norm;
    represent = rep;
  };

  virtual Realvec move(Realvec& position,
		       const Realvec& displacement,
		       boost::shared_ptr<face_base>& p)
  {
    Realvec zero(3);
      zero[0] = 0e0;
      zero[1] = 0e0;
      zero[2] = 0e0;

    std::cout << "face_base->move" << std::endl;
    return zero;
  };

  virtual bool still_in_the_face( const Realvec& position,
				  const Realvec& displacement )
  { return false; };

  int get_id(){ return face_id; };
  Realvec get_normal_vector(){ return normal; };
  Realvec get_represent_vector(){ return represent; };
};

#endif /*FACE_BASE_HPP*/
