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
public:
  face_base( int id ){ face_id = id; };
  virtual Realvec move(Realvec& position, const Realvec& displacement,
		    boost::shared_ptr<face_base>& p)
  {
    Realvec zero(3);
      zero[0] = 0e0;
      zero[1] = 0e0;
      zero[2] = 0e0;

    std::cout << "face_base->move" << std::endl;
    return zero;
  };
  int get_id(){ return face_id; };
};

#endif /*FACE_BASE_HPP*/
