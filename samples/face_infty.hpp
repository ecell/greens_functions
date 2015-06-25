#ifndef FACE_INFTY_HPP
#define FACE_INFTY_HPP
#include "face_base.hpp"

using namespace greens_functions;

class face_infty : public face_base
{
public:
  face_infty( int temp ) : face_base( temp )
  {
    ;
  };
  ~face_infty(){
    std::cout << "face_inf destructer called" << std::endl;
  }
  virtual Realvec move(Realvec& position, const Realvec& displacement,
		    boost::shared_ptr<face_base>& p);
};

#endif /*FACE_INFTY_HPP*/
