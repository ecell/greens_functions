#ifndef FACE_INFTY_HPP
#define FACE_INFTY_HPP
#include "face_base.hpp"

using namespace greens_functions;

class face_infty : public face_base
{
public:
  face_infty( const int& id, const Realvec& norm, const Realvec& rep )
  : face_base( id, norm, rep )
  {
    ;
  };

  ~face_infty(){
    std::cout << "face_inf destructer called" << std::endl;
  };

  virtual Realvec move(Realvec& position,
		       const Realvec& displacement,
		       boost::shared_ptr<face_base>& p)
  {
    bool in_the_infty_plane;
    in_the_infty_plane = still_in_the_face(position, displacement);
    THROW_UNLESS(std::invalid_argument, in_the_infty_plane);

    Realvec temp(position + displacement);
    return temp;
    
  };

  virtual bool still_in_the_face( const Realvec& position,
				  const Realvec& displacement )
  {
    return true;
  };

};

#endif /*FACE_INFTY_HPP*/
