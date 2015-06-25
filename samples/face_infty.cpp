#include "face_infty.hpp"

Realvec face_infty::move(Realvec& position, const Realvec& displacement,
		      boost::shared_ptr<face_base>& p)
{
  Realvec temp(3);
  temp = position + displacement;
  return temp;
//   boost::shared_ptr<face_base> sptr(this); // double free or corruption (fasttop)
//   p = sptr; // but successfully called twice
}

