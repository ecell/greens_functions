#ifndef FACE_INF_HPP
#define FACE_INF_HPP
#include "FaceBase.hpp"

class FaceInf : public FaceBase
{
private:
  Realvec para_origin;
  
public:
  FaceInf( const int& id, const Realvec& norm, const Realvec& rep)
  : FaceBase( id, norm, rep )
  {
    Realvec zero(0e0, 0e0, 0e0);
    para_origin = zero;
  }

  FaceInf( const int& id, const Realvec& norm, const Realvec& rep, const Realvec& origin)
  : FaceBase( id, norm, rep ), para_origin(origin)
  {
  }

  ~FaceInf(){ };

  virtual std::pair<Realvec, boost::shared_ptr<FaceBase> >
  apply_displacement(const Realvec& position, const Realvec& displacement,
		     const boost::shared_ptr<FaceBase>& p );


  virtual bool in_face( const Realvec& position, const Realvec& displacement, const Real tol = 1e-12 );
//   virtual bool in_face( const std::pair<Real, Real>& parameters, const Real tol = 1e-12 );

  virtual Realvec get_para_origin();
  virtual void print_class_name();
};

std::pair<Realvec, boost::shared_ptr<FaceBase> >
FaceInf::apply_displacement(const Realvec& position, const Realvec& displacement, const FaceBase_sptr& p )
{
  if( !in_face( position, displacement ) )
    throw std::invalid_argument("newpos is not on this face");

  std::pair<Realvec, FaceBase_sptr> retval(position + displacement, p);
  return retval;
}

bool FaceInf::in_face( const Realvec& position, const Realvec& displacement, const Real tol )
{
  Realvec p( position + displacement - para_origin );
  return ( fabs(dot_product(p, normal) ) < tol );
}

Realvec FaceInf::get_para_origin()
{
  return para_origin;
};

void FaceInf::print_class_name()
{
  std::cout << "class: FaceInf" << std::endl;
}

#endif /*FACE_INF_HPP*/
