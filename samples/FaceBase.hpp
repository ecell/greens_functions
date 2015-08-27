#ifndef FACE_BASE_HPP
#define FACE_BASE_HPP
#include <iostream>
#include <utility>
#include <vector>
#include <cmath>
#include <cassert>
#include <boost/shared_ptr.hpp>
#include <boost/array.hpp>
#include "Defs.hpp"
#include "Vector3.hpp"

class Polygon;
typedef Vector3<Real> Realvec;

class FaceBase
{
protected:
  int face_id;
  Realvec normal;
  Realvec represent;

public:
  FaceBase( const int id, const Realvec& norm, const Realvec& rep )
  : face_id(id), normal( norm/length(norm) ), represent( rep/length(rep) )
  {
    assert( dot_product(norm, rep) == 0e0 );
    assert( length( norm ) != 0e0 );
    assert( length( rep ) != 0e0 );
  }

  virtual void set_poly_ptr( boost::shared_ptr<Polygon>& p_sptr );

  virtual std::pair<Realvec, boost::shared_ptr<FaceBase> >
  apply_displacement(const Realvec& position, const Realvec& displacement,
		     const boost::shared_ptr<FaceBase>& p );

  virtual bool in_face(const std::pair<Real, Real>& parameters, const Real tol = 1e-12);

  virtual bool on_vertex(const std::pair<Real, Real>& parameters, const Real tol = 1e-12);

  //return whitch edge the displacement(newpos-pos) goes through.
  virtual int  through_edge(const std::pair<Real, Real>& position,
                            const std::pair<Real, Real>& newposition, const Real tol=1e-12);

  //return the ratio(0~1) of a displacement segment not cross the edge.
  virtual Real cross_ratio(const std::pair<Real, Real>& position,
			   const std::pair<Real, Real>& displacement, const int& edge_num );
  
  //return whether the place is on the edge and rewrite edge_num to the edge.
  //  if false, edge_num = -1.
  virtual bool on_edge(const std::pair<Real, Real>& position, int& edge_num,
		       const Real tol = 1e-12);

  virtual bool on_edge(const std::pair<Real, Real>& position, const Real tol = 1e-12);

  // return pair of parameters of parametric expression of vector pos.
  // pair.first = alpha, pair.second = beta
  virtual std::pair<Real, Real> to_parametric( const Realvec& pos );
  
  // return absolute expression translated from parametric expression.
  virtual Realvec to_absolute( const std::pair<Real, Real>& parameters );
 
  virtual void set_near_vertexs();

  virtual Realvec get_another_vertex( const Realvec& edge );

  virtual Real get_max_a(const Realvec& position, bool& vertex_involve_flag);

  virtual Real get_minimum_height(const Realvec& neighbors_edge);

  virtual Real get_right_angle( const Realvec& neighbors_edge );
  virtual Real get_left_angle( const Realvec& neighbors_edge );

  virtual bool is_gate_at(int edge_id);

  int get_id(){ return face_id; }

  Realvec get_normal_vector(){ return normal; }

  Realvec get_represent_vector(){ return represent; }

  virtual Realvec get_para_origin();

  virtual Realvec get_para_a();

  virtual Realvec get_para_b();

  virtual std::pair<Real, Real> get_para_a_neighbor_at(const int i);

  virtual std::pair<Real, Real> get_para_b_neighbor_at(const int i);

  virtual std::pair<Real, Real> get_ori_vec_neighbor_at(const int i);

  virtual void print_class_name();

  Real smaller_angle(const Realvec& v1, const Realvec& v2);
};


std::pair<Realvec, boost::shared_ptr<FaceBase> >
FaceBase::apply_displacement(const Realvec& position, const Realvec& displacement,
  const boost::shared_ptr<FaceBase>& p )
{
  print_class_name();
  throw std::invalid_argument("apply_displacement: this is not supported");
}

int
FaceBase::through_edge(const std::pair<Real, Real>& position,
  const std::pair<Real, Real>& newposition, const Real tol )
{
  print_class_name();
  throw std::invalid_argument("through_edge: this is not supported");
}

Real
FaceBase::cross_ratio(const std::pair<Real, Real>& position,
  const std::pair<Real, Real>& displacement, const int& edge_num )
{
  print_class_name();
  throw std::invalid_argument("cross_ratio: this is not supported");
}

bool
FaceBase::on_edge(const std::pair<Real, Real>& position, int& edge_num, const Real tol)
{
  print_class_name();
  throw std::invalid_argument("on_edge(pair, int, Real): this is not supported");
}

bool
FaceBase::on_edge(const std::pair<Real, Real>& position, const Real tol)
{
  print_class_name();
  throw std::invalid_argument("on_edge(pair, Real): this is not supported");
}

// return pair of parameters of parametric expression of vector pos.
// pair.first = alpha, pair.second = beta
std::pair<Real, Real>
FaceBase::to_parametric( const Realvec& pos )
{
  print_class_name();
  throw std::invalid_argument("to_parametric: this is not supported");
}

// return absolute expression translated from parametric expression.
Realvec
FaceBase::to_absolute( const std::pair<Real, Real>& parameters )
{
  print_class_name();
  throw std::invalid_argument("to_absolute: this is not supported");
}

bool
FaceBase::in_face(const std::pair<Real, Real>& parameters, const Real tol)
{
// if(on_vertex(parameters))
//   throw std::invalid_argument("in_face::particle is on vertex");

  Real alpha( parameters.first );
  Real beta( parameters.second );
  
  bool alpha_in_range( -tol <= alpha && alpha <= 1e0 + tol );
  bool beta_in_range( -tol <= beta  && beta  <= 1e0 + tol );
  bool sum_in_range( -tol < alpha+beta && alpha+beta <= 1e0 + tol );

  return ( (alpha_in_range && beta_in_range) && sum_in_range);
}

bool
FaceBase::on_vertex(const std::pair<Real, Real>& parameters, const Real tol)
{
  Real alpha( parameters.first );
  Real beta( parameters.second );
  
  bool vtx0( (-tol < alpha && alpha < tol) && (-tol < beta && beta < tol) );
  bool vtx1( (-tol < alpha && alpha < tol) && (1e0-tol < beta && beta < 1e0+tol) );
  bool vtx2( (1e0-tol < alpha && alpha < 1e0+tol) && (-tol < beta && beta < tol) );

  return ( vtx0 && vtx1 && vtx2 );
}

void 
FaceBase::set_poly_ptr( boost::shared_ptr<Polygon>& p_sptr )
{
  print_class_name();
  throw std::invalid_argument("set_poly_ptr: this is not supported");
}

void 
FaceBase::set_near_vertexs()
{
  print_class_name();
  throw std::invalid_argument("set_near_vertex: this is not supported");
}

Realvec 
FaceBase::get_another_vertex( const Realvec& edge )
{
  print_class_name();
  throw std::invalid_argument("get_another_vertex: this is not supported");
}

Real FaceBase::get_max_a(const Realvec& position, bool& vertex_involve_flag)
{
  print_class_name();
  throw std::invalid_argument("get_max_a: this is not supported");
}

Real FaceBase::get_minimum_height(const Realvec& neighbors_edge)
{
  print_class_name();
  throw std::invalid_argument("get_minimum_height: this is not supported");
}

Real FaceBase::FaceBase::get_left_angle( const Realvec& neighbors_edge )
{
  print_class_name();
  throw std::invalid_argument("get_left_angle: this is not supported");
}

Real FaceBase::get_right_angle( const Realvec& neighbors_edge )
{
  print_class_name();
  throw std::invalid_argument("get_right_angle: this is not supported");
}

void FaceBase::print_class_name()
{
  std::cout << "class: FaceBase" << std::endl;
  return;
}

bool FaceBase::is_gate_at(int edge_id)
{
  print_class_name();
  throw std::invalid_argument("is_gate_at: this is not supported");
}

std::pair<Real, Real> FaceBase::get_para_a_neighbor_at(const int i)
{
  print_class_name();
  throw std::invalid_argument("get_para_a_neighbor_at: this is not supported");
}

std::pair<Real, Real> FaceBase::get_para_b_neighbor_at(const int i)
{
  print_class_name();
  throw std::invalid_argument("get_para_b_neighbor_at: this is not supported");
}

std::pair<Real, Real> FaceBase::get_ori_vec_neighbor_at(const int i)
{
  print_class_name();
  throw std::invalid_argument("get_ori_vec_neighbor_at: this is not supported");
}

Realvec FaceBase::get_para_origin()
{
  print_class_name();
  throw std::invalid_argument("get_para_origin: this is not supported");
}

Realvec FaceBase::get_para_a()
{
  print_class_name();
  throw std::invalid_argument("get_para_a: this is not supported");
}

Realvec FaceBase::get_para_b()
{
  print_class_name();
  throw std::invalid_argument("get_para_b: this is not supported");
}


Real FaceBase::smaller_angle(const Realvec& v1, const Realvec& v2)
{
  Real len1( length(v1) );
  Real len2( length(v2) );
  Real inner( dot_product( v1, v2 ) );

  Real angle( acos( inner / len1 / len2 ) );
  return angle;
}


//**************************** std::pair ****************************************
std::pair<Real, Real>
sum( const std::pair<Real, Real>& lhs, const std::pair<Real, Real>& rhs )
{
  std::pair<Real, Real> retpair( lhs.first + rhs.first, lhs.second + rhs.second );
  return retpair;
}

std::pair<Real, Real>
subtract(const std::pair<Real, Real>& lhs, const std::pair<Real, Real>& rhs)
{
  std::pair<Real, Real> retpair( lhs.first - rhs.first, lhs.second - rhs.second );
  return retpair;
}

std::pair<Real, Real>
multiple( const Real lhs, const std::pair<Real, Real>& rhs )
{
  std::pair<Real, Real> retpair( lhs * rhs.first, lhs * rhs.second );
  return retpair;
}

std::pair<Real, Real>
multiple( const std::pair<Real, Real>& lhs, const Real rhs )
{
  std::pair<Real, Real> retpair( rhs * lhs.first, rhs * lhs.second );
  return retpair;
}


typedef boost::shared_ptr<FaceBase> FaceBase_sptr;
#endif /*FACE_BASE_HPP*/
