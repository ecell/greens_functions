#ifndef FACE_BASE_HPP
#define FACE_BASE_HPP
#include <iostream>
#include <utility>
#include <vector>
#include <cmath>
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

  //parametric
  Realvec para_origin;
  Realvec para_a;
  Realvec para_b;
  // neighbors
  //para_a = para_a_neighbor.first * neighbor->para_a + para_a_neighbor.second * neighbor->para_b;
  //para_b = para_b_neighbor.first * neighbor->para_a + para_b_neighbor.second * neighbor->para_b;
  //p(a, b) -> p'( a*para_a_neighbor.first + b*para_b_neighbor.first,
  //		   a*para_a_neighbor.second + b*para_b_neighbor.second)
  std::vector< std::pair<Real, Real> > para_a_neighbor;
  std::vector< std::pair<Real, Real> > para_b_neighbor;
  // para_origin = ori_vec_neighbor.first * neighbor->para_a + 
  // 		   ori_vec_neighbor.second * neighbor->para_b;
  std::vector< std::pair<Real, Real> > ori_vec_neighbor;


public:
  FaceBase( const int& id, const Realvec& norm, const Realvec& rep )
  : face_id(id), normal( norm/length(norm) ), represent( rep/length(rep) ),
    para_a_neighbor(3), para_b_neighbor(3), ori_vec_neighbor(3)
  {
    THROW_UNLESS( std::invalid_argument, dot_product(norm, rep) == 0 );
    THROW_UNLESS( std::invalid_argument, length( norm ) != 0 );
    THROW_UNLESS( std::invalid_argument, length( rep ) != 0 );
  };

  virtual Realvec
  renew_position(const Realvec& position, const Realvec& displacement,
		 boost::shared_ptr<FaceBase>& p )
  {
    print_class_name();
    throw std::invalid_argument("renew_position: this is not supported");
  }

  virtual std::pair<Realvec, boost::shared_ptr<FaceBase> >
  apply_displacement(const Realvec& position, const Realvec& displacement,
		 const boost::shared_ptr<FaceBase>& p )
  {
    print_class_name();
    throw std::invalid_argument("apply_displacement: this is not supported");
  }

  //return whitch edge the displacement(newpos-pos) goes through.
  virtual int through_edge(const std::pair<Real, Real>& position,
			   const std::pair<Real, Real>& newposition, const Real tol=1e-12)
  {
    print_class_name();
    throw std::invalid_argument("through_edge: this is not supported");
  }

  //return the ratio(0~1) of a displacement segment not cross the edge.
  virtual Real cross_ratio(const std::pair<Real, Real>& position,
		   const std::pair<Real, Real>& displacement, const int& edge_num )
  {
    print_class_name();
    throw std::invalid_argument("cross_ratio: this is not supported");
  }

  //return whether the place is on the edge and rewrite edge_num to the edge.
  //  if false, edge_num = -1.
  virtual bool on_edge(const std::pair<Real, Real>& position, int& edge_num, const Real tol = 1e-12)
  {
    print_class_name();
    throw std::invalid_argument("on_edge(pair, int, Real): this is not supported");
  }

  virtual bool on_edge(const std::pair<Real, Real>& position, const Real tol = 1e-12)
  {
    print_class_name();
    throw std::invalid_argument("on_edge(pair, Real): this is not supported");
  }


  // return pair of parameters of parametric expression of vector pos.
  // pair.first = alpha, pair.second = beta
  virtual std::pair<Real, Real> to_parametric( const Realvec& pos )
  {
    print_class_name();
    throw std::invalid_argument("to_parametric: this is not supported");
  }
  // return absolute expression translated from parametric expression.
  virtual Realvec to_absolute( const std::pair<Real, Real>& parameters )
  {
    print_class_name();
    throw std::invalid_argument("to_absolute: this is not supported");
  }

  virtual bool in_the_face( const Realvec& position, const Realvec& displacement )
  {
    print_class_name();
    throw std::invalid_argument("in_the_face: this is not supported");
  }

  virtual bool in_face(const std::pair<Real, Real> parameters, const Real tol = 1e-12)
  {
    Real alpha( parameters.first );
    Real beta( parameters.second );
    
    bool alpha_in_range( -tol <= alpha && alpha <= 1e0 + tol );
    bool beta_in_range( -tol <= beta  && beta  <= 1e0 + tol );
    bool sum_in_range( -tol < alpha+beta && alpha+beta <= 1e0 + tol );

    return ( (alpha_in_range && beta_in_range) && sum_in_range);
  }

  virtual Realvec get_vertex()
  {
    print_class_name();
    throw std::invalid_argument("get_vertex: this is not supported");
  }

  virtual void set_poly_ptr( boost::shared_ptr<Polygon>& p_sptr )
  {
    print_class_name();
    throw std::invalid_argument("set_poly_ptr: this is not supported");
  };

  virtual void set_near_vertexs()
  {
    print_class_name();
    throw std::invalid_argument("set_near_vertex: this is not supported");
  }

  virtual Realvec get_another_vertex( const Realvec& edge )
  {
    print_class_name();
    throw std::invalid_argument("get_another_vertex: this is not supported");
  };

  virtual Real get_max_a(const Realvec& position, bool& vertex_involve_flag)
  {
    print_class_name();
    throw std::invalid_argument("get_max_a: this is not supported");
  };

  virtual Real get_minimum_height(const Realvec& neighbors_edge)
  {
    print_class_name();
    throw std::invalid_argument("get_minimum_height: this is not supported");
  };

  virtual Real get_left_angle( const Realvec& neighbors_edge )
  {
    print_class_name();
    throw std::invalid_argument("get_left_angle: this is not supported");
  }

  virtual Real get_right_angle( const Realvec& neighbors_edge )
  {
    print_class_name();
    throw std::invalid_argument("get_right_angle: this is not supported");
  }

  virtual void print_class_name()
  {
    std::cout << "class: FaceBase" << std::endl;
    return;
  }

  int get_id(){ return face_id; };
  Realvec get_normal_vector(){ return normal; };
  Realvec get_represent_vector(){ return represent; };
  Realvec get_para_origin(){return para_origin;}
  Realvec get_para_a(){return para_a;}
  Realvec get_para_b(){return para_b;}
  std::pair<Real, Real> get_para_a_neighbor_at(const int i)
  {
    return para_a_neighbor.at(i);
  }
  std::pair<Real, Real> get_para_b_neighbor_at(const int i)
  {
    return para_b_neighbor.at(i);
  }
  std::pair<Real, Real> get_ori_vec_neighbor_at(const int i)
  {
    return ori_vec_neighbor.at(i);
  }

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

std::pair<Real, Real> sum( const std::pair<Real, Real>& lhs, const std::pair<Real, Real>& rhs )
{
  std::pair<Real, Real> retpair( lhs.first + rhs.first, lhs.second + rhs.second );
  return retpair;
}

std::pair<Real, Real> subtract(const std::pair<Real, Real>& lhs, const std::pair<Real, Real>& rhs)
{
  std::pair<Real, Real> retpair( lhs.first - rhs.first, lhs.second - rhs.second );
  return retpair;
}

std::pair<Real, Real> multiple( const Real lhs, const std::pair<Real, Real>& rhs )
{
  std::pair<Real, Real> retpair( lhs * rhs.first, lhs * rhs.second );
  return retpair;
}

std::pair<Real, Real> multiple( const std::pair<Real, Real>& lhs, const Real rhs )
{
  std::pair<Real, Real> retpair( rhs * lhs.first, rhs * lhs.second );
  return retpair;
}


typedef boost::shared_ptr<FaceBase> FaceBase_sptr;
#endif /*FACE_BASE_HPP*/
