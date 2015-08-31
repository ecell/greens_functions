#ifndef FACE_AllSYM
#define FACE_AllSYM
#include <utility>
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include "FaceBase.hpp"
#include "Polygon.hpp"
#include "Defs.hpp"

// const int RENEWLOOP_UPPER_LIMIT(100);

class FaceAllSym : public FaceBase
{
  std::vector<Realvec> vertexs;
  std::vector<Realvec> edges;

 //parametric
  Realvec para_origin;
  Realvec para_a;
  Realvec para_b;

public:
  FaceAllSym(const int& id, const Realvec& vtx0, const Realvec& vtx1, const Realvec& vtx2)
  : FaceBase( id, cross_product(vtx1-vtx0, vtx2-vtx0), vtx1-vtx0 ), vertexs(3), edges(3)
  {
    vertexs.at(0) = vtx0;
    vertexs.at(1) = vtx1;
    vertexs.at(2) = vtx2;

    edges.at(0) = vtx1-vtx0;
    edges.at(1) = vtx2-vtx1;
    edges.at(2) = vtx0-vtx2;

    para_origin = vertexs.at(0);
    para_a = edges.at(0);
    para_b = edges.at(2) * (-1e0);
  }
 
  FaceAllSym(const int& id, const Realvec& vtx0, const Realvec& vtx1, const Realvec& vtx2, const Realvec& norm)
  : FaceBase( id, norm, vtx1-vtx0 ), vertexs(3), edges(3)
  {
    bool normal_vec_is_oritented_orthogonally_to_the_edges(
		   dot_product(norm, vtx1-vtx0) == 0 &&
		   dot_product(norm, vtx2-vtx1) == 0 );
    THROW_UNLESS(std::invalid_argument, normal_vec_is_oritented_orthogonally_to_the_edges );

    vertexs.at(0) = vtx0;
    vertexs.at(1) = vtx1;
    vertexs.at(2) = vtx2;

    edges.at(0) = vtx1 - vtx0;
    edges.at(1) = vtx2 - vtx1;
    edges.at(2) = vtx0 - vtx2;

    para_origin = vertexs.at(0);
    para_a = edges.at(0);
    para_b = edges.at(2) * (-1e0);
  } 

  virtual std::pair<Realvec, FaceBase_sptr>
  apply_displacement( const Realvec& position, const Realvec& displacement,
		      const FaceBase_sptr& ptr);

  virtual bool in_face(const std::pair<Real, Real>& parameters, const Real tol = 1e-12);

  virtual int through_edge(const std::pair<Real, Real>& position,
			   const std::pair<Real, Real>& newposition, const Real tol = 1e-12);

  virtual Real cross_ratio(const std::pair<Real, Real>& position,
			   const std::pair<Real, Real>& displacement, const int& edge_num );

  virtual bool on_edge(const std::pair<Real, Real>& position, int& edge_num,
		       const Real tol = 1e-12);
  virtual bool on_edge(const std::pair<Real, Real>& position, const Real tol = 1e-12);

  virtual std::pair<Real, Real> to_parametric( const Realvec& pos );
 
  virtual Realvec to_absolute( const std::pair<Real, Real>& parameters );

  virtual Realvec get_para_origin(){return para_origin;}

  virtual Realvec get_para_a(){return para_a;}

  virtual Realvec get_para_b(){return para_b;}

  virtual void print_class_name();

private:

  std::pair<Real, Real> reverse( const std::pair<Real, Real>& dis, const int edge_num );

  virtual Realvec get_vertex();

  virtual Realvec get_another_vertex(const Realvec& edge);
};

std::pair<Realvec, FaceBase_sptr>
FaceAllSym::apply_displacement(const Realvec& position, const Realvec& displacement,
				const FaceBase_sptr& ptr )
{
// debug
//   std::cout << position[0] << " " << position[1] << " " << position[2] << " ";
//   std::cout << displacement[0] <<" "<< displacement[1] << " " << displacement[2] << std::endl;
 
  if( fabs( dot_product( position - para_origin, normal ) ) > 1e-12 )
    throw std::invalid_argument("apply_displacement: position is not on the plane");

  std::pair<Real, Real> pos_para( to_parametric( position - para_origin ) );
  std::pair<Real, Real> dis_para( to_parametric( displacement ) );
  const FaceBase_sptr face_ptr(ptr); //never change

//never go neighbor face
  for(int num_renew(0); num_renew < RENEWLOOP_UPPER_LIMIT; ++num_renew )
  {
    if(!in_face(pos_para))
      throw std::invalid_argument("apply_displacement: position is not in the face");

    std::pair<Real, Real> newpos( sum(pos_para, dis_para) );

    if( in_face( newpos ) )
    {
      Realvec renewed_pos( face_ptr->get_para_origin() + face_ptr->to_absolute( newpos ) );
      std::pair<Realvec, FaceBase_sptr> retpos( renewed_pos, face_ptr );
      return retpos;
    }

    const int edge_id( face_ptr->through_edge(pos_para, newpos) );

    const Real ratio( face_ptr->cross_ratio( pos_para, dis_para, edge_id) );

    std::pair<Real, Real> temppos( sum( pos_para, multiple(ratio, dis_para) ) );
    std::pair<Real, Real> tempdis( multiple( (1e0 - ratio), dis_para ) );
    assert( on_edge(temppos) );
    tempdis = reverse(tempdis, edge_id);

    pos_para = temppos;
    dis_para = tempdis;
  }

  throw std::invalid_argument("face ptr renewed over 100 times");
}

bool
FaceAllSym::in_face(const std::pair<Real, Real>& parameters, const Real tol )
{
  Real alpha( parameters.first );
  Real beta( parameters.second );

  bool just_on_vertex( (alpha== 0e0 || alpha == 1e0) &&
		       (beta == 0e0 ||  beta == 1e0) );
  if(just_on_vertex)
  {
    std::cout << "alpha: " << std::setprecision(16) << alpha << ", beta: " << beta << std::endl;
    throw std::invalid_argument("in_face : particle is just on a vertex");
  }

  bool alpha_in_range( -tol <= alpha && alpha <= 1e0 + tol );
  bool beta_in_range( -tol <= beta  && beta  <= 1e0 + tol );
  bool sum_in_range( -tol < alpha+beta && alpha+beta <= 1e0 + tol );

  return ( (alpha_in_range && beta_in_range) && sum_in_range);
}

std::pair<Real, Real> FaceAllSym::reverse( const std::pair<Real, Real>& dis, const int edge_num )
{
  switch(edge_num)
  {
    case 0:
      {
        std::pair<Real, Real> reversed(dis.first, -dis.second);
	return reversed;
      }
      break;

    case 1:
      {
	std::pair<Real, Real> reversed(-dis.second, -dis.first);
	return reversed;
      }
      break;

    case 2:
      {
	std::pair<Real, Real> reversed(-dis.first, dis.second);
	return reversed;
      }
      break;

    default:
      throw std::invalid_argument("reverse: invalid edge_num");
      break;
  }
  throw std::invalid_argument("reverse: switch passed");
}

int FaceAllSym::through_edge( const std::pair<Real, Real>& position,
	  		       const std::pair<Real, Real>& newposition, const Real tol)
{
  Real pos_alpha(position.first);
  Real pos_beta(position.second);
  Real newpos_alpha(newposition.first);
  Real newpos_beta(newposition.second);

  //case::first position is on the edge
  int edge_num(-1);
  if( on_edge( position, edge_num ) )
  {
    switch(edge_num)
    {
      //edges.at(0) <=> beta==0
      case 0:
	if(newpos_beta < 0e0)
	  return 0;
	break;
      case 1:
	if(newpos_alpha + newpos_beta > 1e0)
	  return 1;
	break;
      case 2:
	if(newpos_alpha < 0e0)
	  return 2;
	break;
      default:
	throw std::invalid_argument("on_edge returns invalid edge number");
	break;
    }
  }
/* case
//    \ 5/
//     \/
//   2 /\ 1
//____/__\___
// 3 / 0  \ 4
//  /      \
*/
  //through one edge: case {0,1,2}
  if( 0e0 <= newpos_alpha && newpos_beta < 0e0 && newpos_alpha + newpos_beta <= 1e0 )
    return 0;
  if( 0e0 < newpos_alpha && 0e0 < newpos_beta && 1e0 < newpos_alpha + newpos_beta )
    return 1;
  if( 0e0 <= newpos_beta && newpos_alpha < 0e0 && newpos_alpha + newpos_beta <= 1e0 )
    return 2;

  Real dis_alpha(newpos_alpha - pos_alpha);
  Real dis_beta(newpos_beta - pos_beta);
  std::pair<Real, Real> displacement( dis_alpha, dis_beta );
  //through two edge: case {3,4,5}
  if( newpos_alpha < 0e0 && newpos_beta < 0e0 )
  {
  // case 3
  //throwgh edges.at(0) and at(2)
    Real ratio0( cross_ratio(position, displacement, 0) );
    Real ratio2( cross_ratio(position, displacement, 2) );
    if( ratio0 < ratio2 )
    {
      return 0;
    }else if(ratio0 > ratio2)
    {
      return 2;
    }else{
      throw std::invalid_argument( "particle goes through vertex!" );
    }
  }
  else if( newpos_alpha > 0e0 && newpos_beta < 0e0 && newpos_alpha + newpos_beta > 1e0)
  {
  // case 4
  //throwgh edges.at(0) and at(1)
    Real ratio0( cross_ratio(position, displacement, 0) );
    Real ratio1( cross_ratio(position, displacement, 1) );
    if( ratio0 < ratio1 )
    {
      return 0;
    }else if(ratio0 > ratio1)
    {
      return 1;
    }else{
      throw std::invalid_argument( "particle goes through vertex!" );
    }
  }
  else if( newpos_alpha < 0e0 && newpos_beta > 0e0 && newpos_alpha + newpos_beta > 1e0)
  {
  // case 5
  //throwgh edges.at(1) and at(2)
    Real ratio1( cross_ratio(position, displacement, 1) );
    Real ratio2( cross_ratio(position, displacement, 2) );
    if( ratio1 < ratio2 )
    {
      return 1;
    }else if(ratio1 > ratio2)
    {
      return 2;
    }else{
      throw std::invalid_argument( "particle goes through vertex!" );
    }
  }
  throw std::invalid_argument("through_edge::particle does not go out of this face");
  return -1;
}

Real
FaceAllSym::cross_ratio( const std::pair<Real, Real>& position,
			  const std::pair<Real, Real>& displacement,
			  const int& edge_num )
{
  if(!in_face(position))
    throw std::invalid_argument("cross: position is not in the face");
  Real pos_alpha( position.first );
  Real pos_beta( position.second );
  Real dis_alpha( displacement.first );
  Real dis_beta( displacement.second );

  if( dis_alpha == 0e0 && dis_beta == 0e0 )
    throw std::invalid_argument("ratio: displacement is zero vector");

  Real ratio(-1e0);
  switch(edge_num)
  {
  case 0:
    if( dis_beta == 0e0 )
      throw std::invalid_argument("cannot intersect this edge");
    ratio = - pos_beta / dis_beta;
    break;
  case 1:
    if( dis_alpha + dis_beta == 0e0 )
      throw std::invalid_argument("cannot intersect this edge");
    ratio = ( 1e0 - pos_alpha - pos_beta ) / ( dis_alpha + dis_beta );
    break;
  case 2:
    if( dis_alpha == 0e0 )
      throw std::invalid_argument("cannot intersect this edge");
    ratio = - pos_alpha / dis_alpha;
    break;
  default:
    throw std::invalid_argument("invalid edge_num");
  }

  return ratio;
}

/*     2
//    /\
// b /  \
// ^/____\
// 0 ->a  1
*/
bool
FaceAllSym::on_edge(const std::pair<Real, Real>& position,
		     int& edge_num, const Real tol)
{
  Real alpha(position.first);
  Real beta(position.second);

  if( alpha == 0e0 && (-tol <= beta && beta <= 1e0+tol) )
  {
    edge_num = 2;
    return true;
  }
  else if( beta == 0e0 && ( -tol <= alpha && alpha <= 1e0+tol ) )
  {
    edge_num = 0;
    return true;
  }
  else if( fabs( alpha + beta - 1e0 ) < tol && (-tol < alpha && -tol < beta) )
  {
    edge_num = 1;
    return true;
  }else
  {
    edge_num = -1;
    return false;
  }
}

bool
FaceAllSym::on_edge(const std::pair<Real, Real>& position,
		     const Real tol)
{
  Real alpha(position.first);
  Real beta(position.second);

  if( (-tol <= alpha && alpha <= tol) && (-tol <= beta && beta <= 1e0+tol) )
  {
    return true;
  }
  else if( (-tol <= beta && beta <= tol) && ( -tol <= alpha && alpha <= 1e0+tol ) )
  {
    return true;
  }
  else if( fabs( alpha + beta - 1e0 ) < tol && (-tol < alpha && -tol < beta) )
  {
    return true;
  }else
  {
    return false;
  }
}

std::pair<Real, Real> FaceAllSym::to_parametric( const Realvec& pos )
{
  Real alpha, beta;
  Realvec parametric_pos( pos );

  if( (para_a[0]*para_b[1] - para_b[0]*para_a[1]) != 0e0 )
  {
  // matrix( a0 b0 )
  //       ( a1 b1 )
    alpha = ( para_b[1] * parametric_pos[0] - para_b[0] * parametric_pos[1] )
	    / (para_a[0]*para_b[1] - para_b[0]*para_a[1]);

     beta = ( para_a[0] * parametric_pos[1] - para_a[1] * parametric_pos[0] )
	    / (para_a[0]*para_b[1] - para_b[0]*para_a[1]);
 
    //confirm (alpha beta) * (a2 b2) = (pos2)
    if( fabs(alpha * para_a[2] + beta * para_b[2] - parametric_pos[2]) > 1e-12 ) 
      throw std::invalid_argument( "function to_parametric: invalid solution" );
    
    std::pair<Real, Real> parameters(alpha, beta);
    return parameters;

  }else if( (para_a[1]*para_b[2] - para_b[1]*para_a[2]) != 0e0 )
  {
  // matrix( a1 b1 )
  //       ( a2 b2 )
    alpha = ( para_b[2] * parametric_pos[1] - para_b[1] * parametric_pos[2] )
	    / (para_a[1]*para_b[2] - para_b[1]*para_a[2]);

     beta = ( para_a[1] * parametric_pos[2] - para_a[2] * parametric_pos[1] )
	    / (para_a[1]*para_b[2] - para_b[1]*para_a[2]);
 
    //confirm (alpha beta) * (a0 b0) = (pos0)
    if( fabs(alpha * para_a[0] + beta * para_b[0] - parametric_pos[0]) > 1e-12 )
      throw std::invalid_argument( "function to_parametric: invalid solution" );
      
    std::pair<Real, Real> parameters(alpha, beta);
    return parameters;

  }else if( (para_a[2]*para_b[0] - para_b[2]*para_a[0]) != 0e0 )
  {
  // matrix( a2 b2 )
  //       ( a0 b0 )
    alpha = ( para_b[0] * parametric_pos[2] - para_b[2] * parametric_pos[0] )
	    / (para_a[2]*para_b[0] - para_b[2]*para_a[0]);

     beta = ( para_a[2] * parametric_pos[0] - para_a[0] * parametric_pos[2] )
	    / (para_a[2]*para_b[0] - para_b[2]*para_a[0]);
 
    //confirm (alpha beta) * (a1 b1) = (pos1)
    if( fabs(alpha * para_a[1] + beta * para_b[1] - parametric_pos[1]) > 1e-12 ) 
      throw std::invalid_argument( "function to_parametric: invalid solution" ); 
      
    std::pair<Real, Real> parameters(alpha, beta);
    return parameters;

  }

  throw std::invalid_argument( "function to_parametric: could not parametrise input position" );
}

Realvec FaceAllSym::to_absolute( const std::pair<Real, Real>& parameters )
{
  Realvec absolute_pos( para_a * parameters.first + para_b * parameters.second );
  return absolute_pos;
}

Realvec FaceAllSym::get_vertex()
{
  return para_origin; 
};

Realvec FaceAllSym::get_another_vertex(const Realvec& edge)
{
  for(int i(0); i<3; ++i)
  {
    if( edges.at(i)[0] == edge[0] &&
	edges.at(i)[1] == edge[1] &&
	edges.at(i)[2] == edge[2])
    {
      int vertex_id( (i+2) % 3 );
      return vertexs.at( vertex_id );
    }
  }

  bool get_another_vertex_have_invalid_edge_input(false);
  THROW_UNLESS( std::invalid_argument, get_another_vertex_have_invalid_edge_input== true );

  Realvec zero;
  return zero;
}

void FaceAllSym::print_class_name()
{
  std::cout << "class: FaceAllSym" << std::endl;
}

#endif /*FACE_AllSYM*/
