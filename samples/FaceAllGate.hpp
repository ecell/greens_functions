#ifndef FACE_ALL_GATE
#define FACE_ALL_GATE
#include <vector>
#include <utility>
#include <iostream>
#include <cmath>
#include <cassert>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include "FaceBase.hpp"
#include "Polygon.hpp"
#include "rotation.hpp"

const int RENEWLOOP_UPPER_LIMIT(100);

class FaceAllGate : public FaceBase
{
private:
  std::vector<Realvec> vertexs;
  std::vector<Realvec> edges;
  std::vector<Real> angles;
  
  boost::weak_ptr<Polygon> poly_ptr;
  std::vector<Real> near_vert_height;

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
  FaceAllGate(const int id, const Realvec& vtx0, const Realvec& vtx1, const Realvec& vtx2)
  : FaceBase( id, cross_product(vtx1-vtx0, vtx2-vtx0), vtx1-vtx0 ),
    vertexs(3), edges(3), angles(3), para_a_neighbor(3), para_b_neighbor(3), ori_vec_neighbor(3)
  {
    vertexs.at(0) = vtx0;
    vertexs.at(1) = vtx1;
    vertexs.at(2) = vtx2;

    //edges.at(i): vertexs.at(i) -> vertexs.at(i+1)
    edges.at(0) = vtx1 - vtx0;
    edges.at(1) = vtx2 - vtx1;
    edges.at(2) = vtx0 - vtx2;
  
    angles.at(0) = smaller_angle( edges.at(0), edges.at(2)*(-1e0) );
    angles.at(1) = smaller_angle( edges.at(1), edges.at(0)*(-1e0) );
    angles.at(2) = smaller_angle( edges.at(2), edges.at(1)*(-1e0) );

    // define parametric expression
    para_origin = vertexs.at(0);
    para_a = edges.at(0);
    para_b = edges.at(2) * (-1e0);
  }

  FaceAllGate(const int& id, const Realvec& vtx0, const Realvec& vtx1, const Realvec& vtx2,
	      const Realvec& norm)
  : FaceBase( id, norm, vtx1-vtx0 ), vertexs(3), edges(3), angles(3),
    para_a_neighbor(3), para_b_neighbor(3), ori_vec_neighbor(3)
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

    angles.at(0) = smaller_angle( edges.at(0), edges.at(2)*(-1e0) );
    angles.at(1) = smaller_angle( edges.at(1), edges.at(0)*(-1e0) );
    angles.at(2) = smaller_angle( edges.at(2), edges.at(1)*(-1e0) );

    // define parametric expression
    para_origin = vertexs.at(0);
    para_a = edges.at(0);
    para_b = edges.at(2) * (-1e0);
  }

  virtual void set_poly_ptr( boost::shared_ptr<Polygon>& p_sptr)
  {
    poly_ptr = p_sptr;
  }

  virtual std::pair<Realvec, FaceBase_sptr>
  apply_displacement( const Realvec& position, const Realvec& displacement,
		      const FaceBase_sptr& ptr);

  //return whether the position is in a face ( in the range or not ).
  virtual bool in_face(const std::pair<Real, Real>& parameters, const Real tol = 1e-12);

  //return whitch edge the displacement(newpos-pos) goes through.
  virtual int through_edge(const std::pair<Real, Real>& position,
			   const std::pair<Real, Real>& newposition, const Real tol = 1e-12);

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

//set near_vert_height.
  virtual void set_near_vertexs();

  //return the vertex such that input edge does not include.
  virtual Realvec get_another_vertex(const Realvec& edge);

  //find max shell size (for greens function)
  virtual Real get_max_a(const Realvec& position, bool& vertex_involve_flag);
  
  //call this func from neighbor face
  virtual Real get_minimum_height( const Realvec& neighbors_edge );

  //return own angle.
  virtual Real get_right_angle( const Realvec& neighbors_edge );
  virtual Real get_left_angle( const Realvec& neighbors_edge );

  virtual Realvec get_para_origin(){return para_origin;}

  virtual Realvec get_para_a(){return para_a;}

  virtual Realvec get_para_b(){return para_b;}

  virtual std::pair<Real, Real> get_para_a_neighbor_at(const int i)
  {return para_a_neighbor.at(i);}

  virtual std::pair<Real, Real> get_para_b_neighbor_at(const int i)
  {return para_b_neighbor.at(i);}

  virtual std::pair<Real, Real> get_ori_vec_neighbor_at(const int i)
  {return ori_vec_neighbor.at(i);}

  virtual void print_class_name();

private:
  std::pair<Real, Real> 
  translate_pos(std::pair<Real, Real> pos, int edge, FaceBase_sptr face_ptr);

  std::pair<Real, Real> 
  translate_dis(std::pair<Real, Real> dis, int edge, FaceBase_sptr face_ptr);

  // parametrize para_a & para_b using neighbor's apara_b
  void set_neighbors_edge();
  void set_neighbors_ori_vtx();

  //return vertex id that is not on the edge. 
  //  when neighboring face execute this function through ptr, 
  //  this func return the vertex that is not shared with the neighbor.
  int get_another_vertex_id(const Realvec& edge);

  //return neighbor's non-shared vertex vector
  Realvec get_another_vertex_right( const Realvec& neighbors_edge);
  Realvec get_another_vertex_left( const Realvec& neighbors_edge);

  //return own edge id
  int get_another_edge_id_right( const Realvec& neighbors_edge);
  int get_another_edge_id_left( const Realvec& neighbors_edge);

  // return neighbor's right angle from pointer pulled from polygon class.
  Real pull_neighbors_right_angle( const Realvec& neighbors_edge );
  Real pull_neighbors_left_angle( const Realvec& neighbors_edge );

};

std::pair<Realvec, FaceBase_sptr>
FaceAllGate::apply_displacement(const Realvec& position, const Realvec& displacement,
				const FaceBase_sptr& ptr )
{
// debug
//   std::cout << position[0] << " " << position[1] << " " << position[2] << " ";
//   std::cout << displacement[0] <<" "<< displacement[1] << " " << displacement[2] << std::endl;
 
  if( fabs( dot_product( position - para_origin, normal ) ) > 1e-12 )
    throw std::invalid_argument("renew_position: position is not on the plane");

  std::pair<Real, Real> pos_para( to_parametric( position - para_origin ) );
  std::pair<Real, Real> dis_para( to_parametric( displacement ) );
  FaceBase_sptr face_ptr(ptr);

  for(int num_renewface(0); num_renewface < RENEWLOOP_UPPER_LIMIT; ++num_renewface )
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

    int gate( face_ptr->through_edge(pos_para, newpos) );

    FaceBase_sptr next_face( poly_ptr.lock()->get_neighbor(face_ptr->get_id(), gate) );

    const Real ratio( face_ptr->cross_ratio( pos_para, dis_para, gate) );

    std::pair<Real, Real> temppos( sum( pos_para, multiple(ratio, dis_para) ) );
    std::pair<Real, Real> tempdis( multiple( (1e0 - ratio), dis_para ) );
    assert( on_edge(temppos) );

    std::pair<Real, Real> neighbor_pos( translate_pos(temppos, gate, face_ptr) );
    std::pair<Real, Real> neighbor_dis( translate_dis(tempdis, gate, face_ptr) );
    assert( on_edge(neighbor_pos) );

    pos_para = neighbor_pos;
    dis_para = neighbor_dis;
    face_ptr = next_face;
  }

  throw std::invalid_argument("face ptr renewed over 100 times");
}

std::pair<Real, Real> 
FaceAllGate::translate_pos(std::pair<Real, Real> pos, const int edge_id, FaceBase_sptr face_ptr)
{
  std::pair<Real, Real> n_a_vec( face_ptr->get_para_a_neighbor_at(edge_id) );
  std::pair<Real, Real> n_b_vec( face_ptr->get_para_b_neighbor_at(edge_id) );
  std::pair<Real, Real> n_ori_vec( face_ptr->get_ori_vec_neighbor_at(edge_id) );

  Real nbr_pos_a( pos.first * n_a_vec.first + pos.second * n_b_vec.first + n_ori_vec.first );
  Real nbr_pos_b( pos.first * n_a_vec.second+ pos.second * n_b_vec.second+ n_ori_vec.second);
  
  std::pair<Real, Real> neighbor_pos(nbr_pos_a, nbr_pos_b);
 
  return neighbor_pos;
}

std::pair<Real, Real>
FaceAllGate::translate_dis(std::pair<Real, Real> dis, const int edge_id, FaceBase_sptr face_ptr)
{
  std::pair<Real, Real> n_a_vec( face_ptr->get_para_a_neighbor_at(edge_id) );
  std::pair<Real, Real> n_b_vec( face_ptr->get_para_b_neighbor_at(edge_id) );

  Real nbr_dis_a( dis.first * n_a_vec.first + dis.second * n_b_vec.first);
  Real nbr_dis_b( dis.first * n_a_vec.second+ dis.second * n_b_vec.second);
		
  std::pair<Real, Real> neighbor_dis(nbr_dis_a, nbr_dis_b);
 
  return neighbor_dis;
}

bool
FaceAllGate::in_face(const std::pair<Real, Real>& parameters, const Real tol )
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

int FaceAllGate::through_edge( const std::pair<Real, Real>& position,
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
FaceAllGate::cross_ratio( const std::pair<Real, Real>& position,
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
FaceAllGate::on_edge(const std::pair<Real, Real>& position,
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
FaceAllGate::on_edge(const std::pair<Real, Real>& position,
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


std::pair<Real, Real> FaceAllGate::to_parametric( const Realvec& pos )
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
  std::pair<Real, Real> ret(0e0, 0e0);
  return ret;
}

Realvec FaceAllGate::to_absolute( const std::pair<Real, Real>& parameters )
{
  Realvec absolute_pos( para_a * parameters.first + para_b * parameters.second );
  return absolute_pos;
}

/**************************************** renew position ***************************************/

void FaceAllGate::set_near_vertexs()
{
  for(int i(0); i<3; ++i)
  {
    FaceBase_sptr ptr(poly_ptr.lock()->get_neighbor(face_id, i));
    near_vert_height.push_back( ptr->get_minimum_height( edges.at(i) ) );
  }
  set_neighbors_edge();
  set_neighbors_ori_vtx();
  return;
}


void FaceAllGate::set_neighbors_edge()
{
  for(int i(0); i<3; ++i)
  {
    FaceBase_sptr ptr(poly_ptr.lock()->get_neighbor(face_id, i));

    Realvec neighbor_para_a( ptr->get_para_a() );
    Realvec neighbor_para_b( ptr->get_para_b() );

    Realvec neighbor_normal( ptr->get_normal_vector() );
    Real rot_angle( (-1e0) * acos( dot_product(normal, neighbor_normal) ) );

    Realvec axis( cross_product(normal, neighbor_normal) );
    if( length(axis) == 0e0 )
    {
      //do nothing
    }
    else
    {
      neighbor_para_a = rotation( rot_angle, axis/length(axis), neighbor_para_a );
      neighbor_para_b = rotation( rot_angle, axis/length(axis), neighbor_para_b );
    }

//
    if( fabs( dot_product( neighbor_para_a, normal ) ) > 1e-12 )
      throw std::invalid_argument("rotated neighbor edge but it is not on this face");
    if( fabs( dot_product( neighbor_para_b, normal ) ) > 1e-12 )
      throw std::invalid_argument("rotated neighbor edge but it is not on this face");

    std::pair<Real, Real> neighbor_a( to_parametric( neighbor_para_a ) );
    std::pair<Real, Real> neighbor_b( to_parametric( neighbor_para_b ) );

    // ( neighbor_a.first neighbor_a.second ) (para_a) = (neighbor_para_a)
    // ( neighbor_b.first neighbor_b.second ) (para_b)   (neighbor_para_b)
    Real determinant(neighbor_a.first * neighbor_b.second - neighbor_a.second * neighbor_b.first);
    if(determinant == 0e0)
      throw std::invalid_argument("set_neighbors_edge: determinant is zero");

    // ( neighbor_b.second -neighbor_a_second)
    // ( -neighbor_b.first neighbor_a.first )
    std::pair<Real, Real> nbr_a( neighbor_b.second / determinant, -1e0 * neighbor_a.second / determinant );
    std::pair<Real, Real> nbr_b( -1e0 * neighbor_b.first / determinant, neighbor_a.first / determinant );

    para_a_neighbor.at(i) = nbr_a;
    para_b_neighbor.at(i) = nbr_b;
  }
  return;
}

void FaceAllGate::set_neighbors_ori_vtx()
{
  for(int i(0); i<3; ++i)
  {
    FaceBase_sptr ptr(poly_ptr.lock()->get_neighbor(face_id, i));
    Realvec vi_nbr_ori( ptr->get_para_origin() - vertexs.at(i) );

    Realvec neighbor_normal( ptr->get_normal_vector() );
    Real rot_angle( (-1e0) * acos( dot_product(normal, neighbor_normal) ) );

    Realvec axis( cross_product(normal, neighbor_normal) );
    if( length(axis) == 0e0 )
    {
      //do nothing
    }
    else
    {
      vi_nbr_ori = rotation( rot_angle, axis/length(axis), vi_nbr_ori );
    }

    Realvec nbrori_ori( para_origin - vertexs.at(i) - vi_nbr_ori );
    std::pair<Real, Real> neighbor_ori_this_ori( to_parametric( nbrori_ori ) );

    Real nbr_alpha( neighbor_ori_this_ori.first * para_a_neighbor.at(i).first 
		  + neighbor_ori_this_ori.second * para_b_neighbor.at(i).first );
    Real nbr_beta( neighbor_ori_this_ori.first * para_a_neighbor.at(i).second
		  + neighbor_ori_this_ori.second * para_b_neighbor.at(i).second );
    std::pair<Real, Real> ori( nbr_alpha, nbr_beta );

    ori_vec_neighbor.at(i) = ori;
  }
  return;
}



/*    <- return this vertex
 *  /\
 * /->\
 * ---- <-input edge(neighbors).
 * \<-/
 *  \/ <-call this func using ptr
*/
Realvec FaceAllGate::get_another_vertex(const Realvec& edge)
{
  Realvec tempedge( edge * (-1e0) );

  for(int i(0); i<3; ++i)
  {
    if( edges.at(i)[0] == tempedge[0] &&
	edges.at(i)[1] == tempedge[1] &&
	edges.at(i)[2] == tempedge[2])
    {
      int vertex_id( (i+2) % 3 );
      return vertexs.at( vertex_id );
    }
  }
  throw std::invalid_argument("get_another_vertex: cannot find edge shared by neighboring face");

  Realvec zero;
  return zero;
}

//may redundant.
int FaceAllGate::get_another_vertex_id(const Realvec& edge)
{
  Realvec tempedge( edge * (-1e0) );

  for(int i(0); i<3; ++i)
  {
    if( edges.at(i)[0] == tempedge[0] &&
	edges.at(i)[1] == tempedge[1] &&
	edges.at(i)[2] == tempedge[2])
    {
      return (i+2) % 3;
    }
  }
  throw std::invalid_argument("get_another_vertex_id: cant find edge shared by neighboring face");
  return -1;
}

/*    
 *  /\ <- return this edge id
 * /->\
 * ---- <-input edge(neighbors). so reversed.
 * \<-/
 *  \/ <-neighboring face call this func using ptr
*/
int FaceAllGate::get_another_edge_id_right(const Realvec& neighbors_edge)
{
  Realvec tempedge( neighbors_edge * (-1e0) );

  for(int i(0); i<3; ++i)
  {
    if( edges.at(i)[0] == tempedge[0] &&
	edges.at(i)[1] == tempedge[1] &&
	edges.at(i)[2] == tempedge[2])
    {
      return ( (i+1) % 3 );
    }
  }
  throw std::invalid_argument( "get_another_edge_right: cant find edge shared by neighboring face." );
  return -1;
}

int FaceAllGate::get_another_edge_id_left(const Realvec& neighbors_edge)
{
  Realvec tempedge( neighbors_edge * (-1e0) );

  for(int i(0); i<3; ++i)
  {
    if( edges.at(i)[0] == tempedge[0] &&
	edges.at(i)[1] == tempedge[1] &&
	edges.at(i)[2] == tempedge[2])
    {
      return ((i+2) % 3);
    }
  }
  throw std::invalid_argument( "get_another_edge_left: cant find edge shared by neighboring face." );
  return -1;
}

Realvec FaceAllGate::get_another_vertex_right( const Realvec& neighbors_edge )
{
  int gateway_id( get_another_edge_id_right(neighbors_edge) );
  FaceBase_sptr ptr( poly_ptr.lock()->get_neighbor(face_id, gateway_id) );
  Realvec ret_vertex( ptr->get_another_vertex( edges.at(gateway_id) ) );
  return ret_vertex;
}

Realvec FaceAllGate::get_another_vertex_left( const Realvec& neighbors_edge )
{
  int gateway_id( get_another_edge_id_left(neighbors_edge) );
  FaceBase_sptr ptr( poly_ptr.lock()->get_neighbor(face_id, gateway_id) );
  Realvec ret_vertex( ptr->get_another_vertex( edges.at(gateway_id) ) );
  return ret_vertex;
}

Real FaceAllGate::get_max_a(const Realvec& position, bool& vertex_involve_flag)
{
  vertex_involve_flag = false;
  int size( vertexs.size() );

  Realvec vertvec( vertexs.at(0) - position );
  Real min_distance( length(vertvec) );

  for(int i(1); i<size; ++i)
  {
    vertvec = vertexs.at(i) - position;
    if( min_distance > length( vertvec ) ) min_distance = length( vertvec );
  }
  
  if( !near_vert_height.empty() )
  {
    size = near_vert_height.size();

    for(int i(0); i<size; ++i)
    {
      if( min_distance > near_vert_height.at(i) )
	min_distance = near_vert_height.at(i);
    }
  }

  //TODO
  // if minimal distance from particle to vertex is smaller than this threshold
  // this allows the shell to involve only one vertex.
//   if(min_distance < 1e-4)
//   {
//     vertex_involve_flag = true;
//     Real second_distance;
//     size = vertexs.size();
//
//     for(int i(0); i < size; ++i)
//     {
//       if(i == nearest_face_vertex) continue;
//       vertvec = vertexs.at(i) - position;
//       second_distance = length(vertvec);
//     }
//
//     for(int i(0); i < size; ++i)
//     {
//       if(i == nearest_face_vertex) continue;
//       vertvec = vertexs.at(i) - position;
//       if(second_distance < length(vertvec)) second_distance = length(vertvec);
//     }
//
//     if(near_vertexs.empty()) return second_distance;
//
//     size = near_vertexs.size();
//     for(int i(0); i < size; ++i)
//     {
//       if(i == nearest_neighbor_vertex) continue;
//       vertvec = vertexs.at(i) - position;
//       if(second_distance < length(vertvec)) second_distance = length(vertvec);
//     }
//
//     THROW_UNLESS(std::invalid_argument, second_distance > 1e-4);
//
//     return second_distance;
//   }

  return min_distance;
}

Real FaceAllGate::get_left_angle( const Realvec& neighbors_edge )
{
  Realvec tempedge( neighbors_edge * (-1e0) );
  int number_of_edge(edges.size());

  for(int i(0); i<number_of_edge; ++i)
  {
    if( edges.at(i)[0] == tempedge[0] &&
	edges.at(i)[1] == tempedge[1] &&
	edges.at(i)[2] == tempedge[2])
    {
      return angles.at(i);
    }
  }
  throw std::invalid_argument( "get_left_angle has invalid_input" );
  return 0e0;
}

Real FaceAllGate::get_right_angle( const Realvec& neighbors_edge )
{
  Realvec tempedge( neighbors_edge * (-1e0) );
  int number_of_edge(edges.size());

  for(int i(0); i<number_of_edge; ++i)
  {
    if( edges.at(i)[0] == tempedge[0] &&
	edges.at(i)[1] == tempedge[1] &&
	edges.at(i)[2] == tempedge[2])
    {
      return angles.at( (i+1) % 3 );
    }
  }
  throw std::invalid_argument( "get_right_angle has invalid_input" );
  return 0e0;
}

Real FaceAllGate::pull_neighbors_left_angle( const Realvec& neighbors_edge )
{
  int gateway_id( get_another_edge_id_left(neighbors_edge) );
  FaceBase_sptr ptr( poly_ptr.lock()->get_neighbor(face_id, gateway_id) );
  Real neighbor_left_angle( ptr->get_left_angle( edges.at(gateway_id) ) );
  return neighbor_left_angle;
}

Real FaceAllGate::pull_neighbors_right_angle( const Realvec& neighbors_edge )
{
  int gateway_id( get_another_edge_id_right(neighbors_edge) );
  FaceBase_sptr ptr( poly_ptr.lock()->get_neighbor(face_id, gateway_id) );
  Real neighbor_right_angle( ptr->get_right_angle( edges.at(gateway_id) ) );
  return neighbor_right_angle;
}

Real FaceAllGate::get_minimum_height( const Realvec& neighbors_edge )
{
  Real perpendicular(M_PI * 0.5);
  Real min_height(-1e0);
  //own angle
  Real left_angle( get_left_angle( neighbors_edge ) );
  Real right_angle( get_right_angle( neighbors_edge ) );

  int another_vertex_id( get_another_vertex_id(neighbors_edge) );
  if(left_angle < perpendicular && right_angle < perpendicular)
  {
    Real area( length(cross_product(neighbors_edge, edges.at(another_vertex_id) ) ) );
    min_height =  area / length(neighbors_edge);
  }
  
   left_angle += pull_neighbors_left_angle( neighbors_edge );
  right_angle += pull_neighbors_right_angle( neighbors_edge );

  if(left_angle < perpendicular)
  {
     Realvec left_neighbor_vertex( get_another_vertex_left(neighbors_edge) );
     Realvec edge_leftvert( left_neighbor_vertex - vertexs.at( (another_vertex_id + 2) % 3 ) );
     Real left_length( length( edge_leftvert ) );
     Real left_height( left_length * sin(left_angle) );
     if( min_height > left_height || min_height < 0e0 )
       min_height = left_height;
  }

  if(right_angle < perpendicular )
  {
    Realvec right_neighbor_vertex( get_another_vertex_right(neighbors_edge) );
    Realvec edge_rightvert( right_neighbor_vertex - vertexs.at( (another_vertex_id + 2) % 3 ) );
    Real right_length( length( edge_rightvert ) );
    Real right_height( right_length * sin(right_angle) );
    if( min_height > right_height || min_height < 0e0 )
      min_height = right_height;
  }

  if(min_height == -1e0)
    min_height = length(neighbors_edge) * 0.5;

  return min_height;
}

void FaceAllGate::print_class_name()
{
  std::cout << "class: FaceAllGate" << std::endl;
}

#endif /*FACE_TWO_GATE*/
