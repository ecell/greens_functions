#ifndef FACE_ONE_GATE
#define FACE_ONE_GATE
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

class FaceOneGate : public FaceBase
{
  std::vector<Realvec> vertexs;
  std::vector<Realvec> near_vertexs;
  std::vector<Realvec> edges;
  std::vector<Real> angles;
  std::vector<bool> is_gate;
  boost::weak_ptr<Polygon> poly_ptr;


public:
  FaceOneGate(const int& id, const Realvec& vtx0, const Realvec& vtx1, const Realvec& vtx2, const int& gate)
  : FaceBase( id, cross_product(vtx1-vtx0, vtx2-vtx0), vtx1-vtx0 ),
    vertexs(3), edges(3), angles(3), is_gate(3)
  {
    vertexs.at(0) = vtx0;
    vertexs.at(1) = vtx1;
    vertexs.at(2) = vtx2;

    edges.at(0) = vtx1 - vtx0;
    edges.at(1) = vtx2 - vtx1;
    edges.at(2) = vtx0 - vtx2;

    for(int i(0); i<3; ++i)
      is_gate.at(i) = (i == gate);

    angles.at(0) = smaller_angle( edges.at(0), edges.at(2)*(-1e0) );
    angles.at(1) = smaller_angle( edges.at(1), edges.at(0)*(-1e0) );
    angles.at(2) = smaller_angle( edges.at(2), edges.at(1)*(-1e0) );

    para_origin = vertexs.at(0);
    para_a = edges.at(0);
    para_b = edges.at(2) * (-1e0);
  }

  FaceOneGate(const int& id, const Realvec& vtx0, const Realvec& vtx1, const Realvec& vtx2, const Realvec& norm, const int& gate)
  : FaceBase( id, norm, vtx1-vtx0 ), is_gate(3)
  {
    bool normal_vec_is_oritented_orthogonally_to_the_edges( dot_product(norm, vtx1-vtx0) == 0 && dot_product(norm, vtx2-vtx1) == 0 );
    THROW_UNLESS(std::invalid_argument, normal_vec_is_oritented_orthogonally_to_the_edges );
  
    vertexs.at(0) = vtx0;
    vertexs.at(1) = vtx1;
    vertexs.at(2) = vtx2;

    edges.at(0) = vtx1 - vtx0;
    edges.at(1) = vtx2 - vtx1;
    edges.at(2) = vtx0 - vtx2;

    for(int i(0); i<3; ++i)
      is_gate.at(i) = (i == gate);

    angles.at(0) = smaller_angle( edges.at(0), edges.at(2)*(-1e0) );
    angles.at(1) = smaller_angle( edges.at(1), edges.at(0)*(-1e0) );
    angles.at(2) = smaller_angle( edges.at(2), edges.at(1)*(-1e0) );

    para_origin = vertexs.at(0);
    para_a = edges.at(0);
    para_b = edges.at(2) * (-1e0);
  }

  virtual void set_poly_ptr( boost::shared_ptr<Polygon> p_sptr)
  {
    poly_ptr = p_sptr;
  };

//   virtual Realvec renew_position( const Realvec& position, const Realvec& displacement, boost::shared_ptr<FaceBase>& ptr);

  virtual std::pair<Realvec, FaceBase_sptr>
  apply_displacement( const Realvec& position, const Realvec& displacement,
		      const FaceBase_sptr& ptr);

  virtual std::pair<Real, Real> to_parametric( const Realvec& pos );

  virtual bool in_face(const std::pair<Real, Real>& parameters, const Real tol = 1e-12);

  virtual int through_edge(const std::pair<Real, Real>& position,
			   const std::pair<Real, Real>& newposition, const Real tol = 1e-12);

  virtual Real cross_ratio(const std::pair<Real, Real>& position,
			   const std::pair<Real, Real>& displacement, const int& edge_num );

  virtual bool on_edge(const std::pair<Real, Real>& position, int& edge_num,
		       const Real tol = 1e-12);
  virtual bool on_edge(const std::pair<Real, Real>& position, const Real tol = 1e-12);

  virtual Realvec to_absolute( const std::pair<Real, Real>& parameters );

  virtual bool in_the_face( const Realvec& position, const Realvec& displacement );

  virtual void set_near_vertexs();

  virtual Realvec get_another_vertex(const Realvec& edge);

  virtual Realvec get_vertex(){ return para_origin; };

  virtual Real get_max_a(const Realvec& position, bool& vertex_involve_flag);

  virtual void print_class_name();

private:

  std::pair<Real, Real> reverse( const std::pair<Real, Real>& dis, const int edge_num );

};

/*
Realvec FaceOneGate::renew_position
  (const Realvec& position, const Realvec& displacement, boost::shared_ptr<FaceBase>& ptr)
{
  Realvec temppos(position);
  Realvec tempdis(displacement);
//   std::cout << temppos[0] << " " << temppos[1] << " " << temppos[2] << " ";
//   std::cout << tempdis[0] << " " << tempdis[1] << " " << tempdis[2] << std::endl;
  
  THROW_UNLESS(std::invalid_argument, fabs(dot_product(position-para_origin, normal))<1e-12);

  Real newpos_alpha(0e0), newpos_beta(0e0);
  to_parametric(temppos + tempdis - para_origin, newpos_alpha, newpos_beta);

  if( in_this_face( newpos_alpha, newpos_beta ) )
    return para_origin + to_absolute(newpos_alpha, newpos_beta);

  Real pos_alpha(0e0), pos_beta(0e0);
  to_parametric(temppos - para_origin, pos_alpha, pos_beta);

  Real dis_alpha(0e0), dis_beta(0e0);
  to_parametric(tempdis, dis_alpha, dis_beta);
  //dis_alpha( newpos_alpha - pos_alpha ), dis_beta( newpos_beta - pos_beta ) ?

  int crossing_edge( through_edge(pos_alpha, pos_beta, newpos_alpha, newpos_beta) );

  if( is_gate.at(crossing_edge) )
  {
    ptr = poly_ptr.lock()->get_neighbor_ptr_from_gateway( face_id, crossing_edge );
    
    const Real ratio( cross_ratio( pos_alpha, pos_beta, dis_alpha, dis_beta, crossing_edge ) );
    temppos = para_origin 
	    + to_absolute( pos_alpha + ratio*dis_alpha, pos_beta + ratio*dis_beta );
    tempdis = to_absolute( (1e0 - ratio) * dis_alpha, (1e0 - ratio) * dis_beta );

    Realvec neighbor_norm( ptr->get_normal_vector() );

    Real rot_angle( acos( dot_product(normal, neighbor_norm) ) );

    Realvec axis( cross_product(normal, neighbor_norm) );
    Real axis_length( length( axis ) );
    if( axis_length == 0e0 )
    {
      axis = edges.at(crossing_edge) / length( edges.at(crossing_edge) );
    }
    else
    {
      axis = axis / axis_length;
    }

    Realvec new_disp( rotation( rot_angle, axis, tempdis ) );

    if( length(new_disp) >= length(displacement) + 1e-12 )
    {
      std::cout << "ratio: " << ratio << std::endl;
      std::cout << "displacement: " << displacement << " length: " << length(displacement) << std::endl;
      std::cout << "new_disp: " << new_disp << " length: " << length(new_disp) << std::endl;
      throw std::invalid_argument("new displacement is not less than first displacement");
    }
    
  //   if( length(new_disp) < 1e-16*length(temppos) )
  //     return temppos;

    temppos = ptr->renew_position( temppos, new_disp, ptr );
    return temppos;
  }
  else
  {
      Real ratio( cross_ratio( pos_alpha, pos_beta, dis_alpha, dis_beta, crossing_edge ) );

      Real newdis_alpha( (1e0 - ratio) * dis_alpha ), newdis_beta( (1e0 - ratio) * dis_beta  );
      Real newpos_alpha( pos_alpha + ratio * dis_alpha ), newpos_beta( pos_beta + ratio * dis_beta );

      reverse( newdis_alpha, newdis_beta, crossing_edge );

      if( in_this_face( newpos_alpha + newdis_alpha, newpos_beta + newdis_beta ) )
      {
	return to_absolute(newpos_alpha + newdis_alpha, newpos_beta + newdis_beta) + para_origin;
      }
      else
      {
	temppos = to_absolute(newpos_alpha, newpos_beta) + para_origin;
	tempdis = to_absolute(newdis_alpha, newdis_beta);
	temppos = renew_position( temppos, tempdis, ptr );
	return temppos;
      }
  }

  throw std::invalid_argument("renew: all if block did not return");
  Realvec zero;
  return zero;
}*/

// TODO
std::pair<Realvec, FaceBase_sptr>
FaceOneGate::apply_displacement(const Realvec& position, const Realvec& displacement,
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

bool FaceOneGate::in_face(const std::pair<Real, Real>& parameters, const Real tol )
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


int
FaceOneGate::through_edge( const std::pair<Real, Real>& position,
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
FaceOneGate::cross_ratio( const std::pair<Real, Real>& position,
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
FaceOneGate::on_edge(const std::pair<Real, Real>& position,
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
FaceOneGate::on_edge(const std::pair<Real, Real>& position,
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

bool FaceOneGate::in_the_face( const Realvec& position, const Realvec& displacement )
{
  std::pair<Real, Real> para(to_parametric( position+displacement ) );
  Real alpha(para.first), beta(para.second);

  bool just_on_vertex( (alpha == 0e0 || alpha == 1e0) &&
		     ( beta == 0e0 ||  beta == 1e0) );
  if( just_on_vertex )
    throw std::invalid_argument( "in_the_face : particle is just on a vertex." );

  bool alpha_in_the_range( 0e0 <= alpha && alpha <= 1e0 );
  bool  beta_in_the_range( 0e0 <= beta  && beta  <= 1e0 );
  bool sum_in_the_range( 0e0 < alpha + beta && alpha + beta <= 1e0 );

  return (alpha_in_the_range && beta_in_the_range) && sum_in_the_range;
};


std::pair<Real, Real> FaceOneGate::reverse( const std::pair<Real, Real>& dis, const int edge_num )
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
};

std::pair<Real, Real> FaceOneGate::to_parametric( const Realvec& pos )
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

Realvec FaceOneGate::to_absolute( const std::pair<Real, Real>& parameters )
{
  Realvec absolute_pos( para_a * parameters.first + para_b * parameters.second );
  return absolute_pos;
}


/************ renew_position *************/

Real FaceOneGate::get_max_a(const Realvec& position, bool& vertex_involve_flag)
{
  vertex_involve_flag = false;
  int nearest_face_vertex(-1);
  int nearest_neighbor_vertex(-1);

  int size( vertexs.size() );

  Realvec vertvec( vertexs.at(0) - position );
  nearest_face_vertex = 0;

  Real min_distance( length(vertvec) );

  for(int i(1); i<size; ++i)
  {
    vertvec = vertexs.at(i) - position;
    if( min_distance > length( vertvec ) )
    {
      min_distance = length( vertvec );
      nearest_face_vertex = i;
    }
  }
  
  if( !near_vertexs.empty() )
  {
    size = near_vertexs.size();
    for(int i(0); i<size; ++i)
    {
      vertvec = near_vertexs.at(i) - position;
      // length(vertvec) is distance of the straight line between particle position and near vertex
      // but it is lesser equal to distance along the surface
      if( min_distance > length( vertvec ) )
      {
	min_distance = length( vertvec );
	nearest_neighbor_vertex = i;
      }
    }
  }

  // if minimal distance from particle to vertex is smaller than this threshold
  // this allows the shell to involve only one vertex.
  if(min_distance < 1e-4)
  {
    vertex_involve_flag = true;
    Real second_distance;
    size = vertexs.size();

    for(int i(0); i < size; ++i)
    {
      if(i == nearest_face_vertex) continue;
      vertvec = vertexs.at(i) - position;
      second_distance = length(vertvec);
    }

    for(int i(0); i < size; ++i)
    {
      if(i == nearest_face_vertex) continue;
      vertvec = vertexs.at(i) - position;
      if(second_distance < length(vertvec)) second_distance = length(vertvec);
    }

    if(near_vertexs.empty()) return second_distance;

    size = near_vertexs.size();
    for(int i(0); i < size; ++i)
    {
      if(i == nearest_neighbor_vertex) continue;
      vertvec = vertexs.at(i) - position;
      if(second_distance < length(vertvec)) second_distance = length(vertvec);
    }

    THROW_UNLESS(std::invalid_argument, second_distance > 1e-4);
  
    return second_distance;
  }

  return min_distance;
}



void FaceOneGate::set_near_vertexs()
{
  for(int i(0); i<3; ++i)
  {
    if( !is_gate.at(i) ) continue;

    Realvec midpoint( ( vertexs.at(i) + vertexs.at( (i+1)%3 ) ) / 2e0 );
    Real r( length( edges.at(i) ) / 2e0 );

    //edge(i) is gateway so i is gateway-id
    boost::shared_ptr<FaceBase> ptr( poly_ptr.lock()->get_neighbor_ptr_from_gateway(face_id, i) );
    Realvec near_vertex( ptr->get_another_vertex(edges.at(i)) );
    
    Real l( length( near_vertex - midpoint ) );

    if(l <= r)
    {
      near_vertexs.push_back( near_vertex );
    }

  }
}

Realvec FaceOneGate::get_another_vertex(const Realvec& edge)
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
  THROW_UNLESS( std::invalid_argument, get_another_vertex_have_invalid_edge_input );

  Realvec zero;
  return zero;
}

void FaceOneGate::print_class_name()
{
  std::cout << "class: FaceOneGate" << std::endl;
  return;
}

#endif /*FACE_ONE_GATE*/
