#ifndef FACE_ALL_GATE
#define FACE_ALL_GATE
#include <vector>
#include <cmath>
#include "Polygon.hpp"
#include "rotation.hpp"

using namespace greens_functions;

class FaceAllGate : public FaceBase
{
private:
  std::vector<Realvec> vertexs;
  std::vector<Realvec> edges;
  std::vector<bool> is_gate;
  std::vector<Real> angles;
  
  boost::shared_ptr<Polygon> belonging_polygon;
  std::vector<Real> near_vert_height;
  //parametric
  Realvec represent_vertex;
  Realvec a_vec;
  Realvec b_vec;

public:
  FaceAllGate(const int& id, const Realvec& vtx0, const Realvec& vtx1, const Realvec& vtx2)
  : FaceBase( id, cross_product(vtx1-vtx0, vtx2-vtx0), vtx1-vtx0 ),
    vertexs(3), edges(3), is_gate(3,true), angles(3)
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
    represent_vertex = vertexs.at(0);
    a_vec = edges.at(0);
    b_vec = edges.at(2) * (-1e0);
  }

  FaceAllGate(const int& id, const Realvec& vtx0, const Realvec& vtx1, const Realvec& vtx2, const Realvec& norm)
  : FaceBase( id, norm, vtx1-vtx0 ), vertexs(3), edges(3), is_gate(3,true), angles(3)
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
    represent_vertex = vertexs.at(0);
    a_vec = edges.at(0);
    b_vec = edges.at(2) * (-1e0);

  }

  virtual void set_belonging_polygon( boost::shared_ptr<Polygon> p_sptr)
  {
    belonging_polygon = p_sptr;
  };

  //set near_vert_height.
  virtual void set_near_vertexs();

  virtual Realvec renew_position( const Realvec& position, const Realvec& displacement, boost::shared_ptr<FaceBase>& ptr);

  //though input is absolute vector representation, using parametric representation.
  virtual bool still_in_the_face( const Realvec& position, const Realvec& displacement );//legacy

  //return the vertex such that input edge does not include.
  virtual Realvec get_another_vertex(const Realvec& edge);

  //return represent_vertex.
  virtual Realvec get_vertex(){ return vertexs.at(0); };

  //find max shell size (for greens function)
  virtual Real get_max_a(const Realvec& position, bool& vertex_involve_flag);
  
  //call this func from neighbor face
  virtual Real get_minimum_height( const Realvec& neighbors_edge );

  //return own angle.
  virtual Real get_right_angle( const Realvec& neighbors_edge );
  virtual Real get_left_angle( const Realvec& neighbors_edge );

  virtual void print_class_name();

private:
  //return whether the position represented parameter a&b is in this face.
  bool in_this_face(const Real& alpha, const Real& beta);

  //return whitch edge the displacement(newpos-pos) goes through.
  int through_edge(const Real& pos_alpha, const Real& pos_beta, const Real& newpos_alpha, const Real& newpos_beta);

  //return the ratio(0~1) of a displacement segment not cross the edge.
  Real intersection_ratio(const Real& pos_alpha, const Real& pos_beta, const Real& dis_alpha, const Real& dis_beta, const int& edge_num );

  //return whether the place is on the edge and rewrite edge_num to the edge.
  //  if false, edge_num = -1.
  bool on_edge(const Real& alpha, const Real& beta, int& edge_num);

  // rewrite alpha and beta as parameters of parametric expression of vector pos.
  // this doesnt subtract represent_vertex from pos automatically.
  void to_parametric( const Realvec& pos, Real& alpha, Real& beta );
 
  // return absolute expression translated from parametric expression.
  Realvec to_absolute( const Real& alpha, const Real& beta );

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

Realvec
FaceAllGate::renew_position(const Realvec& position, const Realvec& displacement,
			    boost::shared_ptr<FaceBase>& ptr)
{
  Realvec temppos(position);
  Realvec tempdis(displacement);
// debug
//   std::cout << temppos[0] << " " << temppos[1] << " " << temppos[2] << " ";
//   std::cout << tempdis[0] << " " << tempdis[1] << " " << tempdis[2] << std::endl;
 
  THROW_UNLESS(std::invalid_argument,fabs(dot_product(position-represent_vertex, normal) )<1e-12);

  //parametric
  Real newpos_alpha(0e0), newpos_beta(0e0);
  to_parametric(temppos + tempdis - represent_vertex, newpos_alpha, newpos_beta);

  //after move, particle is in this face
  if( in_this_face( newpos_alpha, newpos_beta ) )
    return represent_vertex + to_absolute(newpos_alpha, newpos_beta);

  //find polygon particle gone
  Real pos_alpha(0e0), pos_beta(0e0);
  to_parametric(temppos - represent_vertex, pos_alpha, pos_beta);

  Real dis_alpha(0e0), dis_beta(0e0);
  to_parametric(tempdis, dis_alpha, dis_beta);
  //dis_alpha( newpos_alpha - pos_alpha ), dis_beta( newpos_beta - pos_beta ) ?

  int gate( through_edge(pos_alpha, pos_beta, newpos_alpha, newpos_beta) );

  ptr = belonging_polygon->get_neighbor_ptr_from_gateway( face_id, gate );
  
  const Real ratio( intersection_ratio( pos_alpha, pos_beta, dis_alpha, dis_beta, gate ) );
  temppos = represent_vertex 
	  + to_absolute( pos_alpha + ratio*dis_alpha, pos_beta + ratio*dis_beta );
  tempdis = to_absolute( (1e0 - ratio) * dis_alpha, (1e0 - ratio) * dis_beta );

  Realvec neighbor_norm( ptr->get_normal_vector() );

  Real rot_angle( acos( dot_product(normal, neighbor_norm) ) );

  Realvec axis( cross_product(normal, neighbor_norm) );
  Real axis_length( length( axis ) );
  if( axis_length == 0e0 )
  {
    axis = edges.at(gate) / length( edges.at(gate) );
  }
  else
  {
    axis = axis / axis_length;
  }

  Realvec new_disp( rotation( rot_angle, axis, tempdis ) );

  if( length(new_disp) > length(displacement) )
  {
    std::cout << "ratio: " << ratio << " 1-ratio: " << 1e0-ratio << std::endl;
    std::cout << "new_disp: (" << new_disp[0] << ", " << new_disp[1] << ", " << new_disp[2] << ")" << std::endl;
    std::cout << "displacement: (" << displacement[0] << ", " << displacement[1] << ", " << displacement[2] << ")" << std::endl;
    throw std::invalid_argument("new displacement is not less than first displacement");
  }
  
//   if( length(new_disp) < 1e-16*length(temppos) )
//     return temppos;

  temppos = ptr->renew_position( temppos, new_disp, ptr );
  return temppos;
}

bool
FaceAllGate::in_this_face(const Real& alpha, const Real& beta)
{
  bool just_on_vertex( (alpha== 0e0 || alpha == 1e0) &&
		       (beta == 0e0 ||  beta == 1e0) );
  if(just_on_vertex)
  {
    std::cout << "alpha: " << std::setprecision(16) << alpha << ", beta: " << beta << std::endl;
    throw std::invalid_argument("in_this_face : particle is just on a vertex");
  }

  bool alpha_in_range( 0e0 <= alpha && alpha <= 1e0 );
  bool beta_in_range( 0e0 <= beta  && beta  <= 1e0 );
  bool sum_in_range( 0e0 < alpha+beta && alpha+beta <= 1e0 );

  return ( (alpha_in_range && beta_in_range) && sum_in_range);
}

int
FaceAllGate::through_edge( const Real& pos_alpha, const Real& pos_beta,
			   const Real& newpos_alpha, const Real& newpos_beta)
{
  //case::first position is on the edge
  int edge_num(-1);
  if( on_edge( pos_alpha, pos_beta, edge_num ) )
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
  if( 0e0 <= newpos_alpha && newpos_beta < 0e0 && newpos_alpha + newpos_beta<=1e0 )
    return 0;
  if( 0e0 < newpos_alpha && 0e0 < newpos_beta && 1e0 < newpos_alpha + newpos_beta )
    return 1;
  if( 0e0 <= newpos_beta && newpos_alpha < 0e0 && newpos_alpha + newpos_beta <= 1e0 )
    return 2;

  Real dis_alpha(newpos_alpha - pos_alpha);
  Real dis_beta(newpos_beta - pos_beta);
  //through two edge: case {3,4,5}
  if( newpos_alpha < 0e0 && newpos_beta < 0e0 )
  {
  // case 3
  //throwgh edges.at(0) and at(2)
    Real ratio0( intersection_ratio(pos_alpha, pos_beta, dis_alpha, dis_beta, 0) );
    Real ratio2( intersection_ratio(pos_alpha, pos_beta, dis_alpha, dis_beta, 2) );
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
    Real ratio0( intersection_ratio(pos_alpha, pos_beta, dis_alpha, dis_beta, 0) );
    Real ratio1( intersection_ratio(pos_alpha, pos_beta, dis_alpha, dis_beta, 1) );
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
    Real ratio1( intersection_ratio(pos_alpha, pos_beta, dis_alpha, dis_beta, 1) );
    Real ratio2( intersection_ratio(pos_alpha, pos_beta, dis_alpha, dis_beta, 2) );
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
FaceAllGate::intersection_ratio( const Real& pos_alpha, const Real& pos_beta,
				 const Real& dis_alpha, const Real& dis_beta,
				 const int& edge_num )
{
  Real ratio;
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
    ratio = 0e0;
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
FaceAllGate::on_edge(const Real& alpha, const Real& beta,
		     int& edge_num)
{
  Real eps(1e-10);

  if( alpha == 0e0 && (0e0 <= beta && beta <= 1e0) )
  {
    edge_num = 2;
    return true;
  }
  else if( beta == 0e0 && ( 0e0 <= alpha && alpha <= 1e0 ) )
  {
    edge_num = 0;
    return true;
  }
  else if( (1e0-eps < alpha+beta && alpha+beta < 1e0+eps) &&
	   (0e0 < alpha && 0e0 < beta) )
  {
    edge_num = 1;
    return true;
  }else
  {
    edge_num = -1;
    return false;
  }
}


void FaceAllGate::to_parametric( const Realvec& pos, Real& alpha, Real& beta )
{
  Realvec parametric_pos( pos );

  if( (a_vec[0]*b_vec[1] - b_vec[0]*a_vec[1]) != 0e0 )
  {
  // matrix( a0 b0 ) inverse( b1 -b0 ) pos( pos0 )
  //       ( a1 b1 )        ( -a1 a0 )    ( pos1 )
    alpha = ( b_vec[1] * parametric_pos[0] - b_vec[0] * parametric_pos[1] )
	    / (a_vec[0]*b_vec[1] - b_vec[0]*a_vec[1]);

     beta = ( a_vec[0] * parametric_pos[1] - a_vec[1] * parametric_pos[0] )
	    / (a_vec[0]*b_vec[1] - b_vec[0]*a_vec[1]);
 
    //confirm (alpha beta) * (a2 b2) = (pos2)
    if( fabs(alpha * a_vec[2] + beta * b_vec[2] - parametric_pos[2]) > 1e-12 ) 
      throw std::invalid_argument( "function to_parametric: invalid solution" );

  }else if( (a_vec[1]*b_vec[2] - b_vec[1]*a_vec[2]) != 0e0 )
  {
  // matrix( a1 b1 ) inverse( b2 -b1 ) pos( pos1 )
  //       ( a2 b2 )        ( -a2 a1 )    ( pos2 )
    alpha = ( b_vec[2] * parametric_pos[1] - b_vec[1] * parametric_pos[2] )
	    / (a_vec[1]*b_vec[2] - b_vec[1]*a_vec[2]);

     beta = ( a_vec[1] * parametric_pos[2] - a_vec[2] * parametric_pos[1] )
	    / (a_vec[1]*b_vec[2] - b_vec[1]*a_vec[2]);
 
    //confirm (alpha beta) * (a0 b0) = (pos0)
    if( fabs(alpha * a_vec[0] + beta * b_vec[0] - parametric_pos[0]) > 1e-12 )
      throw std::invalid_argument( "function to_parametric: invalid solution" );

  }else if( (a_vec[2]*b_vec[0] - b_vec[2]*a_vec[0]) != 0e0 )
  {
  // matrix( a2 b2 ) inverse( b0 -b2 ) pos( pos2 )
  //       ( a0 b0 )        ( -a0 a2 )    ( pos0 )
    alpha = ( b_vec[0] * parametric_pos[2] - b_vec[2] * parametric_pos[0] )
	    / (a_vec[2]*b_vec[0] - b_vec[2]*a_vec[0]);

     beta = ( a_vec[2] * parametric_pos[0] - a_vec[0] * parametric_pos[2] )
	    / (a_vec[2]*b_vec[0] - b_vec[2]*a_vec[0]);
 
    //confirm (alpha beta) * (a1 b1) = (pos1)
    if( fabs(alpha * a_vec[1] + beta * b_vec[1] - parametric_pos[1]) > 1e-12 ) 
      throw std::invalid_argument( "function to_parametric: invalid solution" ); 

  }else{
    throw std::invalid_argument( "function to_parametric: could not parametrise input position" );
  }

  return;
}

Realvec FaceAllGate::to_absolute( const Real& alpha, const Real& beta )
{
  Realvec absolute_pos( a_vec * alpha + b_vec * beta );
  return absolute_pos;
}

//legacy
bool FaceAllGate::still_in_the_face( const Realvec& position, const Realvec& displacement )
{
  Real alpha(0e0), beta(0e0);
  to_parametric( position+displacement, alpha, beta );

  bool just_on_vertex( (alpha == 0e0 || alpha == 1e0) &&
		     ( beta == 0e0 ||  beta == 1e0) );
  if( just_on_vertex )
    throw std::invalid_argument( "still_in_the_face : particle is just on a vertex." );

  bool alpha_in_the_range( 0e0 <= alpha && alpha <= 1e0 );
  bool  beta_in_the_range( 0e0 <= beta  && beta  <= 1e0 );
  bool sum_in_the_range( 0e0 < alpha + beta && alpha + beta <= 1e0 );

  return (alpha_in_the_range && beta_in_the_range) && sum_in_the_range;
};

/**************************************** renew position ***************************************/

void FaceAllGate::set_near_vertexs()
{
  for(int i(0); i<3; ++i)
  {
    boost::shared_ptr<FaceBase> ptr(belonging_polygon->get_neighbor_ptr_from_gateway(face_id, i));
    near_vert_height.push_back( ptr->get_minimum_height( edges.at(i) ) );
  }
}

/*    <- return this vertex
 *  /\
 * /->\
 * ---- <-input edge(neighbors). so reversed.
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
  boost::shared_ptr<FaceBase> ptr( belonging_polygon->get_neighbor_ptr_from_gateway(face_id, gateway_id) );
  Realvec ret_vertex( ptr->get_another_vertex( edges.at(gateway_id) ) );
  return ret_vertex;
}

Realvec FaceAllGate::get_another_vertex_left( const Realvec& neighbors_edge )
{
  int gateway_id( get_another_edge_id_left(neighbors_edge) );
  boost::shared_ptr<FaceBase> ptr( belonging_polygon->get_neighbor_ptr_from_gateway(face_id, gateway_id) );
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
  boost::shared_ptr<FaceBase> ptr( belonging_polygon->get_neighbor_ptr_from_gateway(face_id, gateway_id) );
  Real neighbor_left_angle( ptr->get_left_angle( edges.at(gateway_id) ) );
  return neighbor_left_angle;
}

Real FaceAllGate::pull_neighbors_right_angle( const Realvec& neighbors_edge )
{
  int gateway_id( get_another_edge_id_right(neighbors_edge) );
  boost::shared_ptr<FaceBase> ptr( belonging_polygon->get_neighbor_ptr_from_gateway(face_id, gateway_id) );
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
