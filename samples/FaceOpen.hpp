#ifndef FACE_OPEN
#define FACE_OPEN
#include <boost/weak_ptr.hpp>
#include "FaceBase.hpp"
#include "Polygon.hpp"

class FaceOpen : public FaceBase
{
private:
  std::vector<Realvec> vertexs;
  std::vector<Realvec> near_vertexs;
  std::vector<Realvec> edges;
  std::vector<Real> angles;
  std::vector<bool> is_gate;
  boost::weak_ptr<Polygon> poly_ptr;

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
  FaceOpen(const int id, const Realvec& vtx0, const Realvec& vtx1, const Realvec& vtx2, const std::vector<bool> gate)
  : FaceBase( id, cross_product(vtx1-vtx0, vtx2-vtx0), vtx1-vtx0 ), vertexs(3), edges(3),
    angles(3), is_gate(gate), para_a_neighbor(3), para_b_neighbor(3), ori_vec_neighbor(3)
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

    para_origin = vertexs.at(0);
    para_a = edges.at(0);
    para_b = edges.at(2) * (-1e0);
  }

  FaceOpen(const int id, const Realvec& vtx0, const Realvec& vtx1, const Realvec& vtx2, const Realvec& norm, const std::vector<bool> gate)
  : FaceBase( id, norm, vtx1-vtx0 ), vertexs(3), edges(3), angles(3), is_gate(gate),
    para_a_neighbor(3), para_b_neighbor(3), ori_vec_neighbor(3)
  {
    if( dot_product(norm, vtx1-vtx0) < GLOBAL_TOLERANCE )
      throw std::invalid_argument("FaceOpen: normal vector is not vertical to edge 0");
    if( dot_product(norm, vtx2-vtx1) < GLOBAL_TOLERANCE )
      throw std::invalid_argument("FaceOpen: normal vector is not vertical to edge 1");
  
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

    para_origin = vertexs.at(0);
    para_a = edges.at(0);
    para_b = edges.at(2) * (-1e0);
  }

  virtual void set_poly_ptr( boost::shared_ptr<Polygon>& p_sptr)
  {
    poly_ptr = p_sptr;
  }

  virtual std::pair<Realvec, FaceBase_sptr>
  apply_displacement( const Realvec& position, const Realvec& displacement, const FaceBase_sptr& ptr);

  virtual bool in_face(const std::pair<Real, Real>& parameters, const Real tol = GLOBAL_TOLERANCE);

  virtual int through_edge(const std::pair<Real, Real>& position,
			   const std::pair<Real, Real>& newposition, const Real tol = GLOBAL_TOLERANCE);

  virtual Real cross_ratio(const std::pair<Real, Real>& position,
			   const std::pair<Real, Real>& displacement, const int& edge_num );

  virtual bool on_edge(const std::pair<Real, Real>& position, int& edge_num,
		       const Real tol = GLOBAL_TOLERANCE);

  virtual bool on_edge(const std::pair<Real, Real>& position, const Real tol = GLOBAL_TOLERANCE);

  virtual std::pair<Real, Real> to_parametric( const Realvec& pos );

  virtual Realvec to_absolute( const std::pair<Real, Real>& parameters );

  virtual void set_near_vertexs();

  virtual Realvec get_another_vertex(const Realvec& edge);

  virtual Real get_max_a(const Realvec& position, bool& vertex_involve_flag);

  virtual bool is_gate_at(int edge_id){ return is_gate.at(edge_id); };

  virtual Realvec get_para_origin(){return para_origin;}
  virtual Realvec get_para_a(){return para_a;}
  virtual Realvec get_para_b(){return para_b;}

  virtual std::pair<Real, Real> get_para_a_neighbor_at(const int i){return para_a_neighbor.at(i);}
  virtual std::pair<Real, Real> get_para_b_neighbor_at(const int i){return para_b_neighbor.at(i);}
  virtual std::pair<Real, Real> get_ori_vec_neighbor_at(const int i){return ori_vec_neighbor.at(i);}

  virtual void print_class_name();

private:

  std::pair<Real, Real> 
  translate_pos(std::pair<Real, Real> pos, int edge, FaceBase_sptr face_ptr);

  std::pair<Real, Real> 
  translate_dis(std::pair<Real, Real> dis, int edge, FaceBase_sptr face_ptr);

  std::pair<Real, Real> reverse( const std::pair<Real, Real>& dis, const int edge_num );

  void set_neighbors_edge();

  void set_neighbors_ori_vtx();

};

std::pair<Realvec, FaceBase_sptr>
FaceOpen::apply_displacement(const Realvec& position, const Realvec& displacement,
				const FaceBase_sptr& ptr ) 
{
// debug
//   std::cout << position[0] << " " << position[1] << " " << position[2] << " ";
//   std::cout << displacement[0] <<" "<< displacement[1] << " " << displacement[2] << std::endl;
 
  if( fabs( dot_product( position - para_origin, normal ) ) > 1e-12 )
    throw std::invalid_argument("apply_displacement: position is not on the plane");

  std::pair<Real, Real> pos_para( to_parametric( position - para_origin ) );
  std::pair<Real, Real> dis_para( to_parametric( displacement ) );
  FaceBase_sptr face_ptr(ptr); 

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

    if( face_ptr->is_gate_at(edge_id) )
    {
//       std::cout << "id: " << face_ptr->get_id() << " edge_id: " << edge_id << std::endl;
      FaceBase_sptr next_face( poly_ptr.lock()->get_neighbor(face_ptr->get_id(), edge_id) );

      std::pair<Real, Real> neighbor_pos( translate_pos(temppos, edge_id, face_ptr) );
      std::pair<Real, Real> neighbor_dis( translate_dis(tempdis, edge_id, face_ptr) );
      assert( on_edge(neighbor_pos) );

      pos_para = neighbor_pos;
      dis_para = neighbor_dis;
      face_ptr = next_face;
    }
    else
    {
      tempdis = reverse(tempdis, edge_id);
      pos_para = temppos;
      dis_para = tempdis;
//       std::cout <<"temppos: "<< pos_para.first <<", "<< pos_para.second << " ";
//       std::cout <<"tempdis: "<< dis_para.first <<", "<< dis_para.second << std::endl;
    }
   
  }

  throw std::invalid_argument("face ptr renewed over 100 times");
}

std::pair<Real, Real> 
FaceOpen::translate_pos(std::pair<Real, Real> pos, const int edge_id, FaceBase_sptr face_ptr)
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
FaceOpen::translate_dis(std::pair<Real, Real> dis, const int edge_id, FaceBase_sptr face_ptr)
{
  std::pair<Real, Real> n_a_vec( face_ptr->get_para_a_neighbor_at(edge_id) );
  std::pair<Real, Real> n_b_vec( face_ptr->get_para_b_neighbor_at(edge_id) );

  Real nbr_dis_a( dis.first * n_a_vec.first + dis.second * n_b_vec.first);
  Real nbr_dis_b( dis.first * n_a_vec.second+ dis.second * n_b_vec.second);
		
  std::pair<Real, Real> neighbor_dis(nbr_dis_a, nbr_dis_b);
 
  return neighbor_dis;
}

bool
FaceOpen::in_face(const std::pair<Real, Real>& parameters, const Real tol )
{
  Real alpha(parameters.first);
  Real beta(parameters.second);

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
FaceOpen::through_edge(const std::pair<Real, Real>& position,
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
  if( 0e0 <= newpos_alpha && newpos_beta < 0e0 && newpos_alpha + newpos_beta<=1e0 )
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
FaceOpen::cross_ratio( const std::pair<Real, Real>& position,
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
FaceOpen::on_edge(const std::pair<Real, Real>& position,
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
FaceOpen::on_edge(const std::pair<Real, Real>& position,
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

std::pair<Real, Real> FaceOpen::reverse( const std::pair<Real, Real>& dis, const int edge_num )
{

  switch(edge_num)
  {
    case 0:
      {
	Real k( 2e0 * cos( angles.at(0) ) );
	Real a_new( dis.first + k * dis.second * length(para_b) / length(para_a) );
	Real b_new(-dis.second );
	std::pair<Real, Real> reversed(a_new, b_new);
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
	Real k( 2e0 * cos( angles.at(0) ) );
	Real a_new(-dis.first);
	Real b_new(dis.second + k * dis.first * length(para_a) / length(para_b));
	std::pair<Real, Real> reversed(a_new, b_new);
	return reversed;
      }
      break;

    default:
      throw std::invalid_argument("reverse: invalid edge_num");
      break;
  }
  std::cout << "edge: " << edge_num << " dis: " << dis.first << ", " << dis.second << std::endl;
  throw std::invalid_argument("reverse: switch passed");
}


std::pair<Real, Real> FaceOpen::to_parametric( const Realvec& pos )
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

Realvec FaceOpen::to_absolute( const std::pair<Real, Real>& parameters )
{
  Realvec absolute_pos( para_a * parameters.first + para_b * parameters.second );
  return absolute_pos;
}

/************ end renew_position *************/

void FaceOpen::set_near_vertexs()
{
  for(int i(0); i<3; ++i)
  {
    if( !is_gate.at(i) ) continue;
    FaceBase_sptr ptr(poly_ptr.lock()->get_neighbor(face_id, i));
//     near_vert_height.push_back( ptr->get_minimum_height( edges.at(i) ) );
  }
  set_neighbors_edge();
  set_neighbors_ori_vtx();
  return;
}


void FaceOpen::set_neighbors_edge()
{
  for(int i(0); i<3; ++i)
  {
    if( !is_gate.at(i) ) continue;
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

void FaceOpen::set_neighbors_ori_vtx()
{
  for(int i(0); i<3; ++i)
  {
    if( !is_gate.at(i) ) continue;
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

Realvec FaceOpen::get_another_vertex(const Realvec& edge)
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

//   std::cout << "face id:" << face_id << std::endl;
//   std::cout << tempedge << std::endl;
//   std::cout << edges.at(0) << " ";
//   std::cout << edges.at(1) << " ";
//   std::cout << edges.at(2) << std::endl;

  bool get_another_vertex_have_invalid_edge_input(false);
  THROW_UNLESS( std::invalid_argument, get_another_vertex_have_invalid_edge_input );

  Realvec zero;
  return zero;
}

Real FaceOpen::get_max_a(const Realvec& position, bool& vertex_involve_flag)
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
  if(min_distance < 1e-6)
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
      if(second_distance > length(vertvec)) second_distance = length(vertvec);
    }

    if(near_vertexs.empty()) return second_distance;

    size = near_vertexs.size();
    for(int i(0); i < size; ++i)
    {
      if(i == nearest_neighbor_vertex) continue;
      vertvec = vertexs.at(i) - position;
      if(second_distance > length(vertvec)) second_distance = length(vertvec);
    }

    THROW_UNLESS(std::invalid_argument, second_distance > 1e-6);

    return second_distance;
  }

  return min_distance;
}

void FaceOpen::print_class_name()
{
  std::cout << "class: FaceOpen" << std::endl;
}

#endif /*FACE_OPEN*/
