#ifndef FACE_CLOSE
#define FACE_CLOSE
#include <boost/weak_ptr.hpp>
#include "FaceBase.hpp"
#include "Polygon.hpp"

class FaceClose : public FaceBase
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
  //               a*para_a_neighbor.second + b*para_b_neighbor.second)
  std::vector< std::pair<Real, Real> > para_a_neighbor;
  std::vector< std::pair<Real, Real> > para_b_neighbor;
  // para_origin = ori_vec_neighbor.first * neighbor->para_a + 
  //               ori_vec_neighbor.second * neighbor->para_b;
  std::vector< std::pair<Real, Real> > ori_vec_neighbor;


public:
  FaceClose(const int id, const Realvec& vtx0, const Realvec& vtx1, const Realvec& vtx2)
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

  FaceClose(const int& id, const Realvec& vtx0, const Realvec& vtx1, const Realvec& vtx2,
              const Realvec& norm)
  : FaceBase( id, norm, vtx1-vtx0 ), vertexs(3), edges(3), angles(3),
    para_a_neighbor(3), para_b_neighbor(3), ori_vec_neighbor(3)
  {
    if(fabs(dot_product(norm, vtx1-vtx0) ) > GLOBAL_TOLERANCE ||
       fabs(dot_product(norm, vtx2-vtx1) ) > GLOBAL_TOLERANCE )
      throw std::invalid_argument("constructor : normal vector is not vertical");

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
    if(length(edges.at(0)) >= length(edges.at(1)) )
    {
      if(length(edges.at(2)) >= length(edges.at(1)) )
      {//0&2
        para_a = edges.at(0);
        para_b = edges.at(2) * (-1e0);
        para_origin = vertexs.at(0);
      }else{
       //0&1
        para_a = edges.at(0) * (-1e0);
        para_b = edges.at(1);
        para_origin = vertexs.at(1);
      }
    }else{
      if(length(edges.at(2)) >= length(edges.at(0)) )
      {//1&2
        para_a = edges.at(2);
        para_b = edges.at(1) * (-1e0);
        para_origin = vertexs.at(2);
      }else{
       //1&0
        para_a = edges.at(0) * (-1e0);
        para_b = edges.at(1);
        para_origin = vertexs.at(1);
      }
    }
  }

  virtual void set_poly_ptr( boost::shared_ptr<Polygon>& p_sptr)
  {
    poly_ptr = p_sptr;
  }

  virtual std::pair<Realvec, FaceBase_sptr>
  apply_displacement( const Realvec& position, const Realvec& displacement,
                      const FaceBase_sptr& ptr);

  //return whether the position is in a face ( in the range or not ).
  virtual bool in_face(const std::pair<Real, Real>& parameters, const Real tol = GLOBAL_TOLERANCE);

  //return whitch edge the displacement(newpos-pos) goes through.
  virtual int through_edge(const std::pair<Real, Real>& position,
                           const std::pair<Real, Real>& newposition, const Real tol = GLOBAL_TOLERANCE);

  //return the ratio(0~1) of a displacement segment not cross the edge.
  virtual Real cross_ratio(const std::pair<Real, Real>& position,
                           const std::pair<Real, Real>& displacement, const int& edge_num );

  //return whether the place is on the edge and rewrite edge_num to the edge.
  //  if false, edge_num = -1.
  virtual bool on_edge(const std::pair<Real, Real>& position, int& edge_num,
                       const Real tol = GLOBAL_TOLERANCE);

  virtual bool on_edge(const std::pair<Real, Real>& position, const Real tol = GLOBAL_TOLERANCE);

  // return pair of parameters of parametric expression of vector pos.
  // pair.first = alpha, pair.second = beta
  virtual std::pair<Real, Real> to_parametric( const Realvec& pos );
 
  // return absolute expression translated from parametric expression.
  virtual Realvec to_absolute( const std::pair<Real, Real>& parameters );

  virtual std::pair<Real, Real> projection(const Realvec& pos);

  virtual void set_near_vertexs();

//************************************************************//
  //return the vertex such that input edge does not include.
  virtual Realvec get_another_vertex(const Realvec& neighbors_edge);

  //return an id of edge that matches argument neighbors_edge
  int matching_edge(const Realvec& neighbors_edge );

  //find max shell size (for greens function)
  virtual Real get_max_a(const Realvec& position, bool& vertex_include_flag);
  
  //call this func from neighbor face
  virtual Real get_minimum_height( const Realvec& neighbors_edge );

  //return own angle.
  virtual Real get_right_angle( const Realvec& neighbors_edge );
  virtual Real get_left_angle( const Realvec& neighbors_edge );

  virtual Realvec get_vertex_at(const int i){return vertexs.at(i);}
  virtual Real get_angles_at(int i){return angles.at(i);}
  virtual Realvec get_edges_at(int i){return edges.at(i);}

  virtual Realvec get_center_mass()
  {
    Realvec cm(vertexs.at(0) + vertexs.at(1) + vertexs.at(2));
    return (cm / 3e0);
  }

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

  // parametrize para_a & para_b using neighbor's apara_b
  void set_neighbors_edge();
  void set_neighbors_ori_vtx();

  //return vertex id that is not on the edge. 
  //  when neighboring face execute this function through ptr, 
  //  this func return the vertex that is not shared with the neighbor.
  int get_another_vertex_id(const Realvec& neighbors_edge);

  //return own edge id
  int get_another_edge_id_right( const Realvec& neighbors_edge);
  int get_another_edge_id_left( const Realvec& neighbors_edge);

};

std::pair<Realvec, FaceBase_sptr>
FaceClose::apply_displacement(const Realvec& position, const Realvec& displacement,
                                const FaceBase_sptr& ptr )
{
// debug
//   std::cout << position[0] << " " << position[1] << " " << position[2] << " ";
//   std::cout << displacement[0] <<" "<< displacement[1] << " " << displacement[2] << std::endl;

  if( fabs( dot_product( position - para_origin, normal ) ) > GLOBAL_TOLERANCE )
  {
    std::cout << "dot product: " << dot_product( position - para_origin, normal ) << std::endl;
    throw std::invalid_argument("apply_displacement: position is not on the plane");
  }

  std::pair<Real, Real> pos_para( to_parametric( position - para_origin ) );
  std::pair<Real, Real> dis_para( to_parametric( displacement ) );
  FaceBase_sptr face_ptr(ptr);

  for(int num_renewface(0); num_renewface < RENEWLOOP_UPPER_LIMIT; ++num_renewface )
  {
    if(!in_face(pos_para))
//       throw std::invalid_argument("apply_displacement: position is not in the face");
    {
      std::cout << "in_face(pos_para)" << std::endl;
      std::cout << "now renew " << num_renewface << " times." << std::endl;
      std::cout << "in_face(pos_para):  " << in_face(pos_para) << std::endl;
      std::cout << "pos_para:          " << std::setprecision(16) << "(" << pos_para.first 
                << ", " << pos_para.second << ")" << std::endl;
      std::cout << "dis_para:          " << std::setprecision(16) << "(" << dis_para.first
                << ", " << dis_para.second << ")" << std::endl;
      throw std::invalid_argument("apply_displacement: position is not in the face");
    }

    std::pair<Real, Real> newpos( sum(pos_para, dis_para) );

    if( in_face( newpos ) )
    {
      Realvec renewed_pos( face_ptr->get_para_origin() + face_ptr->to_absolute( newpos ) );
      std::pair<Realvec, FaceBase_sptr> retpos( renewed_pos, face_ptr );
      return retpos;
    }

    int gate( face_ptr->through_edge(pos_para, newpos) );

    FaceBase_sptr next_face( poly_ptr.lock()->get_neighbor(face_ptr->get_id(), gate) );

    Real ratio( face_ptr->cross_ratio( pos_para, dis_para, gate) );

    std::pair<Real, Real> temppos( sum( pos_para, multiple(ratio, dis_para) ) );
    std::pair<Real, Real> tempdis( multiple( (1e0 - ratio), dis_para ) );

    if(on_vertex(temppos))
    {
      std::cout << "on_vertex(temppos)" << std::endl;
      std::cout << "now renew " << num_renewface << " times." << std::endl;
      std::cout << "ratio:              " << std::setw(20) << ratio << std::endl;
      std::cout << "in_face(pos_para):  " << in_face(pos_para) << std::endl;
      std::cout << "pos_para:          " << std::setprecision(16) << "(" << pos_para.first 
                << ", " << pos_para.second << ")" << std::endl;
      std::cout << "dis_para:          " << std::setprecision(16) << "(" << dis_para.first
                << ", " << dis_para.second << ")" << std::endl;
      std::cout << "temppos:           " << std::setprecision(16) << "(" << temppos.first
                << ", " << temppos.second << ")" << std::endl;
      throw std::invalid_argument("shortened displacement to put the end of vector is on the edge but it is on vertex" );
    }
//     assert( on_edge(temppos) );
    if( !on_edge(temppos) )
    {
      std::cout << "!on_edge(temppos)" << std::endl;
      std::cout << "now renew " << num_renewface << " times." << std::endl;
      std::cout << "ratio:              " << ratio << std::endl;
      std::cout << "in_face(pos_para):  " << in_face(pos_para) << std::endl;
      std::cout << "pos_para:          " << "(" << std::setprecision(16) << pos_para.first
                << ", " << std::setprecision(16) << pos_para.second << ")" << std::endl;
      std::cout << "dis_para:          " << "(" << std::setprecision(16) << dis_para.first
                << ", " << std::setprecision(16) << dis_para.second << ")" << std::endl;
      std::cout << "temppos:           " << "(" << std::setprecision(16) << temppos.first 
                << ", "  << std::setprecision(16) << temppos.second << ")" << std::endl;
      throw std::invalid_argument("shortened displacement to put the end of vector is on the edge but it is not on edge");
    }

    std::pair<Real, Real> neighbor_pos( translate_pos(temppos, gate, face_ptr) );
    std::pair<Real, Real> neighbor_dis( translate_dis(tempdis, gate, face_ptr) );
    if( !on_edge(neighbor_pos) )
    {
      std::cout << "neighbor_pos: ";
      std::cout << "(" << std::setprecision(16) << neighbor_pos.first << ", " << neighbor_pos.second << ")" << std::endl;
      std::cout << "face id: " << next_face->get_id() << std::endl;
      std::cout << "vertex 0: " << next_face->get_vertex_at(0) << std::endl;
      std::cout << "vertex 1: " << next_face->get_vertex_at(1) << std::endl;
      std::cout << "vertex 2: " << next_face->get_vertex_at(2) << std::endl;

      std::pair<Real, Real> temp;
      temp = face_ptr->get_para_a_neighbor_at(gate);
      std::cout << " n_a_vec: (" << temp.first <<", "<< temp.second << ")" << std::endl;
      temp = face_ptr->get_para_b_neighbor_at(gate);
      std::cout << " n_b_vec: (" << temp.first <<", "<< temp.second << ")" << std::endl;
      temp = face_ptr->get_ori_vec_neighbor_at(gate);
      std::cout << " n_ori_vec: (" << temp.first <<", "<< temp.second << ")" << std::endl;
      throw std::invalid_argument("neighbor_pos is not on edge");
    }

    pos_para = neighbor_pos;
    dis_para = neighbor_dis;
    face_ptr = next_face;
  }

  throw std::invalid_argument("face ptr renewed over RENEWLOOP_UPPER_LIMIT times");
}

std::pair<Real, Real> 
FaceClose::translate_pos(std::pair<Real, Real> pos, const int edge_id, FaceBase_sptr face_ptr)
{
  std::pair<Real, Real> n_a_vec( face_ptr->get_para_a_neighbor_at(edge_id) );
  std::pair<Real, Real> n_b_vec( face_ptr->get_para_b_neighbor_at(edge_id) );
  std::pair<Real, Real> n_ori_vec( face_ptr->get_ori_vec_neighbor_at(edge_id) );

  Real nbr_pos_a( pos.first * n_a_vec.first + pos.second * n_b_vec.first + n_ori_vec.first );
  Real nbr_pos_b( pos.first * n_a_vec.second+ pos.second * n_b_vec.second+ n_ori_vec.second);
 
  return std::make_pair(nbr_pos_a, nbr_pos_b);
}

std::pair<Real, Real>
FaceClose::translate_dis(std::pair<Real, Real> dis, const int edge_id, FaceBase_sptr face_ptr)
{
  std::pair<Real, Real> n_a_vec( face_ptr->get_para_a_neighbor_at(edge_id) );
  std::pair<Real, Real> n_b_vec( face_ptr->get_para_b_neighbor_at(edge_id) );

  Real nbr_dis_a( dis.first * n_a_vec.first + dis.second * n_b_vec.first);
  Real nbr_dis_b( dis.first * n_a_vec.second+ dis.second * n_b_vec.second);

  return std::make_pair(nbr_dis_a, nbr_dis_b);
}

bool
FaceClose::in_face( const std::pair<Real, Real>& parameters, const Real tol )
{
  Real alpha( parameters.first );
  Real beta( parameters.second );

  bool just_on_vertex( (alpha== 0e0 || alpha == 1e0) && (beta == 0e0 ||  beta == 1e0) );
  if(just_on_vertex)
  {
//     std::cout << "Warning: position is just on the vertex." << std::endl;
    std::cout << "alpha: " << std::setprecision(16) << alpha << ", beta: " << beta << std::endl;
    throw std::invalid_argument("in_face : particle is just on a vertex");
//     return true;
  }

  bool alpha_in_range( -tol <= alpha && alpha <= 1e0 + tol );
  bool beta_in_range( -tol <= beta  && beta  <= 1e0 + tol );
  bool sum_in_range( -tol < alpha+beta && alpha+beta <= 1e0 + tol );

  return ( (alpha_in_range && beta_in_range) && sum_in_range);
}

int FaceClose::through_edge( const std::pair<Real, Real>& position,
                             const std::pair<Real, Real>& newposition, const Real tol)
{
  assert(in_face(position));

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

Real FaceClose::cross_ratio( const std::pair<Real, Real>& position,
                             const std::pair<Real, Real>& displacement,
                             const int& edge_num )
{
  if(!in_face(position))
    throw std::invalid_argument("cross: position is not in the face");
  Real pos_alpha( position.first );
  Real pos_beta( position.second );
  Real dis_alpha( displacement.first );
  Real dis_beta( displacement.second );

  switch(edge_num)
  {
  case 0:
    if( dis_beta == 0e0 )
      throw std::invalid_argument("cannot intersect this edge");
    return (-pos_beta / dis_beta);

  case 1:
    if( dis_alpha + dis_beta == 0e0 )
      throw std::invalid_argument("cannot intersect this edge");
    return ( ( 1e0 - pos_alpha - pos_beta ) / ( dis_alpha + dis_beta ) );

  case 2:
    if( dis_alpha == 0e0 )
      throw std::invalid_argument("cannot intersect this edge");
    return (-pos_alpha / dis_alpha);

  default:
    throw std::invalid_argument("invalid edge_num");
  }
  //ratio must be in (0, 1] range.
//   if( ratio < 0e0 || 1e0 < ratio )
//     throw std::invalid_argument("ratio is not in range");
}

/*     2
//    /\
// b /  \
// ^/____\
// 0 ->a  1
*/
bool
FaceClose::on_edge(const std::pair<Real, Real>& position,
                     int& edge_num, const Real tol)
{
  Real alpha(position.first);
  Real beta(position.second);

  if( fabs(alpha) < tol && (-tol <= beta && beta <= 1e0+tol) )
  {
    edge_num = 2;
    return true;
  }
  else if( fabs(beta) < tol && ( -tol <= alpha && alpha <= 1e0+tol ) )
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
FaceClose::on_edge(const std::pair<Real, Real>& position,
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


std::pair<Real, Real> FaceClose::to_parametric( const Realvec& pos )
{
  Real alpha, beta;
  Realvec parametric_pos( pos );

  Real determinant = para_a[0]*para_b[1] - para_b[0]*para_a[1];
  if( determinant != 0e0 )
  {
  // matrix( a0 b0 )-1 = (b1 -b0)  /
  //       ( a1 b1 )     (-a1 a0) / det
    alpha = ( para_b[1] * parametric_pos[0] - para_b[0] * parametric_pos[1] ) / determinant;
     beta = ( para_a[0] * parametric_pos[1] - para_a[1] * parametric_pos[0] ) / determinant;
 
    //confirm (alpha beta) * (a2 b) = (pos2)
    if( fabs(alpha * para_a[2] + beta * para_b[2] - parametric_pos[2]) > GLOBAL_TOLERANCE )
    {
      for(int i(0); i<3; ++i)
      std::cout <<"vertexs.at(" << i << "): " << vertexs.at(i) << std::endl;
      std::cout << "alpha: " << alpha << " beta: " << beta << std::endl;
      std::cout << "position: " << pos << std::endl;
      std::cout << "value: " << alpha * para_a[2] + beta * para_b[2] - parametric_pos[2] << std::endl;
//       throw std::invalid_argument( "function to_parametric: invalid solution" );
    }
    
    return std::make_pair(alpha, beta);
  }
  
  determinant = para_a[1]*para_b[2] - para_b[1]*para_a[2];
  if( determinant != 0e0 )
  {
  // matrix( a1 b1 )
  //       ( a2 b2 )
    alpha = ( para_b[2] * parametric_pos[1] - para_b[1] * parametric_pos[2] ) / determinant;
     beta = ( para_a[1] * parametric_pos[2] - para_a[2] * parametric_pos[1] ) / determinant;
 
    //confirm (alpha beta) * (a0 b0) = (pos0)
    if( fabs(alpha * para_a[0] + beta * para_b[0] - parametric_pos[0]) > GLOBAL_TOLERANCE )
    {
      for(int i(0); i<3; ++i)
        std::cout <<"vertexs.at(" << i << "): " << vertexs.at(i) << std::endl;
      std::cout << "alpha: " << alpha << " beta: " << beta << std::endl;
      std::cout << "position: " << pos << std::endl;
      std::cout << "value: " << alpha * para_a[0] + beta * para_b[0] - parametric_pos[0] << std::endl;
//       throw std::invalid_argument( "function to_parametric: invalid solution" );
    }
      
    return std::make_pair(alpha, beta);
  }
  
  determinant = para_a[2]*para_b[0] - para_b[2]*para_a[0];
  if( determinant != 0e0 )
  {
  // matrix( a2 b2 )
  //       ( a0 b0 )
    alpha = ( para_b[0] * parametric_pos[2] - para_b[2] * parametric_pos[0] ) / determinant;
     beta = ( para_a[2] * parametric_pos[0] - para_a[0] * parametric_pos[2] ) / determinant;
 
    //confirm (alpha beta) * (a1 b1) = (pos1)
    if( fabs(alpha * para_a[1] + beta * para_b[1] - parametric_pos[1]) > GLOBAL_TOLERANCE) 
    { 
      for(int i(0); i<3; ++i)
      std::cout <<"vertexs.at(" << i << "): " << vertexs.at(i) << std::endl;
      std::cout << "alpha: " << alpha << " beta: " << beta << std::endl;
      std::cout << "position: " << pos << std::endl;
      std::cout << "value: " << alpha * para_a[1] + beta * para_b[1] - parametric_pos[1] << std::endl;
//       throw std::invalid_argument( "function to_parametric: invalid solution" );
    }
     
    return std::make_pair(alpha, beta);
  }

  throw std::invalid_argument( "function to_parametric: could not parametrise input position" );
}

std::pair<Real, Real> FaceClose::projection(const Realvec& pos)
{
  Realvec _pos(pos);
  Real alpha, beta, gamma;
  Real det_term1( para_a[0] * para_b[1] * normal[2] + para_a[1] * para_b[2] * normal[0] + para_a[2] * para_b[0] * normal[1] );
  Real det_term2( para_a[0] * para_b[2] * normal[1] + para_a[1] * para_b[0] * normal[2] + para_a[2] * para_b[1] * normal[0] );
  Real determinant(det_term1 - det_term2);

  if(determinant == 0e0)
    throw std::invalid_argument("projection: determinant is zero.");
  if(isnan(determinant))
  {
    std::cout << "face id: " << get_id() <<std::endl;
    std::cout << "projection: determinant is NaN" << std::endl;
    std::cout << "vertexs.at(0)" << vertexs.at(0) << std::endl;
    std::cout << "vertexs.at(1)" << vertexs.at(1) << std::endl;
    std::cout << "vertexs.at(2)" << vertexs.at(2) << std::endl;
  }
//     throw std::invalid_argument("projection: determinant is NaN.");
 
  alpha = (para_b[1]*normal[2] - normal[1]*para_b[2])*_pos[0] + (para_b[2]*normal[0] - para_b[0]*normal[2])*_pos[1] + (para_b[0]*normal[1] - para_b[1]*normal[0])*_pos[2];
   beta = (para_a[2]*normal[1] - normal[2]*para_a[1])*_pos[0] + (para_a[0]*normal[2] - para_a[2]*normal[0])*_pos[1] + (para_a[1]*normal[0] - para_a[0]*normal[1])*_pos[2];
  gamma = (para_a[1]*para_b[2] - para_a[2]*para_b[1])*_pos[0] + (para_a[2]*para_b[0] - para_a[0]*para_b[2])*_pos[1] + (para_a[0]*para_b[1] - para_a[1]*para_b[0])*_pos[2];

  alpha = alpha / determinant;
   beta = beta / determinant;
  gamma = gamma / determinant;

  if(gamma > 1e-4)
    throw std::invalid_argument("projection: gamma > 1e-4");

  return std::make_pair(alpha, beta);
}

Realvec FaceClose::to_absolute( const std::pair<Real, Real>& parameters )
{
  Realvec absolute_pos( para_a * parameters.first + para_b * parameters.second );
  return absolute_pos;
}

/*************************************** renew position end ***********************************/

void FaceClose::set_near_vertexs()
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


void FaceClose::set_neighbors_edge()
{
  for(int i(0); i<3; ++i)
  {
    FaceBase_sptr ptr(poly_ptr.lock()->get_neighbor(face_id, i));

    Realvec neighbor_para_a( ptr->get_para_a() );
    Realvec neighbor_para_b( ptr->get_para_b() );

    Realvec neighbor_normal( ptr->get_normal_vector() );

    Real rot_angle;
    if(fabs(dot_product(normal, neighbor_normal) <= 1e0) )
    {
      rot_angle = (-1e0) * acos( dot_product(normal, neighbor_normal) );
    }
    else if(dot_product(normal, neighbor_normal) > 0e0)
    {
      std::cout << "Warning: dot product of neighbor normal vectors is : 1 + ";
      std::cout << std::scientific << std::setprecision(16) << (dot_product(normal, neighbor_normal) - 1e0) << std::endl;
      std::cout << "This treats as just 1. acos(dot product) is 0e0." << std::endl;
      std::cout << "face id is " << get_id() << std::endl;
      std::cout << std::endl;
      rot_angle = 0e0;
    }else{
      std::cout << "Warning: dot product of neighbor normal vectors is : -1 + ";
      std::cout << std::scientific << std::setprecision(16) << (dot_product(normal, neighbor_normal) + 1e0) << std::endl;;
      std::cout << "This treats as just -1. acos(dot product) is M_PI." << std::endl;
      std::cout << "face id is " << get_id() << std::endl;
      std::cout << std::endl;
      rot_angle = M_PI;
    }

    Realvec axis( cross_product(normal, neighbor_normal) );
    if( length(axis) == 0e0 || rot_angle == 0e0 )
    {
      ;//do nothing
    }
    else
    {
      neighbor_para_a = rotation( rot_angle, axis/length(axis), neighbor_para_a );
      neighbor_para_b = rotation( rot_angle, axis/length(axis), neighbor_para_b );
    }

//     if( fabs( dot_product( neighbor_para_a, normal ) ) > GLOBAL_TOLERANCE )
//     {
//       std::cout << "dot product: " << fabs(dot_product( neighbor_para_a, normal )) << std::endl;
//       throw std::invalid_argument("rotated neighbor edge but it is not on this face");
//     }
//
//     if( fabs( dot_product( neighbor_para_b, normal ) ) > GLOBAL_TOLERANCE )
//     {
//       std::cout << "dot product: " << fabs(dot_product( neighbor_para_b, normal )) << std::endl;
//       throw std::invalid_argument("rotated neighbor edge but it is not on this face");
//     }

    std::pair<Real, Real> neighbor_a( projection( neighbor_para_a ) );
    std::pair<Real, Real> neighbor_b( projection( neighbor_para_b ) );

    // ( neighbor_a.first neighbor_a.second ) (para_a) = (neighbor_para_a)
    // ( neighbor_b.first neighbor_b.second ) (para_b)   (neighbor_para_b)
    Real determinant(neighbor_a.first * neighbor_b.second - neighbor_a.second * neighbor_b.first);
    if(determinant == 0e0)
      throw std::invalid_argument("set_neighbors_edge: determinant is zero");
    if(isnan( determinant ))
    {
      std::cout << "face id: " << get_id() << std::endl;
      std::cout << "set_neighbors_edge: determinant is NaN" << std::endl;
      std::cout << "neighbor_a: (" << neighbor_a.first << ", " << neighbor_a.second << ")" << std::endl;
      std::cout << "neighbor_b: (" << neighbor_b.first << ", " << neighbor_b.second << ")" << std::endl;
      std::cout << "rotated neighbor_para_a: " << neighbor_para_a << std::endl;
      std::cout << "rotated neighbor_para_b: " << neighbor_para_b << std::endl;
      std::cout << "rotation angle: " << rot_angle << std::endl;
      std::cout << "dot_product(normal, neighbor_normal): " << std::setprecision(16) << dot_product(normal, neighbor_normal) << std::endl;
      std::cout << "length(axis): " << length(axis) << std::endl;
      std::cout << "neighbor_para_a: " << ptr->get_para_a() << std::endl;
      std::cout << "neighbor_para_b: " << ptr->get_para_b() << std::endl;

    }
//       throw std::invalid_argument("set_neighbors_edge: determinant is NaN");

    // ( neighbor_b.second -neighbor_a_second)
    // ( -neighbor_b.first neighbor_a.first )
    std::pair<Real, Real> nbr_a( neighbor_b.second / determinant, -neighbor_a.second / determinant );
    std::pair<Real, Real> nbr_b( -neighbor_b.first / determinant, neighbor_a.first / determinant );

    if( ( isnan(nbr_a.first) || isnan(nbr_a.second) ) || (isnan(nbr_b.first) || isnan(nbr_b.second)) )
    {
      std::cout << "nbr_a: " << nbr_a.first <<", "<< nbr_a.second << std::endl;
      std::cout << "nbr_b: " << nbr_b.first <<", "<< nbr_b.second << std::endl;
      std::cout << "neighbor_a.first" << neighbor_a.first << ", neighbor_a.sec" << neighbor_a.second << std::endl;
      std::cout << "neighbor_b.first" << neighbor_b.first << ", neighbor_b.sec" << neighbor_b.second << std::endl;
      std::cout << "determinant: " << determinant << std::endl;
    }

    para_a_neighbor.at(i) = nbr_a;
    para_b_neighbor.at(i) = nbr_b;
  }
  return;
}

void FaceClose::set_neighbors_ori_vtx()
{
  for(int i(0); i<3; ++i)
  {
    FaceBase_sptr ptr(poly_ptr.lock()->get_neighbor(face_id, i));
    Realvec vi_nbr_ori( ptr->get_para_origin() - vertexs.at(i) );

    Realvec neighbor_normal( ptr->get_normal_vector() );

    Real rot_angle;
    if(fabs(dot_product(normal, neighbor_normal) <= 1e0) )
    {
      rot_angle = (-1e0) * acos( dot_product(normal, neighbor_normal) );
    }else if(dot_product(normal, neighbor_normal) > 0e0){
      std::cout << "Warning: dot product of neighbor normal vectors is : 1 + ";
      std::cout << std::scientific << std::setprecision(16);
      std::cout << (dot_product(normal, neighbor_normal) - 1e0) << std::endl;
      std::cout << "This treats as just 1. acos(dot product) is 0e0." << std::endl;
      std::cout << "face id is " << get_id() << std::endl;
      std::cout << std::endl;
      rot_angle = 0e0;
    }else{
      std::cout << "Warning: dot product of neighbor normal vectors is : -1 + ";
      std::cout << std::scientific << std::setprecision(16);
      std::cout << (dot_product(normal, neighbor_normal) + 1e0) << std::endl;
      std::cout << "This treats as just -1. acos(dot product) is M_PI." << std::endl;
      std::cout << "face id is " << get_id() << std::endl;
      std::cout << std::endl;
      rot_angle = M_PI;
    }


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
    std::pair<Real, Real> neighbor_ori_this_ori( projection( nbrori_ori ) );

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
Realvec FaceClose::get_another_vertex(const Realvec& neighbors_edge)
{
  return vertexs.at( get_another_vertex_id(neighbors_edge) );
}

int FaceClose::get_another_vertex_id(const Realvec& neighbors_edge)
{
  int i(matching_edge(neighbors_edge));
  return (i+2) % 3;
}

int FaceClose::get_another_edge_id_right(const Realvec& neighbors_edge)
{
  int i(matching_edge(neighbors_edge) );
  return ( (i+1) % 3 );
}

int FaceClose::get_another_edge_id_left(const Realvec& neighbors_edge)
{
  int i(matching_edge(neighbors_edge) );
  return ((i+2) % 3);
}

Real FaceClose::get_max_a(const Realvec& position, bool& vertex_include_flag)
{
  vertex_include_flag = false;

  Real vertdist( length( vertexs.at(0) - position ) );
  Real min_distance( vertdist );

  for(int i(1); i<3; ++i)
  {
    vertdist = length( vertexs.at(i) - position );
    if( min_distance > vertdist ) min_distance = vertdist;
  }

//****** near triangles
  if( !near_vert_height.empty() )
  {
    int size( near_vert_height.size() );

    for(int i(0); i<size; ++i)
    {
      if( min_distance > near_vert_height.at(i) )
          min_distance = near_vert_height.at(i);
    }
  }

  // if minimal distance from particle to vertex is smaller than this threshold
  // this allows the shell to include only one vertex.
  if(min_distance < VERTEX_THRESHOLD)
  {
    vertex_include_flag = true;
    Real second_distance(-1e0);

    for(int i(0); i < 3; ++i)
    {
      Real len(length(vertexs.at(i) - position));

      if( fabs(min_distance - len) < GLOBAL_TOLERANCE )
      {
        continue;
      }
      else if(second_distance > len || second_distance < 0e0)
      {
        //len is not minimum and shorter than second_distance
        //.OR.
        //second_distance is not substituted yet then
        second_distance = len;
      }
    }

    if(near_vert_height.empty())
      return second_distance;

    int size(near_vert_height.size());
    for(int i(0); i < size; ++i)
    {
      if(fabs(min_distance - near_vert_height.at(i)) < GLOBAL_TOLERANCE)
        continue;
      if(second_distance > near_vert_height.at(i) )
        second_distance = near_vert_height.at(i);
    }

    THROW_UNLESS(std::invalid_argument, second_distance > VERTEX_THRESHOLD);

    return second_distance;
  }else{
    return min_distance;
  }
}

Real FaceClose::get_left_angle( const Realvec& neighbors_edge )
{
  int id( matching_edge(neighbors_edge) );
  return angles.at(id);
}

Real FaceClose::get_right_angle( const Realvec& neighbors_edge )
{
  int id( matching_edge(neighbors_edge) );
  return angles.at( (id+1) % 3);
}

int FaceClose::matching_edge(const Realvec& neighbors_edge )
{
  Realvec inverse( neighbors_edge * (-1e0) );
  //assert when the edge shared by two faces, direction of the edge is opposit between these faces.
  for(int i(0); i<3; ++i)
  {
    if( fabs(edges.at(i)[0] - inverse[0]) < GLOBAL_TOLERANCE &&
        fabs(edges.at(i)[1] - inverse[1]) < GLOBAL_TOLERANCE &&
        fabs(edges.at(i)[2] - inverse[2]) < GLOBAL_TOLERANCE )
    {
      return i;
    }
  }
  throw std::invalid_argument( "matching_edge could not find edge that matches argument." );
}


Real FaceClose::get_minimum_height( const Realvec& neighbors_edge )
{
  const Real perpendicular(M_PI * 0.5);
  int edge_id( matching_edge(neighbors_edge) );

  //TODO: depends on the directions of normal vectors of each faces... <- true?
  Real left_angle( this->angles.at(edge_id) );
  Real right_angle( this->angles.at( (edge_id+1)%3 ) );

  // parallelogram's area( 2 times larger than that of this face )
  Real area(length(cross_product(edges.at(0), edges.at(2)*(-1e0) ) ) );
  Real min_height(area / length(neighbors_edge) );

  if(min_height <= 0e0) throw std::invalid_argument("min_height is negative or zero");

//****************************************************************************//
  FaceBase_sptr right_face( poly_ptr.lock()->get_neighbor(face_id, (edge_id+1)%3 ) );
  FaceBase_sptr  left_face( poly_ptr.lock()->get_neighbor(face_id, (edge_id+2)%3 ) );
  Realvec right_edge(edges.at( (edge_id + 1)%3 ) );
  Realvec  left_edge(edges.at( (edge_id + 2)%3 ) );
  Real r( (length(neighbors_edge) / 2e0) );

  while(true)
  {
    right_angle += right_face->get_right_angle(right_edge);
    left_angle  +=  left_face->get_left_angle( left_edge );
    
    if( left_angle >= perpendicular && right_angle >= perpendicular )
      break;

    if(left_angle < perpendicular)
    {
      int left_edge_id( (left_face->matching_edge(left_edge) + 2) % 3 );
      left_edge = left_face->get_edges_at(left_edge_id);

      Real max_l(r * sqrt(2e0 * (1e0 + cos(M_PI - 2 * left_angle) ) ) );
      Real l(length(left_edge) );
      if(l < max_l)
      {
        Real height( l * sin(left_angle) );
        if(min_height > height)
          min_height = height;
      }
      left_face = poly_ptr.lock()->get_neighbor( left_face->get_id(), left_edge_id );
    }

    if(right_angle < perpendicular )
    {
      int right_edge_id( (right_face->matching_edge(right_edge) + 1) % 3 );
      right_edge = right_face->get_edges_at(right_edge_id);

      Real max_l(r * sqrt(2e0 * (1e0 + cos(M_PI - 2 * right_angle) ) ) );
      Real l(length(right_edge) );
      if(l < max_l)
      {
        Real height( l * sin(right_angle) );
        if(min_height > height)
          min_height = height;
      }
      right_face = poly_ptr.lock()->get_neighbor( right_face->get_id(), right_edge_id );
    }
  }

  return min_height;
}

void FaceClose::print_class_name()
{
  std::cout << "class: FaceClose" << std::endl;
}

#endif /*FACE_CLOSE*/
