#ifndef FACE_TWO_GATE
#define FACE_TWO_GATE
#include <vector>
#include <cmath>
#include "Polygon.hpp"
#include "rotation.hpp"

using namespace greens_functions;

class FaceTwoGate : public FaceBase
{
private:
  std::vector<Realvec> vertexs;
  //to get the vertex nearest the particle, lookup in (std::vector vertexs) and this near_vertexs.
  std::vector<Realvec> near_vertexs;
  std::vector<Realvec> edges;
  std::vector<int> gateway;// edges.at(gateway.at(i)) == gateway edge
  std::vector<bool> is_gate;
  boost::shared_ptr<Polygon> belonging_polygon;
  enum CROSS_STATUS{
    NOT_CROSS,
    CROSS,
    ON_THE_EDGE,
  };

public:
  FaceTwoGate(const int& id, const Realvec& vtx0, const Realvec& vtx1, const Realvec& vtx2, const int& gate0, const int& gate1)
  : FaceBase( id, cross_product(vtx1-vtx0, vtx2-vtx0), vtx1-vtx0 ), gateway(2), is_gate(3)
  {
    vertexs.push_back( vtx0 );
    vertexs.push_back( vtx1 );
    vertexs.push_back( vtx2 );

    edges.push_back( vtx1 - vtx0 );//edges.at(i): vertexs.at(i) -> vertexs.at(i+1)
    edges.push_back( vtx2 - vtx1 );
    edges.push_back( vtx0 - vtx2 );
  
    gateway.at(0) = gate0;
    gateway.at(1) = gate1;

    for(int i(0); i<3; ++i)
    {
      if(i != gate0 && i != gate1)
      {
	is_gate.at(i) = false;
      }else
      {
	is_gate.at(i) = true;
      }
    }
  }

  FaceTwoGate(const int& id, const Realvec& vtx0, const Realvec& vtx1, const Realvec& vtx2, const Realvec& norm, const int& gate0, const int& gate1)
  : FaceBase( id, norm, vtx1-vtx0 ), gateway(2)
  {
    bool normal_vec_is_oritented_orthogonally_to_the_edges(
 		 dot_product(norm, vtx1-vtx0) == 0 &&
		 dot_product(norm, vtx2-vtx1) == 0 );
    THROW_UNLESS(std::invalid_argument, normal_vec_is_oritented_orthogonally_to_the_edges );
  
    vertexs.push_back( vtx0 );
    vertexs.push_back( vtx1 );
    vertexs.push_back( vtx2 );

    edges.push_back( vtx1 - vtx0 );
    edges.push_back( vtx2 - vtx1 );
    edges.push_back( vtx0 - vtx2 );

    gateway.at(0) = gate0;
    gateway.at(1) = gate1;

    for(int i(0); i<3; ++i)
    {
      if(i != gate0 && i != gate1)
      {
	is_gate.at(i) = false;
      }else
      {
	is_gate.at(i) = true;
      }
    }
  }

  virtual void set_belonging_polygon( boost::shared_ptr<Polygon> p_sptr)
  {
    belonging_polygon = p_sptr;
  };

  virtual Realvec renew_position( const Realvec& position, const Realvec& displacement, boost::shared_ptr<FaceBase>& ptr);

  virtual bool still_in_the_face( const Realvec& position, const Realvec& displacement );

  virtual void set_near_vertexs();
  virtual Realvec get_another_vertex(const Realvec& edge);
  virtual Realvec get_vertex(){ return vertexs.at(0); };

private:

  CROSS_STATUS is_cross( const Realvec& position, const Realvec& displacement, const Realvec& edge, const int& edgebase );

  bool is_on_the_edge( const Realvec& position, const Realvec& edge, const int& edgebase );

  Real cross_point( const Realvec& position, const Realvec& displacement, const Realvec& edge, const int& edgebase );
  
  Realvec reverse( const Realvec& target, const Realvec& edge );

//   int get_gateway( const int& num)
//   {
//     return gateway.at(num);
//   };

//   Realvec get_gateway_edge(const int& g)
//   {
//     return edges.at(g);
//   };

};



Realvec FaceTwoGate::renew_position(const Realvec& position, const Realvec& displacement, boost::shared_ptr<FaceBase>& ptr)
{
  Realvec temppos(position);
  Realvec tempdis(displacement);
  //debug
  std::cout << temppos[0] << " " << temppos[1] << " " << temppos[2] << " ";
  std::cout << tempdis[0] << " " << tempdis[1] << " " << tempdis[2] << std::endl;
  
  bool in_the_face( still_in_the_face(temppos, tempdis) );
  if(in_the_face)
  {
    return temppos + tempdis;
  }
  else
  { //particle will go out of the face
    int size( edges.size() );
    do
    {
      int on_this_edge(-1);

      for(int i(0); i < size; ++i)
      {
	if( !is_on_the_edge( temppos, edges.at(i), i ) ) continue;

	if( !is_gate.at(i) )
	{
	  on_this_edge = i;

	  int num( (i+2) % 3 );
	  
	  Realvec cross1( cross_product( edges.at(i), ( edges.at(num) * (-1e0) ) ) );
	  Realvec cross2( cross_product( edges.at(i), ( temppos + tempdis - vertexs.at(i) ) ) );
	  Real cross_status( dot_product( cross1, cross2 ) );

	  if(cross_status < 0e0)
	  {
	    Realvec reverse_disp( reverse( tempdis, edges.at(i) ) );
	    tempdis = reverse_disp;
	  }else{
	    ;//do nothing
	  }

	  break;

	}else{

	  int gatenum;
	  if(i == gateway.at(0))
	  {
	    gatenum = 0;
	  }else{
	    THROW_UNLESS(std::invalid_argument, i==gateway.at(1));
	    gatenum = 1;
	  }

	  int num( (i+2) % 3 );

	  Realvec cross1( cross_product( edges.at(i), ( edges.at(num) * (-1e0) ) ) );
	  Realvec cross2( cross_product( edges.at(i), ( temppos + tempdis - vertexs.at(i) ) ) );
	  Real cross_status( dot_product( cross1, cross2 ) );

	  // position + displacement is out of the face through edges.at(i)
	  // and the edge is gateway to another face
	  if(cross_status < 0e0)
	  {
	    ptr = belonging_polygon->get_neighbor_ptr_from_gateway(face_id, gateway.at(gatenum));
	    Realvec neighbor_norm( ptr->get_normal_vector() );
	    Real theta( acos( dot_product(normal, neighbor_norm) ) );

	    Realvec rot_disp( rotation( theta, edges.at(gateway.at(gatenum)), tempdis) );
	    tempdis = rot_disp;
	    temppos = ptr->renew_position( temppos, tempdis, ptr );

	    return temppos;

	  }else{
	    on_this_edge = i;
	  }

	  break;
	}
      }

      for(int i(0); i<size; ++i)
      {
	if( !is_gate.at(i) )
	{
	  if( i == on_this_edge ) continue;

	  CROSS_STATUS cross_status( is_cross(temppos, tempdis, edges.at(i), i) );
	  if(cross_status == CROSS)
	  {
	    Real ratio( cross_point(temppos, tempdis, edges.at(i), i) );
	    THROW_UNLESS(std::invalid_argument, 0e0 <= ratio && ratio <= 1e0);

	    temppos = temppos + tempdis * ratio;
	    tempdis = tempdis * (1e0 - ratio);
	    Realvec reverce_disp( reverse( tempdis, edges.at(i) ) );
	    tempdis = reverce_disp;

	    break;
	  }else if(cross_status == ON_THE_EDGE)
	  {
	    temppos = temppos + tempdis;
	    Realvec zero(0e0, 0e0, 0e0);
	    tempdis = zero;

	    break;
	  }else if(cross_status == NOT_CROSS)
	  {
// 	    std::cout << "NOT_CROSS: non-gateway edge " << i << std::endl;//do nothing
	  }else{
	    throw std::invalid_argument("invalid cross_status");
	  }

	}else{//the case of gateway edge
	  if(i == on_this_edge) continue;
	  //when particle is on the edge and go out of the face, it has been treated above
	  //so when particle is on the edge, displacement directs inner region of the face

	  int gatenum;
	  if(i == gateway.at(0))
	  {
	    gatenum = 0;
	  }else{
	    THROW_UNLESS(std::invalid_argument, i==gateway.at(1));
	    gatenum = 1;
	  }

	  CROSS_STATUS cross_status( is_cross(temppos, tempdis, edges.at(i), i) );
	  //std::cout << cross_status << std::endl; 
	  //TODO particle is just on the edge, NOT_CROSS and enter infinit loop
	  if(cross_status == CROSS)
	  {
	    ptr = belonging_polygon->get_neighbor_ptr_from_gateway(face_id, gateway.at(gatenum));
	    Realvec neighbor_norm( ptr->get_normal_vector() );
	    Real theta( acos( dot_product(normal, neighbor_norm) ) );

	    Real ratio( cross_point(temppos, tempdis, edges.at(i), i) );
	    THROW_UNLESS(std::invalid_argument, 0e0 <= ratio && ratio <= 1e0);

	    //put particle on the edge
	    temppos = temppos + tempdis * ratio;
	    tempdis = tempdis * (1e0 - ratio);
  
	    Realvec rot_disp( rotation( theta, edges.at(gateway.at(gatenum)), tempdis ) );
	    tempdis = rot_disp;
	    temppos = ptr->renew_position( temppos, tempdis, ptr );
	    
	    return temppos;
	  }else if(cross_status == ON_THE_EDGE)
	  {
	    temppos = temppos + tempdis;
	    Realvec zero(0e0, 0e0, 0e0);
	    tempdis = zero;

	    break;
	  }else if( cross_status == NOT_CROSS )
	  {
//  	    std::cout << "NOT_CROSS: gateway edge " << i << std::endl;
// 	    std::cout << temppos << " " << tempdis << " " << temppos + tempdis << std::endl;
//  	    std::cout << still_in_the_face(temppos, tempdis) << std::endl;
// 	    std::cout << i << std::endl;
// 	    std::cout << temppos[0] << " " << temppos[1] << " " << temppos[2] << " ";
// 	    std::cout << tempdis[0] << " " << tempdis[1] << " " << tempdis[2] << " " << "NOT_CROSS" << std::endl;
	  }
	}
      }

    }while( !still_in_the_face(temppos, tempdis) );
  }
  return temppos + tempdis; 
};

bool FaceTwoGate::still_in_the_face( const Realvec& position, const Realvec& displacement )
{
  //caution: sometimes rotation andle is nearly equal and smaller than 2pi
  //         (typically in order of 1e-15 ), in that case this func returns 0
  //         to treat this, use epsilon
  //
  //         when distance between position + displacement and edge < 1e-5,
  //         2pi - angle can be laeger than 1e-12 even pos+dis is in the face...
  for(int i(0); i<3; ++i)
  {
    if( is_on_the_edge( position+displacement, edges.at(i), i) ) return true;
  }

  Real epsilon(1e-8);
  Realvec temppos( position + displacement );
  Realvec v0p( vertexs.at(0) - temppos );
  Realvec v1p( vertexs.at(1) - temppos );
  Realvec v2p( vertexs.at(2) - temppos );
  Real lv0p( length(v0p) );
  Real lv1p( length(v1p) );
  Real lv2p( length(v2p) );

  Real angle;
  angle = acos( dot_product( v0p, v1p ) / lv0p / lv1p );
  angle += acos( dot_product( v1p, v2p ) / lv1p / lv2p );
  angle += acos( dot_product( v2p, v0p ) / lv2p / lv0p );

  return (angle > 2 * M_PI - epsilon);
};

FaceTwoGate::CROSS_STATUS FaceTwoGate::is_cross
  ( const Realvec& position, const Realvec& displacement, const Realvec& edge, const int& edgebase )
{
  //assuming that position vector is on the face
  //and displacement is also on the same plane

  //edge line and displacement segment
  Realvec vp( position - vertexs.at(edgebase));
  Realvec vpd( position + displacement - vertexs.at(edgebase));
  Realvec cross1( cross_product(edge, vp) );
  Realvec cross2( cross_product(edge, vpd) );
  Real cross_stat1( dot_product(cross1, cross2) );

//displacement line and edge segment
  Realvec pv( vertexs.at(edgebase) - position );
  Realvec pve( vertexs.at(edgebase) + edge - position );
  Realvec cross3( cross_product(displacement, pv) );
  Realvec cross4( cross_product(displacement, pve) );
  Real cross_stat2( dot_product(cross3, cross4) );

  if( cross_stat1 < 0 && cross_stat2 < 0 )
  {
    return CROSS;
  }else if( cross_stat1 > 0 || cross_stat2 > 0 ){
    return NOT_CROSS;
  }else if( cross_stat1 == 0 && cross_stat2 < 0){
    return ON_THE_EDGE;
  }else{
    std::cout << "caution: through the vertex" << std::endl;
    return ON_THE_EDGE;
  }
};

bool FaceTwoGate::is_on_the_edge( const Realvec& position, const Realvec& edge, const int& edgebase )
{
  Real epsilon(1e-12);
  Realvec vp( position - vertexs.at(edgebase) );
  Real cross( length( cross_product(vp, edge) ) );

  return ( cross < epsilon );
};

Real FaceTwoGate::cross_point( const Realvec& position, const Realvec& displacement, const Realvec& edge, const int& edgebase )
{
//displacement and the edge must cross each other
//use is_cross() for assure
  Realvec v1p( position - vertexs.at(edgebase));
  Realvec v1pd( position + displacement - vertexs.at(edgebase));
  Real d1( length( cross_product(edge, v1p) ) / length(edge) );
  Real d2( length( cross_product(edge, v1pd)) / length(edge) );
  Real ratio( d1 / (d1 + d2) );

  return ratio;
}

Realvec FaceTwoGate::reverse( const Realvec& target, const Realvec& edge )
{
  Realvec edge_n( edge / length( edge ) );

  Real d( 2e0 * dot_product( target, edge_n ) );
  Realvec retvec( edge_n * d - target );

  return retvec;
};

void FaceTwoGate::set_near_vertexs()
{
  for(int i(0); i<3; ++i)
  {
    if( !is_gate.at(i) ) continue;

    Realvec midpoint( ( vertexs.at(i) + vertexs.at( (i+1)%3 ) ) / 2e0 );
    Real r( length( edges.at(i) ) / 2e0 );

    //edge(i) is gateway so i is gateway-id
    boost::shared_ptr<FaceBase> ptr( belonging_polygon->get_neighbor_ptr_from_gateway(face_id, i) );
    Realvec near_vertex( ptr->get_another_vertex(edges.at(i)) );
    
    Real l( length( near_vertex - midpoint ) );

    if(l <= r)
    {
      near_vertexs.push_back( near_vertex );
    }

  }
}

Realvec FaceTwoGate::get_another_vertex(const Realvec& edge)
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



#endif /*FACE_TWO_GATE*/
