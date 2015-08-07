#ifndef FACE_ONE_GATE
#define FACE_ONE_GATE
#include <vector>
#include <cmath>
#include "Polygon.hpp"
#include "rotation.hpp"

using namespace greens_functions;

class FaceOneGate : public FaceBase
{
private:
  std::vector<Realvec> vertexs;
  //to get the vertex nearest the particle, lookup in (std::vector vertexs) and this near_vertexs.
  std::vector<Realvec> near_vertexs;
  std::vector<Realvec> edges;
  int gateway;// edges.at(gateway) == gateway edge
  std::vector<bool> is_gate;
  boost::shared_ptr<Polygon> belonging_polygon;

  enum CROSS_STATUS{
    NOT_CROSS,
    CROSS,
    ON_THE_EDGE,
  };


public:
  FaceOneGate(const int& id, const Realvec& vtx0, const Realvec& vtx1, const Realvec& vtx2, const int& gate)
  : FaceBase( id, cross_product(vtx1-vtx0, vtx2-vtx0), vtx1-vtx0 ), gateway(gate), is_gate(3)
  {
    vertexs.push_back( vtx0 );
    vertexs.push_back( vtx1 );
    vertexs.push_back( vtx2 );

    edges.push_back( vtx1 - vtx0 );
    edges.push_back( vtx2 - vtx1 );
    edges.push_back( vtx0 - vtx2 );

    for(int i(0); i<3; ++i)
    {
      if(i != gate)
      {
	is_gate.at(i) = false;
      }else
      {
	is_gate.at(i) = true;
      }
    }

  }

  FaceOneGate(const int& id, const Realvec& vtx0, const Realvec& vtx1, const Realvec& vtx2, const Realvec& norm, const int& gate)
  : FaceBase( id, norm, vtx1-vtx0 ), gateway(gate), is_gate(3)
  {
    bool normal_vec_is_oritented_orthogonally_to_the_edges( dot_product(norm, vtx1-vtx0) == 0 && dot_product(norm, vtx2-vtx1) == 0 );
    THROW_UNLESS(std::invalid_argument, normal_vec_is_oritented_orthogonally_to_the_edges );
  
    vertexs.push_back( vtx0 );
    vertexs.push_back( vtx1 );
    vertexs.push_back( vtx2 );

    edges.push_back( vtx1 - vtx0 );
    edges.push_back( vtx2 - vtx1 );
    edges.push_back( vtx0 - vtx2 );

    for(int i(0); i<3; ++i)
    {
      if(i != gate)
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

  virtual void set_near_vertexs();

  virtual Realvec renew_position( const Realvec& position, const Realvec& displacement, boost::shared_ptr<FaceBase>& ptr);

  virtual bool still_in_the_face( const Realvec& position, const Realvec& displacement );

  virtual Realvec get_vertex()
  {
    return vertexs.at(0); 
  };

  virtual Realvec get_another_vertex(const Realvec& edge);
  virtual void print_class_name();

//   bool is_gateway_edge( int i )
//   {
//     return is_gate.at(i);
//   }

private:

  CROSS_STATUS is_cross( const Realvec& position, const Realvec& displacement, const Realvec& edge, const int& edgebase );

  bool is_on_the_edge( const Realvec& position, const Realvec& edge, const int& edgebase );

  Real cross_point( const Realvec& position, const Realvec& displacement, const Realvec& edge, const int& edgebase );
  
  Realvec reverse( const Realvec& target, const Realvec& edge );

  Realvec renew_position_on_nongate_edge(const int i, const Realvec& pos, const Realvec& dis);//returns new displacement

  Realvec renew_position_goes_another_face_from_on_edge(const Realvec& pos, const Realvec& dis, boost::shared_ptr<FaceBase>& ptr);//returns new position and renew the FaceBase ptr

  int get_gateway()
  {
    return gateway;
  };

};


Realvec FaceOneGate::renew_position
  (const Realvec& position, const Realvec& displacement, boost::shared_ptr<FaceBase>& ptr)
{
  Realvec temppos(position);
  Realvec tempdis(displacement);
  std::cout << temppos[0] << " " << temppos[1] << " " << temppos[2] << " ";
  std::cout << tempdis[0] << " " << tempdis[1] << " " << tempdis[2] << std::endl;
  
  bool in_the_face( still_in_the_face(temppos, tempdis) );
  if(in_the_face)
  {
    return temppos + tempdis;
  }else
  {
    do
    {
      int size( edges.size() );
      int on_this_edge(-1);

      //particle is on the edge
      for(int i(0); i < size; ++i)
      {
	if( !is_on_the_edge( temppos, edges.at(i), i ) ) continue;

	if( !is_gate.at(i) )
	{
	  on_this_edge = i;
	  tempdis = renew_position_on_nongate_edge(i, temppos, tempdis);
	  break;
	}
	else
	{
	  int num( (i+2) % 3 );

	  Realvec cross1( cross_product( edges.at(i), ( edges.at(num) * (-1e0) ) ) );
	  Realvec cross2( cross_product( edges.at(i), ( temppos + tempdis - vertexs.at(i) ) ) );
	  Real cross_status( dot_product( cross1, cross2 ) );

	  if(cross_status < 0e0)
	  {
	    //out of the face ( renew position and ptr )
	    temppos = renew_position_goes_another_face_from_on_edge( temppos, tempdis, ptr );

	    return temppos;

	  }else{
	    //still in the face
	    on_this_edge = i;
	  }

	  break;
	}
      }

      for(int i(0); i<size; ++i)
      {
	if( !is_gate.at(i) )
	{
	  if(i == on_this_edge) continue;
	  CROSS_STATUS cross_status( is_cross(temppos, tempdis, edges.at(i), i) );
	  if(cross_status == CROSS)
	  {
	    Real ratio( cross_point(temppos, tempdis, edges.at(i), i) );

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
	  }
	}else{
	  if(i == on_this_edge) continue;

	  CROSS_STATUS cross_status( is_cross(temppos, tempdis, edges.at(i), i) );
	  if(cross_status == CROSS)
	  {
	    ptr = belonging_polygon->get_neighbor_ptr_from_gateway(face_id, gateway);
	    Realvec neighbor_norm( ptr->get_normal_vector() );
	    Real theta( acos( dot_product(normal, neighbor_norm) ) );

	    Real ratio( cross_point(temppos, tempdis, edges.at(i), i) );

	    //put particle on the edge
	    temppos = temppos + tempdis * ratio;
	    tempdis = tempdis * (1e0 - ratio);

  	    Realvec norm1( cross_product(normal, neighbor_norm) );
	    if( length(norm1) == 0 ) norm1 = edges.at(gateway);
	    Realvec axis( norm1 / length(norm1) );

	    Realvec rot_disp( rotation( theta, axis, tempdis ) );
	    tempdis = rot_disp;
	    temppos = ptr->renew_position( temppos, tempdis, ptr );
	    
	    return temppos;
	  }else if(cross_status == ON_THE_EDGE)
	  {
	    temppos = temppos + tempdis;
	    Realvec zero(0e0, 0e0, 0e0);
	    tempdis = zero;

	    break;
	  }else
	  {
	    ;//do nothing: particle will go out of the face but doesnot pass through this edge.
	  }
	}
      }

    }while( !still_in_the_face(temppos, tempdis) );
  }
  return temppos + tempdis; 
}

Realvec FaceOneGate::renew_position_on_nongate_edge(const int i, const Realvec& pos, const Realvec& dis)
{
  Realvec temppos(pos), tempdis(dis);

  int num( (i + 2) % 3 );// (i-1) mod 3
  Realvec cross1( cross_product( edges.at(i), ( edges.at(num) * (-1e0) ) ) );
  Realvec cross2( cross_product( edges.at(i), ( temppos + tempdis - vertexs.at(i) ) ) );
  Real cross_status( dot_product( cross1, cross2 ) );

  if( cross_status < 0e0 )
  {
    //out of the face ( renew displacement )
    Realvec reverse_disp( reverse( tempdis, edges.at(i) ) );
    tempdis = reverse_disp;
    return tempdis;
  }else{
    //in the face ( do nothing. then renew position normally )
    return tempdis;
  }
}

Realvec FaceOneGate::renew_position_goes_another_face_from_on_edge
  (const Realvec& pos, const Realvec& dis, boost::shared_ptr<FaceBase>& ptr)
{
  Realvec temppos(pos), tempdis(dis);

  ptr = belonging_polygon->get_neighbor_ptr_from_gateway(face_id, gateway);
  Realvec neighbor_norm( ptr->get_normal_vector() );
  Real theta( acos( dot_product(normal, neighbor_norm) ) );

  //particle is on an edge and 
  //the neighbor face also have the same edge as a member
  Realvec norm1( cross_product(normal, neighbor_norm) );
  if( length(norm1) == 0 ) norm1 = edges.at(gateway);
  Realvec axis( norm1 / length(norm1) );
  
  Realvec rot_disp( rotation( theta, axis, tempdis) );
  Realvec retvec ( ptr->renew_position( temppos, rot_disp, ptr ) );

  return retvec;
}

bool FaceOneGate::still_in_the_face( const Realvec& position, const Realvec& displacement )
{
  //caution: sometimes rotation andle is nearly equal and smaller than 2pi
  //         (typically in order of 1e-15 ), in that case this func returns 0
  //         to treat this, use epsilon
  //
  //         if particle is on the edge, this func returns true.
  Real epsilon(1e-8);
  Realvec temppos(position + displacement);
  Realvec v0p( vertexs.at(0) - temppos );
  Realvec v1p( vertexs.at(1) - temppos );
  Realvec v2p( vertexs.at(2) - temppos );
  Real lv0p( length(v0p) );
  Real lv1p( length(v1p) );
  Real lv2p( length(v2p) );

  Real rot;
  rot = acos( dot_product( v0p, v1p ) / lv0p / lv1p );
  rot += acos( dot_product( v1p, v2p ) / lv1p / lv2p );
  rot += acos( dot_product( v2p, v0p ) / lv2p / lv0p );

  return (rot >= 2 * M_PI - epsilon);
}

FaceOneGate::CROSS_STATUS
FaceOneGate::is_cross( const Realvec& position, const Realvec& displacement, const Realvec& edge, const int& edgebase )
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
}

bool FaceOneGate::is_on_the_edge( const Realvec& position, const Realvec& edge, const int& edgebase )
{
  Real epsilon(1e-14);
  Realvec vp( position - vertexs.at(edgebase) );
  Real cross( length( cross_product(vp, edge) ) );

  return ( cross < epsilon );
}

Real FaceOneGate::cross_point( const Realvec& position, const Realvec& displacement, const Realvec& edge, const int& edgebase )
{
//displacement and the edge must cross each other
//use is_cross() for assure
  Realvec v1p( position - vertexs.at(edgebase));
  Realvec v1pd( position + displacement - vertexs.at(edgebase));
  Real d1( length( cross_product(edge, v1p) ) / length(edge) );
  Real d2( length( cross_product(edge, v1pd)) / length(edge) );
  Real ratio( d1 / (d1 + d2) );

  THROW_UNLESS(std::invalid_argument, 0e0 <= ratio && ratio <= 1e0);
  return ratio;
}

Realvec FaceOneGate::reverse( const Realvec& target, const Realvec& edge )
{
  Realvec edge_n( edge / length( edge ) );

  Real d( 2e0 * dot_product( target, edge_n ) );
  Realvec retvec( edge_n * d - target );

  return retvec;
}

void FaceOneGate::set_near_vertexs()
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
