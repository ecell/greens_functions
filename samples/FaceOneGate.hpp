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
  std::vector<Realvec> edges;
  int gateway;// edges.at(gateway) == gateway edge
  boost::shared_ptr<Polygon> belonging_polygon;

public:
  FaceOneGate(const int& id, const Realvec& vtx0, const Realvec& vtx1, const Realvec& vtx2, const int& gate)
  : FaceBase( id, cross_product(vtx1-vtx0, vtx2-vtx0), vtx1-vtx0 ), gateway(gate)
  {
    vertexs.push_back( vtx0 );
    vertexs.push_back( vtx1 );
    vertexs.push_back( vtx2 );

    edges.push_back( vtx1 - vtx0 );
    edges.push_back( vtx2 - vtx1 );
    edges.push_back( vtx0 - vtx2 );
  }

  FaceOneGate(const int& id, const Realvec& vtx0, const Realvec& vtx1, const Realvec& vtx2, const Realvec& norm, const int& gate)
  : FaceBase( id, norm, vtx1-vtx0 ), gateway(gate)
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
  }

  virtual void set_belonging_polygon( boost::shared_ptr<Polygon> p_sptr)
  {
    belonging_polygon = p_sptr;
  };

  virtual Realvec move(Realvec& position, Realvec& displacement, boost::shared_ptr<FaceBase>& ptr);

  virtual bool still_in_the_face( const Realvec& position, const Realvec& displacement );

  virtual Realvec get_vertex()
  {
    return vertexs.at(0); 
  };

private:

  int is_cross( const Realvec& position, const Realvec& displacement, const Realvec& edge, const int& edgebase );

  bool is_on_the_edge( const Realvec& position, const Realvec& edge, const int& edgebase );

  Real cross_point( const Realvec& position, const Realvec& displacement, const Realvec& edge, const int& edgebase );
  
  Realvec reverse( const Realvec& target, const Realvec& edge );

  int get_gateway()
  {
    return gateway;
  };

};



Realvec FaceOneGate::move(Realvec& position, Realvec& displacement, boost::shared_ptr<FaceBase>& ptr)
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

	if( i != gateway )
	{

	  on_this_edge = i;
	  int num;
	  if( i != 0)
	  {
	    num = i-1;
	  }else{
	    num = 2;
	  }

	  Realvec cross1( cross_product( edges.at(i), ( edges.at(num) * (-1e0) ) ) );
	  Realvec cross2( cross_product( edges.at(i), ( temppos + tempdis - vertexs.at(i) ) ) );
	  Real cross_status( dot_product( cross1, cross2 ) );
	  if(cross_status < 0e0)// disp is out of the face through edge(i)
	  {
	    Realvec reverse_disp( reverse( tempdis, edges.at(i) ) );
	    tempdis = reverse_disp;
	  }else{
	    ;//do nothing
	  }

	  break;
	}
	else
	{
	  int num;
	  if( i != 0)
	  {
	    num = i-1;
	  }else{
	    num = 2;
	  }

	  Realvec cross1( cross_product( edges.at(i), ( edges.at(num) * (-1e0) ) ) );
	  Realvec cross2( cross_product( edges.at(i), ( temppos + tempdis - vertexs.at(i) ) ) );
	  Real cross_status( dot_product( cross1, cross2 ) );
	  if(cross_status < 0e0)// disp is out of the face through edge(i)
	  {
	    ptr = belonging_polygon->get_neighbor_ptr_from_gateway(face_id, gateway);
	    Realvec neighbor_norm( ptr->get_normal_vector() );
	    Real theta( acos( dot_product(normal, neighbor_norm) ) );

	    Realvec rot_disp( rotation( theta, edges.at(gateway), tempdis) );
	    tempdis = rot_disp;
	    temppos = ptr->move( temppos, tempdis, ptr );

	    return temppos;

	  }else{
	    on_this_edge = i;
	  }

	  break;
	}
      }

      //particle is not on the edge.
      for(int i(0); i<size; ++i)
      {
	if( i != gateway )
	{
	  if(i == on_this_edge) continue;
	  int cross_status( is_cross(temppos, tempdis, edges.at(i), i) );
	  if(cross_status == 1)
	  {
	    Real ratio( cross_point(temppos, tempdis, edges.at(i), i) );
	    THROW_UNLESS(std::invalid_argument, 0e0 <= ratio && ratio <= 1e0);

	    temppos = temppos + tempdis * ratio;
	    tempdis = tempdis * (1e0 - ratio);
	    Realvec reverce_disp( reverse( tempdis, edges.at(i) ) );
	    tempdis = reverce_disp;

	    break;
	  }else if(cross_status == -1)
	  {
	    temppos = temppos + tempdis;
	    Realvec zero(0e0, 0e0, 0e0);
	    tempdis = zero;

	    break;
	  }
	}else{
	  if(i == on_this_edge) continue;
	  int cross_status( is_cross(temppos, tempdis, edges.at(i), i) );
	  if(cross_status == 1)
	  {
	    ptr = belonging_polygon->get_neighbor_ptr_from_gateway(face_id, gateway);
	    Realvec neighbor_norm( ptr->get_normal_vector() );
	    Real theta( acos( dot_product(normal, neighbor_norm) ) );

	    Real ratio( cross_point(temppos, tempdis, edges.at(i), i) );
	    THROW_UNLESS(std::invalid_argument, 0e0 <= ratio && ratio <= 1e0);

	    //move to edge
	    temppos = temppos + tempdis * ratio;
	    tempdis = tempdis * (1e0 - ratio);
  
	    Realvec rot_disp( rotation( theta, edges.at(gateway), tempdis ) );
	    tempdis = rot_disp;
	    temppos = ptr->move( temppos, tempdis, ptr );
	    
	    return temppos;
	  }else if(cross_status == -1)
	  {
	    temppos = temppos + tempdis;
	    Realvec zero(0e0, 0e0, 0e0);
	    tempdis = zero;

	    break;
	  }
	}
      }

    }while( !still_in_the_face(temppos, tempdis) );
  }
  return temppos + tempdis; 
};

bool FaceOneGate::still_in_the_face( const Realvec& position, const Realvec& displacement )
{
  //caution: sometimes rotation andle is nearly equal and smaller than 2pi
  //         (typically in order of 1e-15 ), in that case this func returns 0
  //         to treat this, use epsilon
  //
  //         if particle is on the edge, this func returns true.
  Real epsilon(1e-12);
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
};

// 0: not cross
// 1: cross
// -1 on the edge
int FaceOneGate::is_cross( const Realvec& position, const Realvec& displacement, const Realvec& edge, const int& edgebase )
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
    return 1;
  }else if( cross_stat1 > 0 || cross_stat2 > 0 ){
    return 0;
  }else if( cross_stat1 == 0 && cross_stat2 < 0){
    return -1;
  }else{
    std::cout << "caution: through the vertex" << std::endl;
    return -1;
  }
};

bool FaceOneGate::is_on_the_edge( const Realvec& position, const Realvec& edge, const int& edgebase )
{
  Real epsilon(1e-14);
  Realvec vp( position - vertexs.at(edgebase) );
  Real cross( length( cross_product(vp, edge) ) );

  return ( cross < epsilon );
};

Real FaceOneGate::cross_point( const Realvec& position, const Realvec& displacement, const Realvec& edge, const int& edgebase )
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

Realvec FaceOneGate::reverse( const Realvec& target, const Realvec& edge )
{
  Realvec edge_n( edge / length( edge ) );

  Real d( 2e0 * dot_product( target, edge_n ) );
  Realvec retvec( edge_n * d - target );

  return retvec;
};


#endif /*FACE_ONE_GATE*/
