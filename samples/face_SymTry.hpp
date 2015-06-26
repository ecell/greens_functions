#ifndef FACE_SYMTRY
#define FACE_SYMTRY
#include <vector>
#include <cmath>
#include "face_base.hpp"

using namespace greens_functions;

class face_SymTry : public face_base
{
  std::vector<Realvec> vertexs;
public:
  face_SymTry(const int& id, const Realvec& vtx0, 
	      const Realvec& vtx1, const Realvec& vtx2)
  : face_base( id, cross_product(vtx1-vtx0, vtx2-vtx0), vtx1-vtx0 )
  {
    vertexs.push_back( vtx0 );
    vertexs.push_back( vtx1 );
    vertexs.push_back( vtx2 );
 }
 
   face_SymTry(const int& id, const Realvec& vtx0, 
	       const Realvec& vtx1, const Realvec& vtx2,
	       const Realvec& norm)
  : face_base( id, norm, vtx1-vtx0 )
  {
    bool normal_vec_is_oritented_orthogonally_to_the_edges(
 		 dot_product(norm, vtx1-vtx0) == 0 &&
		 dot_product(norm, vtx2-vtx1) == 0 );
    THROW_UNLESS(std::invalid_argument, normal_vec_is_oritented_orthogonally_to_the_edges );

    vertexs.push_back( vtx0 );
    vertexs.push_back( vtx1 );
    vertexs.push_back( vtx2 );
 } 

  virtual Realvec move(Realvec& position,
		       Realvec& displacement,
		       boost::shared_ptr<face_base>& p)
  {
    bool in_the_face( still_in_the_face(position, displacement) );
    if(in_the_face)
    {
      Realvec temp(position + displacement);
      return temp;
    }else{
      do{
	int status_01( is_cross( position, displacement, 0, 1 ) );
	int status_12( is_cross( position, displacement, 1, 2 ) );
	int status_20( is_cross( position, displacement, 2, 0 ) );
	int pos_status( status_01 + 2*status_12 + 3*status_20 );

    //TODO case of particle is on the edge and on the vertex
	int v_num[2];
	switch(pos_status){
	  case 1:
	    v_num[0] = 0;
	    v_num[1] = 1;
	    break;
	  case 2:
	    v_num[0] = 1;
	    v_num[1] = 2;
	    break;
	  case 3:
	    v_num[0] = 2;
	    v_num[1] = 0;
	    break;
	  default:
	    bool position_status(false);
	    THROW_UNLESS(std::invalid_argument, position_status);
	    break;
	}
	Real ratio( cross_point(position, displacement, v_num[0], v_num[1]) );
	THROW_UNLESS(std::invalid_argument, 
		     0e0 <= ratio && ratio <= 1e0);
	position = position + displacement * ratio;
	displacement = displacement * (1 - ratio);
	Realvec rev_dis( reverse( displacement, v_num[0], v_num[1]) );
	displacement = rev_dis;
      }while( !still_in_the_face(position, displacement) );
    }
    return position + displacement; 
  };

  virtual bool still_in_the_face( const Realvec& position,
				  const Realvec& displacement )
  {
    //caution: sometimes rotation andle is nearly equal 2pi
    //         (typically in order of 1e-15 )
    //         but less than 2pi so this func returns 0
    //         to treat this, use epsilon
    Real epsilon(1e-8);
    Realvec temppos(position + displacement);
    Realvec v0p( vertexs.at(0) - temppos );
    Realvec v1p( vertexs.at(1) - temppos );
    Realvec v2p( vertexs.at(2) - temppos );
    Real lv0p( length(v0p) );
    Real lv1p( length(v1p) );
    Real lv2p( length(v2p) );

    double rot;
    rot = acos( dot_product( v0p, v1p ) / lv0p / lv1p );
    rot += acos( dot_product( v1p, v2p ) / lv1p / lv2p );
    rot += acos( dot_product( v2p, v0p ) / lv2p / lv0p );

    return rot >= 2 * M_PI - epsilon;
  };

  // 0: not cross
  // 1: cross
  // -1 on the edge
  int is_cross( const Realvec& position, const Realvec& displacement, int vertex_num1, int vertex_num2 )
  {
  //upon assumption that position vector is on the face
  //and displacement vector is also on the same face
    THROW_UNLESS( std::invalid_argument,
		  vertex_num1<3 && vertex_num2<3 );
    Realvec v12( vertexs.at(vertex_num2) - vertexs.at(vertex_num1) );
    Realvec v1p( position - vertexs.at(vertex_num1));
    Realvec v1pd( position + displacement - vertexs.at(vertex_num1));
    Realvec cross1( cross_product(v12, v1p) );
    Realvec cross2( cross_product(v12, v1pd) );

    Real cross( dot_product(cross1, cross2) );
    if( cross < 0 )
    {
      return 1;
    }else if( cross > 0 ){
      return 0;
    }else{
      return -1;
    }
  }

  Real cross_point( const Realvec& position, const Realvec& displacement, int v_num1, int v_num2 )
  {
  //displacement and the edge must cross each other
  //use is_cross() for assure
    Realvec v12( vertexs.at(v_num2) - vertexs.at(v_num1));
    Realvec v1p( position - vertexs.at(v_num1));
    Realvec v1pd( position + displacement - vertexs.at(v_num1));
    Real d1( length( cross_product(v12, v1p) ) / length(v12) );
    Real d2( length( cross_product(v12, v1pd)) / length(v12) );
    Real ratio( d1 / (d1 + d2) );

    return ratio;
  }
  
  Realvec reverse( const Realvec& target, int v_num1,
		   int v_num2 )
  {
    Realvec v12( vertexs.at(v_num2) - vertexs.at(v_num1) );
    Realvec v12_n( v12 / length( v12 ) );

    Real d( 2e0 * dot_product( target, v12_n ) );
    Realvec retvec( v12_n * d - target );
  
    return retvec;
  };

  virtual Realvec get_vertex()
  {
    return vertexs.at(0); 
  };

};
#endif /*FACE_SYMTRY*/
