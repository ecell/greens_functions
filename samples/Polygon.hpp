#ifndef POLYGON_HPP
#define POLYGON_HPP
#include <map>
#include <cmath>
#include <boost/shared_ptr.hpp>
#include "FaceBase.hpp"

using namespace greens_functions;

struct id_gateway{
  int face_id;
  int gateway;
  id_gateway(int id, int gate): face_id(id), gateway(gate){};
};

bool operator<(const id_gateway& lhs, const id_gateway& rhs)
{
  if(lhs.face_id == rhs.face_id){
    return lhs.gateway < rhs.gateway;
  }
  return lhs.face_id < rhs.face_id;
}

bool operator>(const id_gateway& lhs, const id_gateway& rhs)
{
  if(lhs.face_id == rhs.face_id){
    return lhs.gateway > rhs.gateway;
  }
  return lhs.face_id > rhs.face_id;
}

class Polygon
{
private:
  std::map< int, boost::shared_ptr<FaceBase> > face_map;
  std::map< id_gateway, boost::shared_ptr<FaceBase> > gateway_map;

public:
  Polygon( boost::shared_ptr<FaceBase> ptr )
  {
    int id( ptr->get_id() );
    face_map[id] = ptr;
  }

  void insert( boost::shared_ptr<FaceBase> ptr );
  
  void set_neighbor( const int& gate, const boost::shared_ptr<FaceBase>& ptr0, const boost::shared_ptr<FaceBase>& ptr1 );

// use this function after all faces are inserted and all gateway edges are connected
//   void set_near_vertex();
 
  boost::shared_ptr<FaceBase> get_ptr_from_id(int id);

  boost::shared_ptr<FaceBase> get_neighbor_ptr_from_gateway( const int& id, const int& gate);

// private:

};

void Polygon::insert( boost::shared_ptr<FaceBase> ptr )
{
  int id( ptr->get_id() );

  std::map< int, boost::shared_ptr<FaceBase> >::iterator itr;
  itr = face_map.find(id);
  bool id_duplicate( itr == face_map.end() );
  THROW_UNLESS(std::invalid_argument, id_duplicate);

  face_map[id] = ptr;
  return;
};

void Polygon::set_neighbor( const int& gate, const boost::shared_ptr<FaceBase>& ptr0, const boost::shared_ptr<FaceBase>& ptr1 )
{
  id_gateway ig( ptr0->get_id(), gate );

  std::map< id_gateway, boost::shared_ptr<FaceBase> >::iterator itr;
  itr = gateway_map.find(ig);
  bool id_gateway_duplicate( itr == gateway_map.end() );
  THROW_UNLESS(std::invalid_argument, id_gateway_duplicate);

  gateway_map[ig] = ptr1;
  return;
};

// void Polygon::set_near_vertex()
// {
//   for(std::map< int, boost::shared_ptr<FaceBase> >::iterator itr( face_map.begin() ); itr != face_map.end(); ++itr)
//   {
//     for(int i(0); i<3; ++i)
//     {
//       int id(itr->first);
//       boost::shared_ptr<FaceBase> fb_sptr;
//
//       if( is_gateway_edge( i ) )
//       {
// 	int id( itr->first );
//       }
//     }
//   }
// }

boost::shared_ptr<FaceBase> Polygon::get_ptr_from_id( int id )
{
  std::map< int, boost::shared_ptr<FaceBase> >::iterator itr;
  itr = face_map.find(id);
  bool find_face( itr != face_map.end() );
  THROW_UNLESS(std::invalid_argument, find_face);

  return face_map[id];
};

boost::shared_ptr<FaceBase> Polygon::get_neighbor_ptr_from_gateway( const int& id, const int& gate )
{
  id_gateway ig(id, gate);

  std::map< id_gateway, boost::shared_ptr<FaceBase> >::iterator itr;
  itr = gateway_map.find(ig);
  bool find_neighbor_face( itr != gateway_map.end() );
  THROW_UNLESS(std::invalid_argument, find_neighbor_face);

  return gateway_map[ig];
}


#endif /*POLYGON_HPP*/
