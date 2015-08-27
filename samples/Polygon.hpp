#ifndef POLYGON_HPP
#define POLYGON_HPP
#include <map>
#include <cmath>
#include <boost/shared_ptr.hpp>
#include "FaceBase.hpp"

class Polygon
{
private:
  std::map< int, boost::shared_ptr<FaceBase> > face_map;
  std::map< std::pair<int, int>, boost::shared_ptr<FaceBase> > gateway_map;

public:
  Polygon()
  {
  }

  Polygon( boost::shared_ptr<FaceBase> ptr )
  {
    int id( ptr->get_id() );
    face_map[id] = ptr;
  }

  void insert( boost::shared_ptr<FaceBase> ptr );
  
  void set_neighbor( const int gate, const boost::shared_ptr<FaceBase>& ptr0, const boost::shared_ptr<FaceBase>& ptr1 );

// use this function after all faces are inserted and all gateway edges are connected
  void set_near_vertex();
 
  FaceBase_sptr id_to_faceptr(int id);

  FaceBase_sptr get_neighbor( const int id, const int gate);

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

void Polygon::set_neighbor( const int gate, const boost::shared_ptr<FaceBase>& ptr0, const boost::shared_ptr<FaceBase>& ptr1 )
{
  std::pair<int, int> ig( ptr0->get_id(), gate );

  std::map< std::pair<int, int>, boost::shared_ptr<FaceBase> >::iterator itr;
  itr = gateway_map.find(ig);

  if( itr != gateway_map.end() )
    throw std::invalid_argument("polygon::set_neighbor: id duplicated");

  gateway_map[ig] = ptr1;
  return;
};

void Polygon::set_near_vertex()
{
  for(std::map< int, boost::shared_ptr<FaceBase> >::iterator itr( face_map.begin() ); itr != face_map.end(); ++itr)
  {
    itr->second->set_near_vertexs();
  }
}

FaceBase_sptr Polygon::id_to_faceptr( int id )
{
  std::map< int, boost::shared_ptr<FaceBase> >::iterator itr;
  itr = face_map.find(id);

  if(itr == face_map.end() )
    throw std::invalid_argument("Polygon::id_to_faceptr: no face");

  return face_map[id];
};

FaceBase_sptr Polygon::get_neighbor( const int id, const int gate )
{
  std::pair<int, int> ig(id, gate);

  std::map< std::pair<int, int>, boost::shared_ptr<FaceBase> >::iterator itr;
  itr = gateway_map.find(ig);

  if( itr == gateway_map.end() )
  {
    std::cout << "face_id: " << id << " edge_id: " << gate << std::endl;
    throw std::invalid_argument("no neighbor face");
  }

  return gateway_map[ig];
}


#endif /*POLYGON_HPP*/
