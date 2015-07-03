#ifndef POLYGON_HPP
#define POLYGON_HPP
#include <map>
#include <cmath>
#include <boost/shared_ptr.hpp>
#include "face_base.hpp"

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

class polygon
{
private:
  std::map< int, boost::shared_ptr<face_base> > face_map;
  std::map< id_gateway, boost::shared_ptr<face_base> > gateway_map;

public:
  polygon( boost::shared_ptr<face_base> ptr )
  {
    int id( ptr->get_id() );
    face_map[id] = ptr;
  }

  void insert( boost::shared_ptr<face_base> ptr )
  {
    int id( ptr->get_id() );

    std::map< int, boost::shared_ptr<face_base> >::iterator itr;
    itr = face_map.find(id);
    bool id_duplicate( itr == face_map.end() );
    THROW_UNLESS(std::invalid_argument, id_duplicate);

    face_map[id] = ptr;
    return;
  };

  void set_neighbor( const int& gate, const boost::shared_ptr<face_base>& ptr0, const boost::shared_ptr<face_base>& ptr1 )
  {
    id_gateway ig( ptr0->get_id(), gate );

    std::map< id_gateway, boost::shared_ptr<face_base> >::iterator itr;
    itr = gateway_map.find(ig);
    bool id_gateway_duplicate( itr == gateway_map.end() );
    THROW_UNLESS(std::invalid_argument, id_gateway_duplicate);

    gateway_map[ig] = ptr1;
    return;
  }

  boost::shared_ptr<face_base> get_ptr_from_id(int id)
  {
    std::map< int, boost::shared_ptr<face_base> >::iterator itr;
    itr = face_map.find(id);
    bool find_face( itr != face_map.end() );
    THROW_UNLESS(std::invalid_argument, find_face);

    return face_map[id];
  };

  boost::shared_ptr<face_base> get_neighbor_ptr_from_gateway( const int& id, const int& gate)
  {
    id_gateway ig(id, gate);

    std::map< id_gateway, boost::shared_ptr<face_base> >::iterator itr;
    itr = gateway_map.find(ig);
    bool find_neighbor_face( itr != gateway_map.end() );
    THROW_UNLESS(std::invalid_argument, find_neighbor_face);

    return gateway_map[ig];
  }

};

#endif /*POLYGON_HPP*/
