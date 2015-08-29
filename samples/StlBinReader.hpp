#ifndef STL_BIN_READER
#define STL_BIN_READER
#include <iostream>
#include <fstream>
#include <utility>
#include <string>
#include <vector>
#include <list>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <boost/array.hpp>
#include "Defs.hpp"
#include "Vector3.hpp"
#include "Polygon.hpp"
#include "FaceBase.hpp"
#include "FaceAllGate.hpp"
// include "FaceBorder.hpp"

class StlBinReader
{

  typedef std::pair<Realvec, Realvec> vtxpair;

  struct TriangleBase
  {
    std::vector<bool> is_gate;
    int face_id;
    Realvec normal;
    std::vector<Realvec> vertexs;
    
    TriangleBase(int id): is_gate(3), face_id(id), vertexs(3){}
    vtxpair get_edgepair(const int i)
    { return std::make_pair( vertexs.at(i), vertexs.at( (i+1) % 3 ) ); }
  };

  struct edge_faces
  {
    int first_face_id;
    int second_face_id;
    vtxpair edge;

    edge_faces(const int id, const vtxpair& edgepair)
      : first_face_id(id), second_face_id(-1), edge(edgepair){}
    bool is_pair(){ return second_face_id != -1; }
    void set_pair(int id){ second_face_id = id; }
  };

  bool already_read;
  bool polygon_close;
  unsigned int numTriangle;
  std::ifstream stlfile;
  std::vector<TriangleBase> faces;
  std::vector<edge_faces> connection;

public:

  StlBinReader(std::string filename)
    : already_read(false), stlfile(filename.c_str(), std::ios::in | std::ios::binary)
  {
    if( stlfile.fail() )
      throw std::invalid_argument("file open error");
  }

  void read_file();

  std::pair<boost::shared_ptr<Polygon>, std::vector<FaceBase_sptr> > get_polygon();

private:

  Realvec to_vec(const char* vecbin);

  TriangleBase to_tri(const int id, const char* facetbin);

  //returns the object is close(true) or open(false).
  //write std::vector<edges_faces> connection.
  bool detect_connection();

  bool same_edge(const vtxpair& edge1, const vtxpair& edge2);

  bool is_same_vec(const Realvec& lhs, const Realvec& rhs, const Real tol = 1e-12);

  std::vector<boost::shared_ptr<FaceBase> > get_faces_close();

//   std::pair<boost::shared_ptr<Polygon>, std::vector<FaceBase_sptr> > get_polygon();

  //for debug.
  void dump_connection();
};

void StlBinReader::read_file()
{
  if(already_read)
     throw std::invalid_argument("this object has already read stl file");

  char header[80];
  stlfile.read(header, 80);
  std::cout << "header: " << header << std::endl;

  char numTri[4];
  stlfile.read(numTri, 4);
  numTriangle =  *(unsigned long*)numTri;
  std::cout << "number of triangles: " << numTriangle << std::endl;

  for(unsigned int i(0); i<numTriangle; ++i)
  {
    char facet[50];
    stlfile.read(facet, 50);
    faces.push_back( to_tri( i, facet ) );
  }

  polygon_close = detect_connection();
  std::cout << "polygon close: " << polygon_close << std::endl;

  if(polygon_close)
  {
//     if(numTriangle * 3 / 2 != connection.size() )
//     {
//       std::cout << "number of edges and number of faces are invalid: " << std::endl;
      std::cout << "faces: " << numTriangle << " faces * 1.5: " << numTriangle * 3 / 2 << " edges:  " << connection.size() << std::endl;
//     }
  }

//   dump_connection();

  already_read = true;
  return;
}

bool StlBinReader::detect_connection()
{
  if(faces.empty())
    throw std::invalid_argument("file not read yet");

  edge_faces edge0_0( faces.at(0).face_id, faces.at(0).get_edgepair(0) );
  edge_faces edge0_1( faces.at(0).face_id, faces.at(0).get_edgepair(1) );
  edge_faces edge0_2( faces.at(0).face_id, faces.at(0).get_edgepair(2) );
  edge_faces face0[3] = {edge0_0, edge0_1, edge0_2};

  std::list<edge_faces> edge_face_list(face0, face0+3);
//   std::vector<edge_faces> 

  for(unsigned int face_num(1); face_num<numTriangle; ++face_num)
  {
    for(int edge_num(0); edge_num<3; ++edge_num)
    {
      bool is_new(true);
      for(std::list<edge_faces>::iterator itr = edge_face_list.begin();
          itr != edge_face_list.end(); ++itr)
      {
        if( same_edge(faces.at(face_num).get_edgepair(edge_num), itr->edge) ) 
        {
          itr->set_pair( faces.at(face_num).face_id );
          connection.push_back(*itr);
          is_new = false;
          itr = edge_face_list.erase(itr);
          break;
        }
      }

      if(is_new)
      {
        edge_faces tempedge(faces.at(face_num).face_id, faces.at(face_num).get_edgepair(edge_num) );
        edge_face_list.push_back(tempedge);
      }
    }
  }

  if(!edge_face_list.empty())
  {
    for(std::list<edge_faces>::iterator itr = edge_face_list.begin();
        itr != edge_face_list.end(); ++itr)
    {
      std::cout << "first id: " << itr->first_face_id << " second id: " << itr->second_face_id << std::endl;
      std::cout << "edge: " << itr->edge.first << " -> " << itr->edge.second << std::endl; 
    }
  }

  return edge_face_list.empty();
}

std::pair<boost::shared_ptr<Polygon>, std::vector<FaceBase_sptr> >
StlBinReader::get_polygon()
{
  if(polygon_close)
  {
    std::vector<boost::shared_ptr<FaceBase> > retfaces(get_faces_close() );
    boost::shared_ptr<Polygon> poly_sptr(new Polygon(retfaces.at(0) ) );

    for(unsigned int i(1); i<numTriangle; ++i)
      poly_sptr->insert(retfaces.at(i));

    for(unsigned int i(0); i<numTriangle; ++i)
    {
      for(int edge_num(0); edge_num < 3; ++edge_num)
      {
        //TODO
        vtxpair tempedge( faces.at(i).get_edgepair(edge_num) );
        for( std::vector<edge_faces>::iterator itr = connection.begin(); itr != connection.end();
             ++itr)
        {
          if( same_edge(tempedge, itr->edge ) )
          {
            if( faces.at(i).face_id == itr->first_face_id )
              poly_sptr->set_neighbor(edge_num, retfaces.at(i), poly_sptr->id_to_faceptr(itr->second_face_id) );
            else if( faces.at(i).face_id == itr->second_face_id )
              poly_sptr->set_neighbor(edge_num, retfaces.at(i), poly_sptr->id_to_faceptr(itr->first_face_id) );
            else
              throw std::invalid_argument("internal error.");
          }
        }
      }
    }

    for(unsigned int i(0); i<numTriangle; ++i)
      retfaces.at(i)->set_poly_ptr(poly_sptr);
    
    poly_sptr->set_near_vertex();
    return std::make_pair(poly_sptr, retfaces);
  }
  else
  {
    throw std::invalid_argument("open polygon not supported yet...");
//     std::vector<boost::shared_ptr<FaceBase> > retfaces(get_faces_close() );
//     boost::shared_ptr<Polygon> poly_sptr(new Polygon(retfaces.at(0) ) );
//
//     std::cout << "polygon not closed" << std::endl;
//     //TODO
//     return make_pair(retfaces, poly_sptr);
  }
}

std::vector<boost::shared_ptr<FaceBase> > StlBinReader::get_faces_close()
{
  std::vector<FaceBase_sptr> retFBsptr;

  for(unsigned int i(0); i<numTriangle; ++i)
  {
    FaceBase_sptr tempTriangle(new FaceAllGate( faces.at(i).face_id, faces.at(i).vertexs.at(0), faces.at(i).vertexs.at(1), faces.at(i).vertexs.at(2) ) );
    retFBsptr.push_back(tempTriangle);
  }

  return retFBsptr;
}

//reads 50 bytes
StlBinReader::TriangleBase StlBinReader::to_tri(int id, const char* facetbin)
{
  TriangleBase retTri(id);

  retTri.normal = to_vec(facetbin);
  retTri.vertexs.at(0) = to_vec(facetbin+12);
  retTri.vertexs.at(1) = to_vec(facetbin+24);
  retTri.vertexs.at(2) = to_vec(facetbin+36);

//the last 2 bytes have no meaning
  return retTri;
}

//reads 12 bytes
Realvec StlBinReader::to_vec(const char* vecbin)
{
  char fl0[4] = {vecbin[0], vecbin[1], vecbin[2], vecbin[3] };
  char fl1[4] = {vecbin[4], vecbin[5], vecbin[6], vecbin[7] };
  char fl2[4] = {vecbin[8], vecbin[9], vecbin[10], vecbin[11] };

  float x = *( (float*)fl0 );
  float y = *( (float*)fl1 );
  float z = *( (float*)fl2 );

  Real xr( (double)x );
  Real yr( (double)y );
  Real zr( (double)z );

  Realvec retvec( xr, yr, zr );
//   std::cout << retvec << std::endl;
  return retvec;
}

bool StlBinReader::same_edge(const StlBinReader::vtxpair& edge1, const StlBinReader::vtxpair& edge2)
{
  bool reverse( is_same_vec(edge1.first, edge2.second)&&is_same_vec(edge1.second, edge2.first) );
  if(reverse) return true;
  bool forward( is_same_vec(edge1.first, edge2.first)&&is_same_vec(edge1.second, edge2.second) );
  return forward;
}

bool StlBinReader::is_same_vec(const Realvec& lhs, const Realvec& rhs, const Real tol)
{
  if( fabs(lhs[0] - rhs[0]) > tol ) return false;
  if( fabs(lhs[1] - rhs[1]) > tol ) return false;
  if( fabs(lhs[2] - rhs[2]) > tol ) return false;
  return true;
}

void StlBinReader::dump_connection()
{
  std::cout << "number of edges: " << connection.size() << std::endl;
  for(std::vector<edge_faces>::iterator itr = connection.begin(); itr != connection.end();
      ++itr )
  {
    std::cout << "edge: " << itr->edge.first << " -> " << itr->edge.second << std::endl;
    std::cout << "first id: " << itr->first_face_id << " second id: " << itr->second_face_id << std::endl;
  }
  return;
}


#endif
