#include "../../GreensFunction2DAbsSym.hpp"
#include "../FaceAllGate.hpp"
#include "../singleton.hpp"
#include "../Polygon.hpp"
#include <boost/random.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace greens_functions;
typedef boost::shared_ptr<FaceBase> face_sptr; 

int main()
{
  Realvec v0(0e0, 0e0, 0e0);
  Realvec v1(1e0, 0e0, 0e0);
  Realvec v2(0e0, 1e0, 0e0);
  Realvec v3(0e0, 0e0, 1e0);
  Realvec v4(1e0, 1e0, 0e0);
  Realvec v5(0e0, 1e0, 1e0);
  Realvec v6(1e0, 0e0, 1e0);
  Realvec v7(1e0, 1e0, 1e0);

  face_sptr tryangle0_ptr(  new FaceAllGate( 0, v0, v2, v4) );
  face_sptr tryangle1_ptr(  new FaceAllGate( 1, v0, v4, v1) );
  face_sptr tryangle2_ptr(  new FaceAllGate( 2, v1, v4, v7) );
  face_sptr tryangle3_ptr(  new FaceAllGate( 3, v1, v7, v6) );
  face_sptr tryangle4_ptr(  new FaceAllGate( 4, v6, v7, v5) );
  face_sptr tryangle5_ptr(  new FaceAllGate( 5, v6, v5, v3) );
  face_sptr tryangle6_ptr(  new FaceAllGate( 6, v3, v5, v2) );
  face_sptr tryangle7_ptr(  new FaceAllGate( 7, v3, v2, v0) );
  face_sptr tryangle8_ptr(  new FaceAllGate( 8, v2, v5, v7) );
  face_sptr tryangle9_ptr(  new FaceAllGate( 9, v2, v7, v4) );
  face_sptr tryangle10_ptr( new FaceAllGate(10, v3, v0, v1) );
  face_sptr tryangle11_ptr( new FaceAllGate(11, v3, v1, v6) );

  boost::shared_ptr<Polygon> cube_ptr( new Polygon( tryangle0_ptr ) );
  
  cube_ptr->insert( tryangle1_ptr );
  cube_ptr->insert( tryangle2_ptr );
  cube_ptr->insert( tryangle3_ptr );
  cube_ptr->insert( tryangle4_ptr );
  cube_ptr->insert( tryangle5_ptr );
  cube_ptr->insert( tryangle6_ptr );
  cube_ptr->insert( tryangle7_ptr );
  cube_ptr->insert( tryangle8_ptr );
  cube_ptr->insert( tryangle9_ptr );
  cube_ptr->insert( tryangle10_ptr );
  cube_ptr->insert( tryangle11_ptr );

  cube_ptr->set_neighbor( 0, tryangle0_ptr, tryangle7_ptr );
  cube_ptr->set_neighbor( 1, tryangle0_ptr, tryangle9_ptr );
  cube_ptr->set_neighbor( 2, tryangle0_ptr, tryangle1_ptr );

  cube_ptr->set_neighbor( 0, tryangle1_ptr, tryangle0_ptr );
  cube_ptr->set_neighbor( 1, tryangle1_ptr, tryangle2_ptr );
  cube_ptr->set_neighbor( 2, tryangle1_ptr, tryangle10_ptr );

  cube_ptr->set_neighbor( 0, tryangle2_ptr, tryangle1_ptr );
  cube_ptr->set_neighbor( 1, tryangle2_ptr, tryangle9_ptr );
  cube_ptr->set_neighbor( 2, tryangle2_ptr, tryangle3_ptr );

  cube_ptr->set_neighbor( 0, tryangle3_ptr, tryangle2_ptr );
  cube_ptr->set_neighbor( 1, tryangle3_ptr, tryangle4_ptr );
  cube_ptr->set_neighbor( 2, tryangle3_ptr, tryangle11_ptr );

  cube_ptr->set_neighbor( 0, tryangle4_ptr, tryangle3_ptr );
  cube_ptr->set_neighbor( 1, tryangle4_ptr, tryangle8_ptr );
  cube_ptr->set_neighbor( 2, tryangle4_ptr, tryangle5_ptr );

  cube_ptr->set_neighbor( 0, tryangle5_ptr, tryangle4_ptr );
  cube_ptr->set_neighbor( 1, tryangle5_ptr, tryangle6_ptr );
  cube_ptr->set_neighbor( 2, tryangle5_ptr, tryangle11_ptr );

  cube_ptr->set_neighbor( 0, tryangle6_ptr, tryangle5_ptr );
  cube_ptr->set_neighbor( 1, tryangle6_ptr, tryangle8_ptr );
  cube_ptr->set_neighbor( 2, tryangle6_ptr, tryangle7_ptr );

  cube_ptr->set_neighbor( 0, tryangle7_ptr, tryangle6_ptr );
  cube_ptr->set_neighbor( 1, tryangle7_ptr, tryangle0_ptr );
  cube_ptr->set_neighbor( 2, tryangle7_ptr, tryangle10_ptr );

  cube_ptr->set_neighbor( 0, tryangle8_ptr, tryangle6_ptr );
  cube_ptr->set_neighbor( 1, tryangle8_ptr, tryangle4_ptr );
  cube_ptr->set_neighbor( 2, tryangle8_ptr, tryangle9_ptr );

  cube_ptr->set_neighbor( 0, tryangle9_ptr, tryangle8_ptr );
  cube_ptr->set_neighbor( 1, tryangle9_ptr, tryangle2_ptr );
  cube_ptr->set_neighbor( 2, tryangle9_ptr, tryangle0_ptr );

  cube_ptr->set_neighbor( 0, tryangle10_ptr, tryangle7_ptr );
  cube_ptr->set_neighbor( 1, tryangle10_ptr, tryangle1_ptr );
  cube_ptr->set_neighbor( 2, tryangle10_ptr, tryangle11_ptr );

  cube_ptr->set_neighbor( 0, tryangle11_ptr, tryangle10_ptr );
  cube_ptr->set_neighbor( 1, tryangle11_ptr, tryangle3_ptr );
  cube_ptr->set_neighbor( 2, tryangle11_ptr, tryangle5_ptr );

  tryangle0_ptr->set_belonging_polygon( cube_ptr );
  tryangle1_ptr->set_belonging_polygon( cube_ptr );
  tryangle2_ptr->set_belonging_polygon( cube_ptr );
  tryangle3_ptr->set_belonging_polygon( cube_ptr );
  tryangle4_ptr->set_belonging_polygon( cube_ptr );
  tryangle5_ptr->set_belonging_polygon( cube_ptr );
  tryangle6_ptr->set_belonging_polygon( cube_ptr );
  tryangle7_ptr->set_belonging_polygon( cube_ptr );
  tryangle8_ptr->set_belonging_polygon( cube_ptr );
  tryangle9_ptr->set_belonging_polygon( cube_ptr );
  tryangle10_ptr->set_belonging_polygon( cube_ptr );
  tryangle11_ptr->set_belonging_polygon( cube_ptr );

  cube_ptr->set_near_vertex();

//   Realvec position( 1e0/3e0, 1e0/3e0, 0e0 );
  Realvec position( 1e0/3e0, 2e0/3e0, 0e0 );
//   std::cout << position << std::endl;

  particle mol(0, position, tryangle0_ptr);

  boost::random::mt19937 mt(0);
  boost::random::uniform_real_distribution<Real> rand(0.0, 1.0);

  Real t_end;
//   Real a_max;
//   std::cout << "a_max: ";
//   std::cin >> a_max;
  std::cout << "t_end: ";
  std::cin >> t_end;

  std::string filename;
  std::cout << "filename: ";
  std::cin >> filename;
  std::ofstream fout( filename.c_str() );
//   std::ofstream prof( "profile.dat" );

  Real t(0), dt(0), theta(0), r(0), a(0);

  fout << mol << " " << t << " " << a << std::endl;

//   prof << mol.get_face_id() << " ";
//   prof << t << std::endl;
//   std::vector<Real> face_time(8, 0e0);
  
//   for(int j(0); j<5; ++j)
//   {
//     t_end *= 1e1;
//     for(int i(0); i<8; ++i)
//     {
//       face_time.at(i) = 0e0;
//     }
//     t = 0e0;

    do
    {
      Real D(1e0);
      a = mol.get_max_a();
      GreensFunction2DAbsSym gf(D,a);

      dt = gf.drawTime( rand(mt) );
//       face_time.at(mol.get_face_id()) += dt;
   
      if(t + dt > t_end)
      {
	dt = t_end - t;
	r = gf.drawR( rand(mt), dt );
      }else{
	r = a;
      }

      if( mol.involve_vertex() )
      {
	//TODO
	theta = 2 * M_PI * rand(mt);
      }else{
	theta = 2 * M_PI * rand(mt);
      }

      mol.move( r, theta );

      t += dt;

//       prof << mol.get_face_id() << " ";
//       prof << t << std::endl;
     
      fout << mol << " " << t << " " << a << std::endl;

    }while(t < t_end);

//     for(int i(0); i<8; ++i)
//     {
//       prof << i << " " << face_time.at(i) / t_end * 1e2 << std::endl;
//     }
//     prof << "\n\n";
//     std::cout << j << "end." << std::endl;
//   }
//
 

  return 0;
}
