#include "../../GreensFunction2DAbsSym.hpp"
#include "../FaceTwoGate.hpp"
#include "../singleton.hpp"
#include "../Polygon.hpp"
#include <boost/random.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace greens_functions;

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

  face_sptr tryangle0_ptr( new FaceTwoGate(0, v0, v2, v4, 0, 2) );
  face_sptr tryangle1_ptr( new FaceTwoGate(1, v0, v4, v1, 0, 1) );
  face_sptr tryangle2_ptr( new FaceTwoGate(2, v1, v4, v7, 0, 2) );
  face_sptr tryangle3_ptr( new FaceTwoGate(3, v1, v7, v6, 0, 1) );
  face_sptr tryangle4_ptr( new FaceTwoGate(4, v6, v7, v5, 0, 2) );
  face_sptr tryangle5_ptr( new FaceTwoGate(5, v6, v5, v3, 0, 1) );
  face_sptr tryangle6_ptr( new FaceTwoGate(6, v3, v5, v2, 0, 2) );
  face_sptr tryangle7_ptr( new FaceTwoGate(7, v3, v2, v0, 0, 1) );

  boost::shared_ptr<Polygon> cube_ptr( new Polygon( tryangle0_ptr ) );
  
  cube_ptr->insert( tryangle1_ptr );
  cube_ptr->insert( tryangle2_ptr );
  cube_ptr->insert( tryangle3_ptr );
  cube_ptr->insert( tryangle4_ptr );
  cube_ptr->insert( tryangle5_ptr );
  cube_ptr->insert( tryangle6_ptr );
  cube_ptr->insert( tryangle7_ptr );

  cube_ptr->set_neighbor( 0, tryangle0_ptr, tryangle7_ptr );
  cube_ptr->set_neighbor( 2, tryangle0_ptr, tryangle1_ptr );
  cube_ptr->set_neighbor( 0, tryangle1_ptr, tryangle0_ptr );
  cube_ptr->set_neighbor( 1, tryangle1_ptr, tryangle2_ptr );
  cube_ptr->set_neighbor( 0, tryangle2_ptr, tryangle1_ptr );
  cube_ptr->set_neighbor( 2, tryangle2_ptr, tryangle3_ptr );
  cube_ptr->set_neighbor( 0, tryangle3_ptr, tryangle2_ptr );
  cube_ptr->set_neighbor( 1, tryangle3_ptr, tryangle4_ptr );
  cube_ptr->set_neighbor( 0, tryangle4_ptr, tryangle3_ptr );
  cube_ptr->set_neighbor( 2, tryangle4_ptr, tryangle5_ptr );
  cube_ptr->set_neighbor( 0, tryangle5_ptr, tryangle4_ptr );
  cube_ptr->set_neighbor( 1, tryangle5_ptr, tryangle6_ptr );
  cube_ptr->set_neighbor( 0, tryangle6_ptr, tryangle5_ptr );
  cube_ptr->set_neighbor( 2, tryangle6_ptr, tryangle7_ptr );
  cube_ptr->set_neighbor( 0, tryangle7_ptr, tryangle6_ptr );
  cube_ptr->set_neighbor( 1, tryangle7_ptr, tryangle0_ptr );

  tryangle0_ptr->set_belonging_polygon( cube_ptr );
  tryangle1_ptr->set_belonging_polygon( cube_ptr );
  tryangle2_ptr->set_belonging_polygon( cube_ptr );
  tryangle3_ptr->set_belonging_polygon( cube_ptr );
  tryangle4_ptr->set_belonging_polygon( cube_ptr );
  tryangle5_ptr->set_belonging_polygon( cube_ptr );
  tryangle6_ptr->set_belonging_polygon( cube_ptr );
  tryangle7_ptr->set_belonging_polygon( cube_ptr );

//   Realvec position( 1e0/3e0, 1e0/3e0, 0e0 );
  Realvec position( 1e0/3e0, 2e0/3e0, 0e0 );
//   std::cout << position << std::endl;

  particle mol(0, position, tryangle0_ptr);

  boost::random::mt19937 mt(0);
  boost::random::uniform_real_distribution<Real> rand(0.0, 1.0);

  Real a_max, t_end;
  std::cout << "a_max: ";
  std::cin >> a_max;
  std::cout << "t_end: ";
  std::cin >> t_end;

  std::string filename;
  std::cout << "filename: ";
  std::cin >> filename;
  std::ofstream fout( filename.c_str() );

  Real t(0), dt(0), theta(0), r(0), a(0);

  fout << mol << " " << t << " " << a << std::endl;
  
  do
  {
    Real D(1e0);
    a = a_max * rand(mt);
    GreensFunction2DAbsSym gf(D,a);

    dt = gf.drawTime( rand(mt) );
 
    if(t + dt > t_end)
    {
      dt = t_end - t;
      r = gf.drawR( rand(mt), dt );
    }else{
      r = a;
    }

    theta = 2 * M_PI * rand(mt);

    mol.move( r, theta );

    t += dt;
   
    fout << mol << " " << t << " " << a << std::endl;

  }while(t < t_end);

  return 0;
}
