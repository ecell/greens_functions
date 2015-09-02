#include "../../GreensFunction2DAbsSym.hpp"
#include "../FaceOpen.hpp"
#include "../OneParticle.hpp"
#include "../Polygon.hpp"
#include <boost/random.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace greens_functions;

int main()
{
//   Realvec v1(0e0, 0e0, 0e0);
//   Realvec v2(1e0, 0e0, 0e0);
//   Realvec v3(0e0, 1e0, 0e0);
//   Realvec v4(1e0, 1e0, 0e0);
  Realvec v1(0e0, 1e0, 1e0);
  Realvec v2(0e0, 1e0, 0e0);
  Realvec v3(0e0, 0e0, 0e0);
  Realvec v4(1e0, 1e0, 0e0);

  std::vector<bool> gate(3);
  gate.at(0) = false;
  gate.at(1) = true;
  gate.at(2) = false;

  FaceBase_sptr triangle1_ptr( new FaceOpen(0, v1, v2, v3, gate) );
  FaceBase_sptr triangle2_ptr( new FaceOpen(1, v4, v3, v2, gate) );

  boost::shared_ptr<Polygon> rectangle_ptr( new Polygon( triangle1_ptr ) );
  rectangle_ptr->insert( triangle2_ptr );
  rectangle_ptr->set_neighbor( 1, triangle1_ptr, triangle2_ptr);
  rectangle_ptr->set_neighbor( 1, triangle2_ptr, triangle1_ptr);

  triangle1_ptr->set_poly_ptr( rectangle_ptr );
  triangle2_ptr->set_poly_ptr( rectangle_ptr );

  rectangle_ptr->set_near_vertex();

//   Realvec position( 1e0/3e0, 1e0/3e0, 0e0 );
  Realvec position( 0e0, 2e0/3e0, 1e0/3e0 );

  OneParticle particle(0, position, triangle1_ptr);

  boost::random::mt19937 mt(0);
  boost::random::uniform_real_distribution<Real> rand(0.0, 1.0);

  Real t_end;
  std::cout << "t_end: ";
  std::cin >> t_end;

  std::string filename;
  std::cout << "filename: ";
  std::cin >> filename;
  std::ofstream fout( filename.c_str() );

  Real t(0), dt(0), theta(0), r(0), a(0);

  fout << particle << " " << t << " " << a << std::endl;
  
  do
  {
    Real D(1e0);
    a = particle.get_max_a();
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

    particle.move( r, theta );

    t += dt;
   
    fout << particle << " " << t << " " << a << std::endl;

  }while(t < t_end);

  return 0;
}
