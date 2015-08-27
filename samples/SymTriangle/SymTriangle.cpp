#include "../../GreensFunction2DAbsSym.hpp"
#include "../FaceAllSym.hpp"
#include "../OneParticle.hpp"
#include <boost/random.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace greens_functions;
typedef boost::shared_ptr<FaceBase> face_sptr;

int main()
{
  Realvec v1(1e0, 0e0, 0e0);
  Realvec v2(0e0, 1e0, 0e0);
  Realvec v3(0e0, 0e0, 1e0);
//   Realvec norm( cross_product(v2 - v1, v3 - v1) );

  face_sptr triangle_ptr( new FaceAllSym(0, v1, v2, v3) );

  Realvec e12( v2 - v1 );
  Realvec e13( v3 - v1 );
  Realvec position( v1 + e12/3e0 + e13/3e0 );
//   std::cout << position << std::endl;

  particle mol(0, position, triangle_ptr);

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
