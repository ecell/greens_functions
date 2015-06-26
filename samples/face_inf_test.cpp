#include "../GreensFunction2DAbsSym.hpp"
#include "face_infty.hpp"
#include "singleton.hpp"
#include <boost/random.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace greens_functions;

int main()
{
  Realvec norm(3);
    norm[0] = 0e0;
    norm[1] = 0e0;
    norm[2] = 1e0;
  Realvec rep(3);
    rep[0] = 1e0;
    rep[1] = 0e0;
    rep[2] = 0e0;

//   face_infty infty_plane(0); // double free or corruption (out))
//   boost::shared_ptr<face_infty> ptr( &infty_plane ); // but successfully run for program end
  face_sptr infplane_ptr( new face_infty(0, norm, rep) );

  Realvec position(3);
    position[0] = 0e0;
    position[1] = 0e0;
    position[2] = 0e0;

  particle mol(0, position, infplane_ptr);

  boost::random::mt19937 mt(0);
  boost::random::uniform_real_distribution<Real> rand(0.0, 1.0);

  Real a_max, t_end, dt_end;
  std::cout << "a: ";
  std::cin >> a_max;
  std::cout << "t_end: ";
  std::cin >> t_end;

  std::string filename;
  std::cout << "filename: ";
  std::cin >> filename;
  std::ofstream fout( filename.c_str() );

  Real t(0), dt(0), theta(0), r(0), r_tot(0), a(0);
  
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
   
    fout << mol << t << " " << a << std::endl;

  }while(t < t_end);

  return 0;
}
