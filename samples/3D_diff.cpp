#include "GreensFunction3DAbsSym.hpp"
#include <boost/random.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

const int sample(300);
const int plot(100);

int main()
{
  boost::random::mt19937 mt(0);
  boost::random::uniform_real_distribution<Real> rand(0.0, 1.0);

  Real a_max, t_end, dt_end;
  std::cout << "a: ";
  std::cin >> a_max;
  std::cout << "first_t: ";
  std::cin >> t_end;
  std::cout << "dt_end: ";
  std::cin >> dt_end;

  std::string filename;
  std::cout << "filename: ";
  std::cin >> filename;
  std::ofstream fout( filename.c_str() );

  Real t(0), dt(0), r(0), r_tot(0), a(0);
  Real theta(0), phi(0), x(0), y(0),z(0);
//  fout << x << " " << y << " " << t << std::endl;
// /*
  for(int i(0); i < plot; ++i)
  {
    r_tot = 0;
    for(int j(0); j < sample; ++j)
    {
      t = 0;
      x = 0;
      y = 0;
      z = 0;
// */

      do
      {
	Real D(1e0);
	a = a_max * rand(mt);//[0,a_max)
	greens_functions::GreensFunction3DAbsSym gf(D,a);

	dt = gf.drawTime( rand(mt) );
       
	if(t + dt > t_end)
	{
	  dt = t_end - t;
	  r = gf.drawR( rand(mt), dt );
	}else{
	  r = a;
	}

	theta = 2 * M_PI * rand(mt);
	phi = 2 * M_PI * rand(mt);

	x += r * sin(theta) * cos(phi);
	y += r * sin(theta) * sin(phi);
	z += r * cos(theta);
	t += dt;
       
//        fout << x << " " << y << " " << t << " " << a << std::endl;
	
      }while(t < t_end);
// /*
      r_tot += r;
    }

    fout << r_tot / sample << " " << t << std::endl;
    t_end += dt_end;

  }
// */  
  return 0;
}
