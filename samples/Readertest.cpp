#include <boost/random.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../GreensFunction2DAbsSym.hpp"
#include "StlBinReader.hpp"
#include "FaceAllGate.hpp"
#include "OneParticle.hpp"
#include "Polygon.hpp"

using namespace greens_functions;

int main(int argc, char* argv[])
{
  if(argc != 3)
  {
    std::cout << "please input file type(ASCII(-asc) or binary(-bin) and path to stl file." << std::endl;
    return -1;
  }

  std::string filetype( argv[1] );
  std::string filename( argv[2] );

//   if(filetype == "-bin")
//   {
  if(filetype != "-bin") throw std::invalid_argument("bin stl file only...");
    StlBinReader reader(filename);
    reader.read_file();
//   }
//   else
//   {
//     std::cout << filetype << std::endl;
//   }

  std::pair<boost::shared_ptr<Polygon>, std::vector<FaceBase_sptr> > ptrpair( reader.get_polygon() );
  
  Realvec initial_position( ptrpair.second.at(0)->get_center_mass() );

  OneParticle particle(0, initial_position, ptrpair.second.at(0));

  boost::random::mt19937 mt(0);
  boost::random::uniform_real_distribution<Real> rand(0.0, 1.0);

  Real t_end;
  std::cout << "t_end: ";
  std::cin >> t_end;

  std::string output;
  std::cout << "output: ";
  std::cin >> output;
  std::ofstream fout( output.c_str() );

  Real t(0), dt(0), theta(0), r(0), a(0);

  fout << particle << " " << std::setw(25) << t << " " << std::setw(25) << a << std::endl;

    do
    {
      Real D(1e0);
      a = particle.get_max_a();
      GreensFunction2DAbsSym gf(D,a);

      dt = gf.drawTime( rand(mt) );
//       face_time.at(particle.get_face_id()) += dt;
   
      if(t + dt > t_end)
      {
      	dt = t_end - t;
      	r = gf.drawR( rand(mt), dt );
      }else{
      	r = a;
      }

      if( particle.involve_vertex() )
      {
        //TODO
        theta = 2 * M_PI * rand(mt);
      }else{
        theta = 2 * M_PI * rand(mt);
      }

      particle.move( r, theta );

      t += dt;

      fout << particle << " " << std::setw(25) << t << " " << std::setw(25) << a << std::endl;

    }while(t < t_end);

  return 0;
}