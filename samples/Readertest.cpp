#include <boost/random.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "../GreensFunction2DAbsSym.hpp"
#include "StlFileReader.hpp"
#include "FaceAllGate.hpp"
#include "OneParticle.hpp"
#include "Polygon.hpp"

using namespace greens_functions;

int main(int argc, char* argv[])
{
  if(argc != 4)
  {
    std::cout << "please input file type(ASCII(-asc) or binary(-bin), path to stl file and initial face id." << std::endl;
    return -1;
  }

  std::string filetype( argv[1] );
  std::string filename( argv[2] );
  int initial_face( atoi(argv[3]) );

  StlFileReader reader(filename, filetype);
  reader.read_file();

  std::pair<boost::shared_ptr<Polygon>, std::vector<FaceBase_sptr> > ptrpair( reader.get_polygon() );
  
  //457
  Realvec initial_position( ptrpair.second.at(initial_face)->get_center_mass() );

  OneParticle particle(0, initial_position, ptrpair.second.at(initial_face) );

//   boost::random::mt19937 mt( static_cast<unsigned long>(time(0)) );
  boost::random::mt19937 mt( 0 );
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

  Real index(1e0);
  Real last_time(0e0);
    do
    {
      Real D(1e0);
//       a = particle.get_max_a();
//       a *= 0.5;
      // in order to treat vertex
      a = 1e0;
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

      if( (t < index ) && ( index < t+last_time ) )
      {
        fout << "\n\n";
        fout << particle << " " << std::setw(25) << t << " " << std::setw(25) << a << std::endl;
 
        index = index * 2e0;
      }

      last_time = t;

    }while(t < t_end);

    std::cout << "end." << std::endl;
  return 0;
}
