#include <boost/random.hpp>
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "../GreensFunction2DAbsSym.hpp"
#include "StlFileReader.hpp"
#include "OneParticle.hpp"
#include "Polygon.hpp"

using namespace greens_functions;

const Real RBD_STEP_INTERVAL(1e-7);

struct RandomNumberGen
{
    boost::shared_ptr<gsl_rng> _rng;

    RandomNumberGen(gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937))
        : _rng(rng, &gsl_rng_free)
    {
        ;
    }

    std::pair<Real, Real> RBD_disp(Real D, Real dt)
    {
        Real sigma (std::sqrt(2e0 * D * dt));
        Realvec disp (gaussian(sigma), gaussian(sigma), 0e0);
        return std::make_pair(length(disp), atan2(disp[0], disp[1]));
    }

    Real gaussian(Real sigma, Real mean = 0e0)
    {
        return gsl_ran_gaussian(_rng.get(), sigma) + mean;
    }
};

int main(int argc, char* argv[])
{
    if(argc != 4)
    {
        std::cout << "Usage: ./sim -[stltype] [path/to/stl] <initial face_id>"
                  << std::endl;
        std::cout << "please input file type(ASCII(-asc) or binary(-bin),"
                  << "path to stl file and initial face id." << std::endl;
        return -1;
    }

    std::string filetype( argv[1] );
    std::string filename( argv[2] );
    int initial_face( atoi(argv[3]) );
 
    StlFileReader reader(filename, filetype);
    reader.read_file();
 
    std::pair<boost::shared_ptr<Polygon>, std::vector<FaceBase_sptr> >
        ptrpair( reader.get_polygon() );
    
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
 
    Real t(0);
    Real D = 1e0;
    RandomNumberGen rng;
 
    fout << particle << " " << std::setw(25) << t << std::endl; //a is not defined
 
    do
    {
        Real a = particle.get_max_a();
        Real dt;

        if(a > VERTEX_THRESHOLD)
        {
            //greensfunction
            GreensFunction2DAbsSym gf(D,a);
            Real r;
 
            dt = gf.drawTime( rand(mt) );
 
            if(t + dt > t_end)
            {
              dt = t_end - t;
              r = gf.drawR( rand(mt), dt );
            }
            else
            {
              r = a;
            }
            Real theta(2 * M_PI * rand(mt));
 
            particle.move( r, theta );
        }
        else
        {
            //Reaction Brownian dynamics
            if(t + RBD_STEP_INTERVAL > t_end)
            {
                dt = t_end - t;
            }
            else
            {
                dt = RBD_STEP_INTERVAL;
            }

            std::pair<Real, Real> disp = rng.RBD_disp(D, dt);

            particle.move( disp.first, disp.second );
        }
 
        t += dt;
 
        fout << particle << " " << std::setw(25) << t << " " << std::setw(25) << a << std::endl;
 
    }while(t < t_end);
 
    std::cout << "end." << std::endl;
    return 0;
}
