#include "../GreensFunction2DAbsSym.hpp"
#include "face_infty.hpp"
#include <boost/random.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace greens_functions;

class particle 
{
  Realvec position;
  int particle_id;
  boost::shared_ptr<face_base> face_ptr;
public:
  particle( int id, const Realvec& pos, boost::shared_ptr<face_base> ptr)
  {
    particle_id = id;
    position = pos;
    face_ptr = ptr;
  }
  void move(const Realvec& displacement);
  friend std::ostream& operator<<(std::ostream& os, const particle& part);
};

void particle::move( const Realvec& displacement )
{
  position = face_ptr->move(position, displacement, face_ptr);
  return;
}

std::ostream& operator<<( std::ostream& os, const particle& part)
{
//  os << "particle id: " << part.particle_id << std::endl;
//  os << "position: ";
  os << part.position[0] << " ";
  os << part.position[1] << " ";
  os << part.position[2] << " ";
  return os;
}

int main()
{
//   face_infty infty_plane(0); // double free or corruption (out))
//   boost::shared_ptr<face_infty> ptr( &infty_plane ); // but successfully run for program end
  boost::shared_ptr<face_base> infplane_ptr( new face_infty(0) );

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
  
  Realvec displacement(3);

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

    displacement[0] = r * cos(theta);
    displacement[1] = r * sin(theta);
    displacement[2] = 0;

    fout << "displacement: " << displacement << std::endl;

    mol.move( displacement );

    t += dt;
   
    fout << mol << t << " " << a << std::endl;

  }while(t < t_end);

  return 0;
}
