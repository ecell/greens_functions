#include "GreensFunction3DSym.hpp"
#include <boost/random.hpp>
#include <iostream>
#include <cmath>

int main()
{
    boost::random::mt19937 rng(0);
    boost::random::uniform_real_distribution<Real> rand(0.0, 1.0);

    Real const D(1.0);

    greens_functions::GreensFunction3DSym gf(D);

    Real const t(1e-5);
    std::cout << gf.drawTime(rand(rng)) << std::endl;
    std::cout << gf.drawR(rand(rng), t) << std::endl;
    std::cout << gf.ip_r(GSL_POSINF, t) << std::endl;

    return 0;
}
