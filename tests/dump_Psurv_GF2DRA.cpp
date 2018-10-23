#include "../GreensFunction2DRadAbs.hpp"
#include <boost/random.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

using greens_functions::Real;

void plot_with_ka(const Real ka, const std::string fname)
{
    std::ofstream ofs(fname);
    ofs << std::setprecision(15);

    const Real D  = 1.0;
    const Real a  = 10.0;
    const Real s  = 1.0;
    const Real r0 = 5.0;
    const Real k  = ka * 2 * 3.1416 * s * D;

    const greens_functions::GreensFunction2DRadAbs gf(D, k, r0, s, a);

    const Real dt = 1e-2;
    for(std::size_t i=1; i<=5000; ++i)
    {
        ofs << dt * i << ' ' << gf.p_survival(dt * i) << '\n';
    }

    return;
}

int main()
{
    const static std::size_t N = 1000;

    boost::random::mt19937 mt(123456789);
    boost::random::uniform_01<Real> canonical;

    boost::random::uniform_real_distribution<Real> D_gen(1e-5, 1e-3),
                                                   k_gen(1e-5, 1e-3),
                                                   a_gen(1e-5, 1e-3);


    // ------------------------------------------------------------------ //
    // drawTime                                                           //
    // ------------------------------------------------------------------ //
    plot_with_ka(1e-5, "dump_GreensFunction2DRadAbs_ka_1e-5.dat");
    plot_with_ka(1e-4, "dump_GreensFunction2DRadAbs_ka_1e-4.dat");
    plot_with_ka(1e-3, "dump_GreensFunction2DRadAbs_ka_1e-3.dat");
    plot_with_ka(1e-2, "dump_GreensFunction2DRadAbs_ka_1e-2.dat");
    plot_with_ka(1e-1, "dump_GreensFunction2DRadAbs_ka_1e-1.dat");
    plot_with_ka(1e+0, "dump_GreensFunction2DRadAbs_ka_1e+0.dat");
    plot_with_ka(1e+1, "dump_GreensFunction2DRadAbs_ka_1e+1.dat");
    plot_with_ka(1e+2, "dump_GreensFunction2DRadAbs_ka_1e+2.dat");
    plot_with_ka(1e+3, "dump_GreensFunction2DRadAbs_ka_1e+3.dat");
    plot_with_ka(1e+4, "dump_GreensFunction2DRadAbs_ka_1e+4.dat");
    plot_with_ka(1e+5, "dump_GreensFunction2DRadAbs_ka_1e+5.dat");



    return 0;
}
