#include "../GreensFunction2DRadAbs.hpp"
#include <boost/random.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

int main()
{
    using greens_functions::Real;
    const static std::size_t N = 1000;

    boost::random::mt19937 mt(123456789);
    boost::random::uniform_01<Real> canonical;

    boost::random::uniform_real_distribution<Real> D_gen(1e-5, 1e-3),
                                                   k_gen(1e-5, 1e-3),
                                                   a_gen(1e-5, 1e-3);

    std::ofstream ofs("dump_GreensFunction2DRadAbs.dat");
    ofs << std::setprecision(15);

    // ------------------------------------------------------------------ //
    // drawTime                                                           //
    // ------------------------------------------------------------------ //
    ofs << "drawTime\n";
    {
        for(std::size_t i=0; i<N; ++i)
        {
            const Real D  = D_gen(mt);
            const Real k  = k_gen(mt);
            const Real a  = a_gen(mt);
            boost::random::uniform_real_distribution<Real> s_gen(a * 1e-2, a);
            const Real s  = s_gen(mt);
            boost::random::uniform_real_distribution<Real> r0_gen(s, a);
            const Real r0 = r0_gen(mt);

            const greens_functions::GreensFunction2DRadAbs gf(D, k, r0, s, a);

            ofs << D << ' ' << k << ' ' << a << ' ' << s << ' ' << r0 << ' ' << std::flush;
            ofs << gf.drawTime(canonical(mt)) << '\n';
        }
    }

    // ------------------------------------------------------------------------
    // drawEvent
    // ------------------------------------------------------------------------

    ofs << "drawEvent\n";
    {
        for(std::size_t i=0; i<N; ++i)
        {
            const Real D = D_gen(mt);
            const Real k = k_gen(mt);
            const Real a = a_gen(mt);
            boost::random::uniform_real_distribution<Real> s_gen(a * 1e-2, a);
            const Real s = s_gen(mt);
            boost::random::uniform_real_distribution<Real> r0_gen(s, a);
            const Real r0 = r0_gen(mt);

            const greens_functions::GreensFunction2DRadAbs gf(D, k, r0, s, a);

            const Real rnd1  = canonical(mt);
            const Real time  = gf.drawTime(rnd1);
            const Real rnd2  = canonical(mt);
            const auto event = gf.drawEventType(rnd2, time);

            switch(event)
            {
                case greens_functions::GreensFunction::IV_ESCAPE:
                {
                    ofs << "escape\n"; break;
                }
                case greens_functions::GreensFunction::IV_REACTION:
                {
                    ofs << "reaction\n"; break;
                }
                default:
                {
                    assert(false);
                }
            }
        }
    }

    // ------------------------------------------------------------------------
    // drawR
    // ------------------------------------------------------------------------

    ofs << "drawR\n";
    {
        for(std::size_t i=0; i<N; ++i)
        {
            const Real D = D_gen(mt);
            const Real k = k_gen(mt);
            const Real a = a_gen(mt);
            boost::random::uniform_real_distribution<Real> s_gen(a * 1e-2, a);
            const Real s = s_gen(mt);
            boost::random::uniform_real_distribution<Real> r0_gen(s, a);
            const Real r0 = r0_gen(mt);

            const greens_functions::GreensFunction2DRadAbs gf(D, k, r0, s, a);

            const Real rnd1 = canonical(mt);
            const Real time = gf.drawTime(rnd1);
            const Real rnd2 = canonical(mt);
            const Real R    = gf.drawR(rnd2, time);

            ofs << D << ' ' << k << ' ' << a << ' ' << s << ' ' << r0 << ' '
                      << time << ' ' << R << '\n';
        }
    }

    // ------------------------------------------------------------------------
    // drawTheta
    // ------------------------------------------------------------------------

    ofs << "drawTheta\n";
    {
        for(std::size_t i=0; i<N; ++i)
        {
            const Real D = D_gen(mt);
            const Real k = k_gen(mt);
            const Real a = a_gen(mt);
            boost::random::uniform_real_distribution<Real> s_gen(a * 1e-2, a);
            const Real s = s_gen(mt);
            boost::random::uniform_real_distribution<Real> r0_gen(s, a);
            const Real r0 = r0_gen(mt);

            const greens_functions::GreensFunction2DRadAbs gf(D, k, r0, s, a);

            const Real rnd1  = canonical(mt);
            const Real time  = gf.drawTime(rnd1);
            const Real rnd2  = canonical(mt);
            const Real R     = gf.drawR(rnd2, time);
            const Real rnd3  = canonical(mt);
            const Real Theta = gf.drawTheta(rnd3, R, time);

            ofs << D << ' ' << k << ' ' << a << ' ' << s << ' ' << r0 << ' '
                      << time << ' ' << R << ' ' << Theta << '\n';
        }
    }

    return 0;
}
