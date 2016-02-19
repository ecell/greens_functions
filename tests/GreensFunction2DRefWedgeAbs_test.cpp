#define BOOST_TEST_MODULE "GreensFunction2DRefWedgeAbs_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include "../GreensFunction2DRefWedgeAbs.hpp"
#include <iostream> 

using namespace greens_functions;

// constructable
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_constructor)
{
    Real D(1e-12), r0(5e-8), a(1e-7), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
}

// Draw time returns positive value and the range is from 0 to infty
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_DrawTime)
{
    Real D(1e-12), a(1e-7), r0(5e-8), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
 
    int resolution = 50;
    Real t, drandom(1e0 / resolution), randomnum(drandom);
    for(int i=1; i < resolution; ++i)
    {
        t = gf.drawTime(randomnum);
        BOOST_CHECK(0.0 < t && t < INFINITY);
        randomnum += drandom;
    }

    // boundary value
    t = gf.drawTime(/*rnd = */0.0);
    BOOST_CHECK_EQUAL(t, 0e0);

    t = gf.drawTime(/*rnd = */1.0);
    BOOST_CHECK_EQUAL(t, INFINITY);
}

//if initial r is on abs boundary, particle doesn't move.
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_DrawTime_r0_equal_a)
{
    Real D(1e-12), a(1e-7), r0(a), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);

    int resolution = 50;
    Real t, drandom(1e0 / resolution), randomnum(0e0);
    for(int i=0; i < resolution; ++i)
    {
        t = gf.drawTime(randomnum);
        BOOST_CHECK_EQUAL(0.0, t);
        randomnum += drandom;
    }
}

// return value of drawR is in range 0 to a
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_DrawR)
{
    Real D(1e-12), a(1e-7), r0(2e-8), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
 
    int resolution = 50;
    Real t(1e-3), drandom(1e0 / resolution), randomnum(drandom);
    Real r;
    for(int i=1; i < resolution; ++i)
    {
        r = gf.drawR(randomnum, t);
        BOOST_CHECK(0.0 <= r && r <= a);
        randomnum += drandom;
    }
 
    Real r1(gf.drawR(/*rnd = */0.0, t));
    Real r2(gf.drawR(/*rnd = */1.0, t));
    BOOST_CHECK(0.0 <= r1 && r1 <= a);
    BOOST_CHECK(0.0 <= r2 && r2 <= a);
 
    BOOST_CHECK_CLOSE(r1, 0e0, 0.0001);
    BOOST_CHECK_CLOSE(r2,  a,  0.0001);
}

// after zero sec, particle is on the initial position.
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_DrawR_zerot)
{
    Real D(1e-12), a(1e-7), r0(2e-8), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
 
    int resolution = 50;
    Real t(0e0), drandom(1e0 / resolution), randomnum(drandom);
    Real r;
    for(int i=1; i < resolution; ++i)
    {
        r = gf.drawR(randomnum, t);
        BOOST_CHECK_EQUAL(r0, r);
        randomnum += drandom;
    }
}

// if initial r is near zero or a, drawR return value being positive and lesser than a
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_DrawR_squeezed)
{
    // r0 is near a
    Real D(1e-12), a(1.01e-8), r0(1.005e-8), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
 
    Real t(1e-6);
 
    Real r(gf.drawR(/*rnd = */0.5, t));
    BOOST_CHECK(0.0 <= r && r <= a);
 
//near 0
    r0 = 0.0001e-8;
    GreensFunction2DRefWedgeAbs gf1(D, r0, a, phi);
    r = gf1.drawR(0.5, t);
    BOOST_CHECK(0.0 <= r && r <= a);
 
//near a
    r0 = 1.0099e-8;
    GreensFunction2DRefWedgeAbs gf2(D, r0, a, phi);
    r = gf2.drawR(0.5, t);
    BOOST_CHECK(0.0 <= r && r <= a);
}

// ******************************* Draw Theta *********************************

// draw theta return value that is in range 0 to phi.
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_DrawTheta)
{
    Real D(1e-12), a(1e-7), r0(5e-8), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
 
    Real t(gf.drawTime(0.5));
    Real r(gf.drawR(0.5, t));
    Real theta;
 
    theta = gf.drawTheta(0.5, r, t);
    BOOST_CHECK(0.0 <= theta && theta <= phi);
 
    theta = gf.drawTheta(0.0, r, t);
    BOOST_CHECK(0.0 <= theta && theta <= phi);
    BOOST_CHECK_CLOSE(theta, 0e0, 0.0001);
 
    theta = gf.drawTheta(1.0, r, t);
    BOOST_CHECK(0.0 <= theta && theta <= phi);
    BOOST_CHECK_CLOSE(theta, phi, 0.0001);
}

// test drawtheta for some phi
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_DrawTheta_for_Phi)
{
    Real D(1e-12), a(1e-7), r0(5e-8);

    int resolution = 50;
    Real phi, r, t;
    Real dtheta = M_PI * 2e0 / resolution;

    for(int i=1; i < resolution; ++i)
    {
        for(int j=0; j<10; ++j)
        {
            phi = i * dtheta;
            GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
            t = gf.drawTime(0.5);
            r = gf.drawR(0.5, t);

            Real theta = gf.drawTheta(/*rand = */ 0.5, r, t);
            BOOST_CHECK(0.0 <= theta && theta <= phi);
        }
    }
}

BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_DrawTheta_extreme)
{
    Real D(1e0), a(1e0), r0(5e-1);

    Real phi, r, t, theta;

    // ordinary r, small phi
    phi = 1e-9;
    r = 5e-1;
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
    t = gf.drawTime(0.5);
    theta = gf.drawTheta(/*rand = */ 0.5, r, t);
    BOOST_CHECK(0.0 <= theta && theta <= phi);

    phi = 1.999999 * M_PI;
    r = 5e-1;
    GreensFunction2DRefWedgeAbs gf1(D, r0, a, phi);
    t = gf1.drawTime(0.5);
    theta = gf1.drawTheta(/*rand = */0.5, r, t);
    BOOST_CHECK(0.0 <= theta && theta <= phi);

    // small r
    phi = 1.999999 * M_PI;
    r = 1e-10;
    GreensFunction2DRefWedgeAbs gf2(D, r0, a, phi);
    t = gf2.drawTime(0.5);
    theta = gf2.drawTheta(/*rand = */ 0.5, r, t);
    BOOST_CHECK(0.0 <= theta && theta <= phi);

    phi = 1e-9;
    r = 1e-9;
    GreensFunction2DRefWedgeAbs gf3(D, r0, a, phi);
    t = gf3.drawTime(0.5);
    theta = gf3.drawTheta(/*rand = */0.5, r, t);
    BOOST_CHECK(0.0 <= theta && theta <= phi);

    // r ~ a
    phi = 1.999999 * M_PI;
    r = a - 1e-9;
    GreensFunction2DRefWedgeAbs gf4(D, r0, a, phi);
    t = gf4.drawTime(0.5);
    theta = gf4.drawTheta(/*rand = */ 0.5, r, t);
    BOOST_CHECK(0.0 <= theta && theta <= phi);

    phi = 1e-9;
    r = a - 1e-9;
    GreensFunction2DRefWedgeAbs gf5(D, r0, a, phi);
    t = gf5.drawTime(0.5);
    theta = gf5.drawTheta(/*rand = */0.5, r, t);
    BOOST_CHECK(0.0 <= theta && theta <= phi);

// extreme t
{
    phi = 1e-9;
    r = 1e-9;
    GreensFunction2DRefWedgeAbs gf6(D, r0, a, phi);
    t = 1e-9;
    theta = gf6.drawTheta(/*rand = */ 0.5, r, t);
    BOOST_CHECK(0.0 <= theta && theta <= phi);

    phi = 1e-9;
    r = 1e-9;
    GreensFunction2DRefWedgeAbs gf7(D, r0, a, phi);
    t = 1e9;
    theta = gf7.drawTheta(/*rand = */0.5, r, t);
    BOOST_CHECK(0.0 <= theta && theta <= phi);

    phi = 1.999999 * M_PI;
    r = 1e-9;
    GreensFunction2DRefWedgeAbs gf8(D, r0, a, phi);
    t = 1e-9;
    theta = gf8.drawTheta(/*rand = */ 0.5, r, t);
    BOOST_CHECK(0.0 <= theta && theta <= phi);

    phi = 1.999999 * M_PI;
    r = 1e-9;
    GreensFunction2DRefWedgeAbs gf9(D, r0, a, phi);
    t = 1e9;
    theta = gf9.drawTheta(/*rand = */0.5, r, t);
    BOOST_CHECK(0.0 <= theta && theta <= phi);
}

{
    phi = 1e-9;
    r = a - 1e-9;
    GreensFunction2DRefWedgeAbs gf6(D, r0, a, phi);
    t = 1e-9;
    theta = gf6.drawTheta(/*rand = */ 0.5, r, t);
    BOOST_CHECK(0.0 <= theta && theta <= phi);

    phi = 1e-9;
    r = a - 1e-9;
    GreensFunction2DRefWedgeAbs gf7(D, r0, a, phi);
    t = 1e9;
    theta = gf7.drawTheta(/*rand = */0.5, r, t);
    BOOST_CHECK(0.0 <= theta && theta <= phi);

    phi = 1.999999 * M_PI;
    r = a - 1e-9;
    GreensFunction2DRefWedgeAbs gf8(D, r0, a, phi);
    t = 1e-9;
    theta = gf8.drawTheta(/*rand = */ 0.5, r, t);
    BOOST_CHECK(0.0 <= theta && theta <= phi);

    phi = 1.999999 * M_PI;
    r = a - 1e-9;
    GreensFunction2DRefWedgeAbs gf9(D, r0, a, phi);
    t = 1e9;
    theta = gf9.drawTheta(/*rand = */0.5, r, t);
    BOOST_CHECK(0.0 <= theta && theta <= phi);
}

}

// if t == zero, drawtheta returns zero.
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_DrawTheta_zerot)
{
    Real D(1e-12), a(1e-7), r0(5e-8), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
 
    Real t(0.0);
    Real theta(gf.drawTheta(0.5, r0, t));
    BOOST_CHECK_EQUAL(0.0, theta);
}

// if t is small
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_DrawTheta_smallt)
{
    Real D(1e-12), a(1e-7), r0(5e-8), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
 
    Real r(5e-8);
    Real t(1e-4);
    Real theta(gf.drawTheta(0.5, r, t));
    BOOST_CHECK(0.0 <= theta && theta <= phi);
}

// if t is large
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_DrawTheta_larget)
{
    Real D(1e-12), a(1e-7), r0(5e-8), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
 
    Real r(5e-8);
    Real t(0.15);
    for(int i(0); i<15; ++i)
    {
        Real theta(gf.drawTheta(0.5, r, t));
        BOOST_CHECK(0.0 <= theta && theta <= phi);
        t += 0.01;
    }
}

// if initial r is near a drawtheta returns value in range(0, phi)
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_DrawTheta_near_a)
{
    Real D(1e-12), a(1.01e-5), r0(1.0099e-5), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
 
    Real t(1e-5);
    Real theta(gf.drawTheta(0.5, r0, t));
    BOOST_CHECK(0.0 <= theta && theta <= phi);
}

/* When r = a, p_int_theta is always zero. So drawTheta is not defined. */
// BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_DrawTheta_r_equal_a)
// {
//   Real D( 1e-12 ), a( 1e-7 ), r0( 9e-8 ), r( a );
//   GreensFunction2DRefWedgeAbs gf( D, r0, a );
//
//   Real t( 1e-4 );
//   Real theta(gf.drawTheta(0.5, r, t));
//   BOOST_CHECK( 0.0 <= theta && theta <= 2e0 * M_PI );
// }

// BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_DrawTheta_r0_near_a_r_equal_a)
// {
//   Real D( 1e-12 ), a( 1e-7 ), r0( a - 1e-9 ), r( a );
//   GreensFunction2DRefWedgeAbs gf( D, r0, a );
//
//   Real t( 1e-2 );
//   Real theta(gf.drawTheta(0.5, r, t));
//   BOOST_CHECK( 0.0 <= theta && theta <= M_PI );
// }

// ************************** p_survival ****************************//
// p_survival returns value in range 0 to 1
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_p_survival)
{
    Real D(1e-12), a(1e-7), r0(5e-8), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
 
    Real t(1e-3);
    Real pintr(gf.p_survival(t));
    BOOST_CHECK(0.0 <= pintr && pintr <= 1.0);

    // with large t
    t = 1e5;
    pintr = gf.p_survival(t);
    BOOST_CHECK(0.0 <= pintr && pintr <= 1.0);

    //with zero t
    t = 0e0;
    pintr = gf.p_survival(t);
    BOOST_CHECK_EQUAL(pintr, 1e0);
}

// p_survival returns value being positive for any t
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_p_survival_never_negative)
{
    Real D(1e-12), a(1e-7), r0(5e-8), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
 
    Real t(1e0);
    Real psurv(gf.p_survival(0e0));
    int resolution(20);

    // for t = 2^n for n = 1 to 20
    for(int i(0); i<resolution; ++i)
    {
        t *= 2e0;
        psurv = gf.p_survival(t);
        BOOST_CHECK(0e0 <= psurv);
    }
}

// p_survival always decrease when t inclease
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_p_survival_never_increase)
{
    Real D(1e-12), a(1e-7), r0(5e-8), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
 
    Real t(1e0);
    Real psurv(gf.p_survival(0e0));
    Real psurv_before(psurv);
    int resolution(50);

    for(int i(0); i<resolution; ++i)
    {
        t *= 1.1;
        psurv = gf.p_survival(t);
        BOOST_CHECK(psurv <= psurv_before);
        psurv_before = psurv;
    }
}

// ************************** p_int_r ****************************//
// p_int_r returns value in range 0 to 1
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_p_int_r)
{
    Real D(1e-12), a(1e-7), r0(5e-8), r(r0), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
 
    Real t(1e-3);
    Real pintr(gf.p_int_r(r, t));
    BOOST_CHECK(0.0 <= pintr && pintr <= 1.0);
}

// assert p_int_r(a, t) = p_survival(t)
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_p_int_r_at_a_is_p_survival)
{
    Real D(1e-12), a(1e-7), r0(5e-8), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
 
    Real t(1e-3);
    Real pintr(gf.p_int_r(a, t));
    Real psurv(gf.p_survival(t));
    BOOST_CHECK_CLOSE(pintr, psurv, 0.0001);

    int resolution = 50;
    for(int i=0; i<resolution; ++i)
    {
        t *= 1.1;
        pintr = gf.p_int_r(a, t);
        psurv = gf.p_survival(t);
        BOOST_CHECK_CLOSE(pintr, psurv, 0.0001);
    }
}

// where r = 0, p_int_r returns 0
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_p_int_r_at_zero_is_zero)
{
    Real D(1e-12), a(1e-7), r0(5e-8), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
 
    Real t(1e-3);
    Real pintr(gf.p_int_r(/* r = */0.0, t));
 
    BOOST_CHECK_EQUAL(pintr, 0.0);
}

// p_int_r is monotonic increase
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_p_int_r_never_decrease)
{
    Real D(1e-12), a(1e-7), r0(5e-8), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
 
    Real t(1e-3), r(0e0), pintr;
    Real pintr_before(gf.p_int_r(0, t));
    int resolution(50);

    for(int i(1); i<resolution; ++i)
    {
        r = i * a / resolution;
        pintr = gf.p_int_r(r, t);
        BOOST_CHECK(pintr_before <= pintr);
        pintr_before = pintr;
    }
}

// p_int_r always returns positive value
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_p_int_r_never_negative)
{
    Real D(1e-12), a(1e-7), r0(5e-8), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
 
    Real t(1e-3), r(0e0), pintr;
    int resolution(50);

    for(int i(0); i<resolution; ++i)
    {
        r = i * a / resolution;
        pintr = gf.p_int_r(r, t);
        BOOST_CHECK(0e0 <= pintr);
    }
}

// ************************** p_int_theta ****************************//

// p_int_theta returns value in range 0 to phi
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_p_int_theta)
{
    Real D(1e-12), a(1e-7), r0(5e-8), r(r0), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);

    Real t(1e-3), theta(5e-1);
    Real pinttheta = gf.p_int_theta(r, theta, t);
    Real pintphi = gf.p_int_phi(r, t);
    BOOST_CHECK(0e0 <= pinttheta && pinttheta <= pintphi);

    // if theta == 0, return 0e0
    pinttheta = gf.p_int_theta(r, 0e0, t);
    BOOST_CHECK_EQUAL(pinttheta, 0e0);

    // if theta == phi, return p_int_phi.
    pinttheta = gf.p_int_theta(r, phi, t);
    BOOST_CHECK_EQUAL(pinttheta, pintphi);

//     // where r = a, return 0e0??
//     pinttheta = gf.p_int_theta(a, theta, t);
//     BOOST_CHECK_EQUAL(pinttheta, 0e0);
}

BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_p_int_theta_extreme)
{
    Real D(1e0), a(1e0), r0(0.5), phi(1e0);

    for(int i = 1; i < 10; ++i)
    {
        std::cout << "phi = " << phi << std::endl;
        GreensFunction2DRefWedgeAbs gf(D, r0, a, phi * M_PI);
        Real t = gf.drawTime(0.5);
        Real r = gf.drawR(0.5, t);
        for(int j=1; j < 10; ++j)
        {
            Real theta = phi * M_PI * 0.1 * j;
            std::cout << "theta = " << theta << std::endl;
            Real pinttheta = gf.p_int_theta(r, theta, t);
            Real pintphi = gf.p_int_phi(r, t);
            BOOST_CHECK(0e0 < pinttheta && pinttheta < pintphi);
        }
        phi *= 1e-1;
    }
}


// never be negative
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_p_int_theta_never_negative)
{
    Real D(1e-12), a(1e-7), r0(5e-8), r(r0), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
 
    Real t( 1e-3 );
    Real p, theta;
    Real pintphi = gf.p_int_phi(r, t);
    int resolution(50);
 
    for(int i(0); i<resolution; i++)
    {
        theta = i * phi / resolution;
        p = gf.p_int_theta(r, theta, t);
        BOOST_CHECK(0.0 <= p && p <= pintphi);
    }
}

// monomeric inclease
BOOST_AUTO_TEST_CASE(GF2DRefWedgeAbs_p_theta_never_decrease)
{
    Real D(1e-12), a(1e-7), r0(5e-8), r(r0), phi(M_PI);
    GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);
 
    Real t(1e-3);
    Real pint_prev(0e0), pint(0e0), theta;
    int resolution(50);
 
    for(int i(0); i<resolution; i++)
    {
        theta = i * phi / resolution;
        pint = gf.p_int_theta(r, theta, t);
        BOOST_CHECK(pint >= pint_prev);
        pint_prev = pint;
    }
}
