#define BOOST_TEST_MODULE "GreensFunction2DAbs_test"
#define _USE_MATH_DEFINES

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include "../GreensFunction2DAbs.hpp"
#include <iostream>

using namespace greens_functions;

BOOST_AUTO_TEST_CASE(GF2DAbs_constructor)
{
    Real D(1e-12), r0(5e-8), a(1e-7);
    GreensFunction2DAbs gf(D, r0, a);
}

BOOST_AUTO_TEST_CASE(GF2DAbs_DrawTime)
{
    Real D(1e-12), a(1e-7), r0(5e-8);
    GreensFunction2DAbs gf(D, r0, a);
 
    Real t(gf.drawTime(0.5));
    BOOST_CHECK(0.0 < t && t < std::numeric_limits<Real>::infinity());
 
    t = gf.drawTime(0.0);
    BOOST_CHECK_EQUAL(t, 0e0);

    t = gf.drawTime(1.0);
    BOOST_CHECK_EQUAL(t, std::numeric_limits<Real>::infinity());
}

BOOST_AUTO_TEST_CASE(GF2DAbs_DrawTime_r0_equal_a)
{
    Real D(1e-12), a(1e-7), r0(a);
    GreensFunction2DAbs gf(D, r0, a);
 
    Real t(gf.drawTime(0.5));
    BOOST_CHECK_EQUAL(0.0, t);
}

BOOST_AUTO_TEST_CASE(GF2DAbs_DrawR)
{
    Real D(1e-12), a(1e-7), r0(2e-8);
    GreensFunction2DAbs gf(D, r0, a);
 
    Real t(1e-3);
 
    Real r(gf.drawR(0.5, t));
    BOOST_CHECK(0.0 <= r && r <= a);
 
    Real r1(gf.drawR(0.0, t));
    Real r2(gf.drawR(1.0, t));
    BOOST_CHECK(0.0 <= r1 && r1 <= a);
    BOOST_CHECK(0.0 <= r2 && r2 <= a);
 
    BOOST_CHECK(std::abs(r1) <= 1e-15);
    BOOST_CHECK_CLOSE(std::abs(r2 - a), 0.0, 0.0001);
}

BOOST_AUTO_TEST_CASE(GF2DAbs_DrawR_zerot)
{
    Real D(1e-12), a(1e-7), r0(2e-8);
    GreensFunction2DAbs gf(D, r0, a);
 
    Real t(0.0);
 
    Real r(gf.drawR(0.5, t));
    BOOST_CHECK_EQUAL(r0, r);
}

BOOST_AUTO_TEST_CASE(GF2DAbs_DrawR_squeezed)
{
    Real D( 1e-12 ), a( 1.01e-8 ), r0( 1.005e-8 );
    GreensFunction2DAbs gf( D, r0, a );
 
    Real t(1e-6);
 
    Real r(gf.drawR(0.5, t));
    BOOST_CHECK(0.0 <= r && r <= a);
 
//near 0
    r0 = 0.0001e-8;
    GreensFunction2DAbs gf1(D, r0, a);
    r = gf1.drawR(0.5, t);
    BOOST_CHECK(0.0 <= r && r <= a);
 
//near a
    r0 = 1.0099e-8;
    GreensFunction2DAbs gf2(D, r0, a);
    r = gf2.drawR(0.5, t);
    BOOST_CHECK(0.0 <= r && r <= a);
}

BOOST_AUTO_TEST_CASE(GF2DAbs_DrawTheta)
{
    Real D(1e-12), a(1e-7), r0(5e-8);
    GreensFunction2DAbs gf(D, r0, a);
 
    Real t(gf.drawTime(0.5));
    Real r(gf.drawR(0.5, t));
    Real theta;
 
    theta = gf.drawTheta(0.5, r, t);
    BOOST_CHECK( 0.0 <= theta && theta <= 2e0 * M_PI );
 
    theta = gf.drawTheta(0.0, r, t);
    BOOST_CHECK( 0.0 <= theta && theta <= 2e0 * M_PI );
 
    theta = gf.drawTheta(1.0, r, t);
    BOOST_CHECK( 0.0 <= theta && theta <= 2e0 * M_PI );
}

BOOST_AUTO_TEST_CASE(GF2DAbs_DrawTheta_zerot)
{
    Real D( 1e-12 ), a( 1e-7 ), r0( 5e-8 );
    GreensFunction2DAbs gf( D, r0, a );
 
    Real t( 0.0 );
    Real theta(gf.drawTheta(0.5, r0, t));
    BOOST_CHECK_EQUAL( 0.0, theta );
}

BOOST_AUTO_TEST_CASE(GF2DAbs_DrawTheta_smallt)
{
    Real D( 1e-12 ), a( 1e-7 ), r0( 5e-8 );
    GreensFunction2DAbs gf( D, r0, a );
 
    Real r( 5e-8 );
    Real t( 1e-4 );
    Real theta(gf.drawTheta(0.5, r, t));
    BOOST_CHECK( 0.0 <= theta && theta <= 2e0 * M_PI );
}

BOOST_AUTO_TEST_CASE(GF2DAbs_DrawTheta_larget)
{
    Real D( 1e-12 ), a( 1e-7 ), r0( 5e-8 );
    GreensFunction2DAbs gf( D, r0, a );
 
    Real r(5e-8);
    Real t(0.15);
    for(int i(0); i<15; ++i)
    {
        Real theta(gf.drawTheta(0.5, r, t));
        BOOST_CHECK(0.0 <= theta && theta <= 2e0 * M_PI);
        t += 0.01;
    }
}

BOOST_AUTO_TEST_CASE(GF2DAbs_DrawTheta_near_a)
{
    Real D( 1e-12 ), a( 1.01e-5 ), r0( 1.0099e-5 );
    GreensFunction2DAbs gf( D, r0, a );
 
    Real t( 1e-5 );
    Real theta(gf.drawTheta(0.5, r0, t));
    BOOST_CHECK( 0.0 <= theta && theta <= 2e0 * M_PI );
}

/* When r = a, p_int_theta is always zero. So drawTheta is not defined. */
// BOOST_AUTO_TEST_CASE(GF2DAbs_DrawTheta_r_equal_a)
// {
//   Real D( 1e-12 ), a( 1e-7 ), r0( 9e-8 ), r( a );
//   GreensFunction2DAbs gf( D, r0, a );
//
//   Real t( 1e-4 );
//   Real theta(gf.drawTheta(0.5, r, t));
//   BOOST_CHECK( 0.0 <= theta && theta <= 2e0 * M_PI );
// }

// BOOST_AUTO_TEST_CASE(GF2DAbs_DrawTheta_r0_near_a_r_equal_a)
// {
//   Real D( 1e-12 ), a( 1e-7 ), r0( a - 1e-9 ), r( a );
//   GreensFunction2DAbs gf( D, r0, a );
//
//   Real t( 1e-2 );
//   Real theta(gf.drawTheta(0.5, r, t));
//   BOOST_CHECK( 0.0 <= theta && theta <= M_PI );
// }

// ************************** p_survival ****************************//
BOOST_AUTO_TEST_CASE(GF2DAbs_p_survival)
{
    Real D(1e-12), a(1e-7), r0(5e-8);
    GreensFunction2DAbs gf(D, r0, a);
 
    Real t(1e-3);
    Real pintr(gf.p_survival(t));
    BOOST_CHECK(0.0 <= pintr && pintr <= 1.0);

    t = 1e5;
    pintr = gf.p_survival(t);
    BOOST_CHECK(0.0 <= pintr && pintr <= 1.0);

    t = 0e0;
    pintr = gf.p_survival(t);
    BOOST_CHECK_EQUAL(pintr, 1e0);
}

BOOST_AUTO_TEST_CASE(GF2DAbs_p_survival_never_negative)
{
    Real D(1e-12), a(1e-7), r0(5e-8);
    GreensFunction2DAbs gf(D, r0, a);
 
    Real t(1e0);
    Real psurv(gf.p_survival(0e0));
    int resolution(20);

    for(int i(0); i<resolution; ++i)
    {
        t *= 2e0;
        psurv = gf.p_survival(t);
        BOOST_CHECK(0e0 <= psurv);
    }
}

BOOST_AUTO_TEST_CASE(GF2DAbs_p_survival_never_increase)
{
    Real D(1e-12), a(1e-7), r0(5e-8);
    GreensFunction2DAbs gf(D, r0, a);
 
    Real t(1e0);
    Real psurv(gf.p_survival(0e0));
    Real psurv_before(psurv);
    int resolution(20);

    for(int i(0); i<resolution; ++i)
    {
        t *= 2e0;
        psurv = gf.p_survival(t);
        BOOST_CHECK(psurv <= psurv_before);
        psurv_before = psurv;
    }
}

// ************************** p_int_r ****************************//
BOOST_AUTO_TEST_CASE(GF2DAbs_p_int_r)
{
    Real D( 1e-12 ), a( 1e-7 ), r0( 5e-8 ), r( r0 );
    GreensFunction2DAbs gf( D, r0, a );
 
    Real t( 1e-3 );
    Real pintr( gf.p_int_r( r, t ) );
    BOOST_CHECK( 0.0 <= pintr && pintr <= 1.0 );
}

/* p_int_r(a, t) = p_survival(t) */
BOOST_AUTO_TEST_CASE(GF2DAbs_p_int_r_at_a_is_p_survival)
{
    Real D( 1e-12 ), a( 1e-7 ), r0( 5e-8 );
    GreensFunction2DAbs gf( D, r0, a );
 
    Real t( 1e-3 );
    Real pintr( gf.p_int_r( a, t ) );
    Real psurv( gf.p_survival( t ) );
 
    BOOST_CHECK_CLOSE( pintr, psurv, 0.0001 );
}

BOOST_AUTO_TEST_CASE(GF2DAbs_p_int_r_at_zero_is_zero)
{
    Real D( 1e-12 ), a( 1e-7 ), r0( 5e-8 );
    GreensFunction2DAbs gf( D, r0, a );
 
    Real t( 1e-3 );
    Real pintr( gf.p_int_r( 0.0, t ) );
 
    BOOST_CHECK_EQUAL( pintr, 0.0 );
}


BOOST_AUTO_TEST_CASE(GF2DAbs_p_int_r_never_decrease)
{
    Real D( 1e-12 ), a( 1e-7 ), r0( 5e-8 );
    GreensFunction2DAbs gf( D, r0, a );
 
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

BOOST_AUTO_TEST_CASE(GF2DAbs_p_int_r_never_negative)
{
    Real D( 1e-12 ), a( 1e-7 ), r0( 5e-8 );
    GreensFunction2DAbs gf( D, r0, a );
 
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

BOOST_AUTO_TEST_CASE(GF2DAbs_p_int_theta)
{
    Real D(1e-12), a(1e-7), r0(5e-8), r(r0);
    GreensFunction2DAbs gf(D, r0, a);

    Real t(1e-3), theta(5e-1);
    Real pinttheta = gf.p_int_theta(r, theta, t);
    Real pint2pi = gf.p_int_2pi(r, t);
    BOOST_CHECK(0e0 <= pinttheta && pinttheta <= pint2pi);

    pinttheta = gf.p_int_theta(r, 0e0, t);
    BOOST_CHECK_EQUAL(pinttheta, 0e0);

    pinttheta = gf.p_int_theta(r, 2e0 * M_PI, t);
    BOOST_CHECK_EQUAL(pinttheta, pint2pi);

    pinttheta = gf.p_int_theta(a, theta, t);
    BOOST_CHECK_EQUAL(pinttheta, 0e0);
}

BOOST_AUTO_TEST_CASE(GF2DAbs_p_int_theta_never_negative)
{
    Real D( 1e-12 ), a( 1e-7 ), r0( 5e-8 ), r( r0 );
    GreensFunction2DAbs gf( D, r0, a );
 
    Real t( 1e-3 );
    Real pmin(1e0), p, theta;
    int resolution(50);
 
    for(int i(0); i<resolution; i++)
    {
        theta = i * M_PI / resolution;
        p = gf.p_int_theta(r, theta, t);
        pmin = std::min(pmin, p);
    }
 
    BOOST_CHECK(0.0 <= pmin);
}

BOOST_AUTO_TEST_CASE(GF2DAbs_p_theta_never_decrease)
{
    Real D( 1e-12 ), a( 1e-7 ), r0( 5e-8 ), r( r0 );
    GreensFunction2DAbs gf( D, r0, a );
 
    Real t( 1e-3 );
    Real pint_prev(0e0), pint(0e0), theta;
    int resolution(50);
 
    for(int i(0); i<resolution; i++)
    {
        theta = i * M_PI / resolution;
        pint = gf.p_int_theta(r, theta, t);
        BOOST_CHECK(pint >= pint_prev);
        pint_prev = pint;
    }
}
