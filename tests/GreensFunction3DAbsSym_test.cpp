#define BOOST_TEST_MODULE "GreensFunction3DAbsSym_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <math.h>
#include "../GreensFunction.hpp"
#include "../GreensFunction3DAbsSym.hpp"

using namespace greens_functions;

BOOST_AUTO_TEST_CASE(GF3DAbsSym_constructor)
{
  Real D( 1e-12 ), a( 1.0 );
  GreensFunction3DAbsSym gf( D, a );
}

BOOST_AUTO_TEST_CASE(GF3DAbsSym_no_shell)
{
  Real D( 1e-12 ), a( std::numeric_limits<Real>::infinity() );
  GreensFunction3DAbsSym gf( D, a );
}

BOOST_AUTO_TEST_CASE(GF3DAbsSym_zero_shell)
{
  Real D( 1e-12 ), a( 0.0 );
  GreensFunction3DAbsSym gf( D, a );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK_EQUAL(t, 0.0);
  
  Real r( gf.drawR(0.5, t) );
  BOOST_CHECK_EQUAL(r, 0.0);
}

BOOST_AUTO_TEST_CASE(GF3DAbsSym_DrawTime)
{
  Real D( 1e-12 ), a( 1e-7 );
  GreensFunction3DAbsSym gf( D, a );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK(0.0 < t && t < std::numeric_limits<Real>::infinity());

  t = gf.drawTime( 0.0 );
  BOOST_CHECK(0.0 < t && t < std::numeric_limits<Real>::infinity());

  t = gf.drawTime( 1.0 - 1e-16 );
  BOOST_CHECK(0.0 < t && t < std::numeric_limits<Real>::infinity());
}

BOOST_AUTO_TEST_CASE(GF3DAbsSym_DrawR)
{
  Real D( 1e-12 ), a( 1e-7 );
  GreensFunction3DAbsSym gf( D, a );

  Real t( gf.drawTime( 0.5 ) );

  Real r( gf.drawR( 0.0, t ) );
  r = gf.drawR( 0.5, t );
  r = gf.drawR( 1.0 - 1e-16, t );
}

BOOST_AUTO_TEST_CASE(GF3DAbsSym_DrawR_zerot)
{
  Real D( 1e-12 ), a( 1e-8 ), t( 0.0 );
  GreensFunction3DAbsSym gf( D, a );

  Real r( gf.drawR( 0.5, t ) );
  BOOST_CHECK_EQUAL( 0.0, r );
}

BOOST_AUTO_TEST_CASE(GF3DAbsSym_DrawR_smallt)
{
  Real D( 1e-12 ), a( 1e-8 );
  GreensFunction3DAbsSym gf( D, a );

  Real t( gf.drawTime(0.5) ), r;
  while( t > 1e-30)
  {
    t *= 1e-4;
    r = gf.drawR( 0.0, t );
    BOOST_CHECK(0.0 <= r && r <= a);
  }
}

BOOST_AUTO_TEST_CASE(GF3DAbsSym_DrawR_large_shell)
{
  Real D( 1e-12 ), a( 1e-3 );
  GreensFunction3DAbsSym gf( D, a );

  Real t( 1e-10 );
  Real r( gf.drawR(0.5, t) );
  BOOST_CHECK(0.0 < r && r <= a);
}

BOOST_AUTO_TEST_CASE(GF3DAbsSym_DrawR_large_t)
{
  Real D( 1e-12 ), a( 1e-6 );
  GreensFunction3DAbsSym gf( D, a );

  Real t( gf.drawTime( 0.0 ) );
  Real r( gf.drawR(0.5, t) );
  BOOST_CHECK(0.0 < r && r <= a);
}

BOOST_AUTO_TEST_CASE(GF3DAbsSym_DrawR_regression1)
{
  Real D( 1e-12 ), a( 3.30588e-08 ), rnd( 0.944502 );
  GreensFunction3DAbsSym gf( D, a );

  Real dt( 0.00411832 );
  Real r( gf.drawR( rnd, dt ) );
  BOOST_CHECK(0.0 < r && r <= a);
}

BOOST_AUTO_TEST_CASE(GF3DAbsSym_p_int_r_is_p_int_r_free_with_large_shell)
{
  Real D( 1e-12 ), a( 1e-6 );
  GreensFunction3DAbsSym gf( D, a );

  Real t( 1e-6 );
  Real r( 1e-9 );
  Real p( gf.p_int_r(r, t) );
  Real pfree( gf.p_int_r_free(r, t) );
  BOOST_CHECK_CLOSE( p, pfree, 0.0001);
}

BOOST_AUTO_TEST_CASE(GF3DAbsSym_p_int_r_at_a_is_p_survival)
{
  Real D( 1e-12 ), a( 1e-8 );
  GreensFunction3DAbsSym gf( D, a );

  Real t( 1e-5 );
  Real p( gf.p_int_r(a, t) );
  Real psurv( gf.p_survival( t ) );
  BOOST_CHECK_CLOSE( p, psurv, 0.0001);
}

BOOST_AUTO_TEST_CASE(GF3DAbsSym_p_int_r_at_a_is_p_survival2)
{
  Real D( 1e-12 ), a( 1e-9 );
  GreensFunction3DAbsSym gf( D, a );

  Real t( 1e-3 );
  Real p( gf.p_int_r(a, t) );
  Real psurv( gf.p_survival( t ) );
  BOOST_CHECK_CLOSE( p, psurv, 0.0001);
}


