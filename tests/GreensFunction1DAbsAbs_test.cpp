#define BOOST_TEST_MODULE "GreensFunction1DAbsAbs_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include "../GreensFunction1DAbsAbs.hpp"

using namespace greens_functions;

BOOST_AUTO_TEST_CASE(GF1DAbsAbs_constructor)
{
  Real D( 1e-12 ), v( 0.0 ), L( 2e-7 ), r0( L / 2.0 );
  GreensFunction1DAbsAbs gf( D, v, r0, L );
}

BOOST_AUTO_TEST_CASE(GF1DAbsABs_DrawTime)
{
  Real D( 1e-12 ), v( 0.0 ), L( 2e-7 ), r0( 5e-8 );
  GreensFunction1DAbsAbs gf( D, v, r0, L );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK(0.0 < t < INFINITY);
  std::cout << std::endl;
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawTime : t = " << t << std::endl;

  Real t( gf.drawTime( 0.0 ) );
  BOOST_CHECK(0.0 < t < INFINITY);
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawTime : t = " << t << std::endl;

  Real t( gf.drawTime( 1e-16 ) );
  BOOST_CHECK(0.0 < t < INFINITY);
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawTime : t = " << t << std::endl;

  Real t( gf.drawTime( 1 - 1e-16 ) );
  BOOST_CHECK(0.0 < t < INFINITY);
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawTime : t = " << t << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DAbsAbs_DrawTime_a_equal_sigma)
{
  Real D( 1e-12 ), v( 0.0 ), L( 0.0 ), r0( L );
  GreensFunction1DAbsAbs gf( D, v, r0, L );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK_EQUAL( t, 0.0 );
  std::cout << std::endl;
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawTime_a_equal_sigma : t = " << t << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DAbsAbs_DrawTime_a_near_sigma)
{
  Real D( 1e-12 ), v( 0.0 ), L( 2e-14 ), r0( L / 2 );
  GreensFunction1DAbsAbs gf( D, v, r0, L );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK(0.0 < t < INFINITY);
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawTime_a_near_sigma : t = " << t << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DAbsAbs_DrawTime_r0_equal_a)
{
  Real D( 1e-12 ), v( 0.0 ), L( 2e-7 ), r0( L );
  GreensFunction1DAbsAbs gf( D, v, r0, L );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK_EQUAL( t , 0.0 );
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawTime_r0_equal_a : t = " << t << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DAbsAbs_DrawTime_r0_equal_sigma)
{
  Real D( 1e-12 ), v( 0.0 ), L( 1e-7 ), r0( L );
  GreensFunction1DAbsAbs gf( D, v, r0, L );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK( 0.0 < t < INFINITY );
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawTime_r0_equal_sigma : t = " << t << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DAbsAbs_DrawEventTipe)
{
  Real D( 1e-12 ), v( 0.0 ), L( 2e-7 ), r0( L / 2 );
  GreensFunction1DAbsAbs gf( D, v, r0, L );

  Real t( gf.drawTime( 0.5 ) );
  EventKind eventType(gf.drawEventType(0.5, t));

//  BOOST_CHECK( 0.0 < t < INFINITY );
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawType = " << eventType << std::endl;
}
















