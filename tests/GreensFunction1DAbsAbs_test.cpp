#define BOOST_TEST_MODULE "GreensFunction1DAbsAbs_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include "../GreensFunction.hpp"
#include "../GreensFunction1DAbsAbs.hpp"

using namespace greens_functions;


BOOST_AUTO_TEST_CASE(GF1DAbsAbs_constructor1)
{
  Real D( 1e-12 ), r0( 1e-7 ), sigma( 1e-7 ), a( 2e-7 );
  GreensFunction1DAbsAbs gf( D, r0, sigma, a );
}

BOOST_AUTO_TEST_CASE(GF1DAbsAbs_constructor2)
{
  Real D( 1e-12 ), v(0.0), r0( 1e-7 ), sigma( 1e-7 ), a( 2e-7 );
  GreensFunction1DAbsAbs gf( D, v, r0, sigma, a );
}

BOOST_AUTO_TEST_CASE(GF1DAbsABs_DrawTime)
{
  Real D( 1e-12 ), v(0.0), r0( 5e-8 ), sigma(2e-7), a(4e-7); //L = a - sigma = 2e-7
  GreensFunction1DAbsAbs gf( D, v, r0, sigma, a );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK(0.0 < t < INFINITY);
  std::cout << std::endl;
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawTime : t = " << t << std::endl;

  t = gf.drawTime( 0.0 );
  BOOST_CHECK(0.0 < t < INFINITY);
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawTime : t = " << t << std::endl;

  t = gf.drawTime( 1e-16 ) ;
  BOOST_CHECK(0.0 < t < INFINITY);
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawTime : t = " << t << std::endl;

  t = gf.drawTime( 1 - 1e-16 );
  BOOST_CHECK(0.0 < t < INFINITY);
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawTime : t = " << t << std::endl;
  std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DAbsAbs_DrawTime_a_equal_sigma)
{
  Real D( 1e-12 ), v( 0.0 ), a( 1e-7 ), r0( 0.0 ), sigma( 1e-7 );
  GreensFunction1DAbsAbs gf( D, v, r0, sigma, a );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK_EQUAL( t, 0.0 );
  std::cout << std::endl;
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawTime_a_equal_sigma : t = " << t << std::endl;
  std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DAbsAbs_DrawTime_a_near_sigma)
{
  Real D( 1e-12 ), v( 0.0 ), r0( 1e-14 ), a( 4e-14 ), sigma( 2e-14 );//L( 2e-14 )
  GreensFunction1DAbsAbs gf( D, v, r0, sigma, a );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK(0.0 < t < INFINITY);
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawTime_a_near_sigma : t = " << t << std::endl;
  std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DAbsAbs_DrawTime_r0_equal_a)
{
  Real D( 1e-12 ), v( 0.0 ), r0( 2e-7 ), sigma( 0.0 ), a( 2e-7 );//L( 2e-7 ), r0( L )
  GreensFunction1DAbsAbs gf( D, v, r0, sigma, a );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK_EQUAL( t , 0.0 );
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawTime_r0_equal_a : t = " << t << std::endl;
  std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DAbsAbs_DrawTime_r0_equal_sigma)
{
  Real D( 1e-12 ), v( 0.0 ), r0( 1e-7 ), a( 2e-7 ), sigma( 1e-7 ); //L(1e-7)
  GreensFunction1DAbsAbs gf( D, v, r0, sigma, a );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK( 0.0 < t < INFINITY );
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawTime_r0_equal_sigma : t = " << t << std::endl;
  std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DAbsAbs_DrawEventType)
{
  Real D( 1e-12 ), v( 0.0 ), r0( 1e-7 ), sigma( 1e-8 ), a( 2.1e-7 ); //L(2e-7)
  GreensFunction1DAbsAbs gf( D, v, r0, sigma, a );

  Real t( gf.drawTime( 0.5 ) );
//  std::cout << "DrawEventType gf.drawTime(0.5): " << t << std::endl;
  EventKind eventType( gf.drawEventType(0.5, t) );

  BOOST_CHECK( ( eventType == IV_ESCAPE || eventType == IV_REACTION ) );
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawType(0.5, t) = " << eventType << std::endl;

//when r0 = 1e-7, sigma = 2e-7, a = 4e-7, fluxratio > 1 and error detected.
  eventType = gf.drawEventType(0.999999, t);
  BOOST_CHECK_EQUAL(eventType, IV_ESCAPE); //ESCAPE
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawType(0.999999,t) = " << eventType << std::endl;

  eventType = gf.drawEventType(0.0, t);
  BOOST_CHECK_EQUAL(eventType, IV_REACTION); //REACTION
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawType(0,t) = " << eventType << std::endl;
  std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DAbsAbs_DrawEventType_smallt)
{
  Real D( 1e-12 ), v( 0.0 ), r0( 1e-6 ), sigma(1e-7), a(2.1e-6); //L(2e-6), r0(L/2)
  GreensFunction1DAbsAbs gf( D, v, r0, sigma, a );

  Real t( gf.drawTime( 0.001 ) );
  EventKind eventType(gf.drawEventType(0.5, t));

  BOOST_CHECK( eventType == IV_ESCAPE || eventType == IV_REACTION );
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawType = " << eventType << std::endl;

  eventType = gf.drawEventType(0.999999, t);
  BOOST_CHECK_EQUAL(eventType, IV_ESCAPE); //ESCAPE
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawType(1,t) = " << eventType << std::endl;

  eventType = gf.drawEventType(0.0, t);
  BOOST_CHECK_EQUAL(eventType, IV_REACTION); //REACTION
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawType(0,t) = " << eventType << std::endl;
  std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DAbsAbs_DrawR)
{
  Real D( 1e-12 ), v( 0.0 ), r0( 1e-7 ), sigma(0.0), a(2e-7);//L(2e-7), r0(L/2)
  Real L( a - sigma );
  GreensFunction1DAbsAbs gf( D, v, r0, sigma, a );

  Real t( 1e-3 );
  Real r(gf.drawR(0.5, t));
  BOOST_CHECK( 0 <= r && r <= L);
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawR : r = " << r << std::endl;

  Real r1 ( gf.drawR(0.0, t) );
  Real r2 ( gf.drawR(0.999999999999, t));
  BOOST_CHECK( r1 == 0 ); //(sigma != 0) => failuar
  BOOST_CHECK( 0 <= r2 && r2 <= L);

  std::cout << "GreensFunction1DAbsAbs_test : test_DrawR : r1 = " << r1 << std::endl;
  std::cout << "GreensFunction1DAbsAbs_test : test_DrawR : r2 = " << r2 << std::endl;
  std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DAbsAbs_DrawR_zerot)
{
  Real D( 1e-12 ), v( 0.0 ), r0( 0.5e-7 ), sigma( 1e-8 ), a( 2.1e-7 ); //L( 1e-7 ), r0(L/2)
  GreensFunction1DAbsAbs gf( D, v, r0, sigma, a );

  Real t( 0.0 );
  Real r(gf.drawR(0.5, t));

  BOOST_CHECK_EQUAL( r0, r );

  std::cout << "GreensFunction1DAbsAbs_test : test_DrawR_zerot : r = " << r << std::endl;
  std::cout << std::endl;
}


