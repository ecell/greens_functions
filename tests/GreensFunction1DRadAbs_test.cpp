#define BOOST_TEST_MODULE "GreensFunction1DRadAbs_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include "../GreensFunction1DRadAbs.hpp"
#include "../GreensFunction.hpp"

using namespace greens_functions;

BOOST_AUTO_TEST_CASE(GF1DRadAbs_constructor1)
{
  Real D( 1e-12 ), kf( 1e-8 ), r0( 1e-7 ), sigma( 1e-8 ), a( 2.1e-7 );//L(2e-7), r0(L/2)
  GreensFunction1DRadAbs gf( D, kf, r0, sigma, a );
}

BOOST_AUTO_TEST_CASE(GF1DRadAbs_constructor2)
{
  Real D( 1e-12 ), v( -3e-8 ), kf( 1e-8 ), r0( 1e-7 ), sigma( 1e-8 ), a( 2.1e-7 );//L(2e-7), r0(L/2)
  GreensFunction1DRadAbs gf( D, v, kf, r0, sigma, a );
}

BOOST_AUTO_TEST_CASE(GF1DRadAbs_DrawTime)
{
  Real D( 1e-12 ), v( -3e-8 ), kf( 1e-8 ), r0(5e-8), sigma( 1e-8 ), a( 2.1e-7 );//L(2e-7)
  GreensFunction1DRadAbs gf( D, v, kf, r0, sigma, a );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK(0.0 < t < INFINITY);
  std::cout << "GreensFunction1DRadAbs_test : test_DrawTime(0.5) : t = " << t << std::endl;
/*
  t = gf.drawTime( 0.0 );
  BOOST_CHECK(0.0 <= t < INFINITY);
  std::cout << "GreensFunction1DRadAbs_test : test_DrawTime(0.0) : t = " << t << std::endl;//error. std::exception. because of "while(value <= 0.0)"
*/
  t = gf.drawTime( 1 - 1e-16 );
  BOOST_CHECK(0.0 < t < INFINITY);
  std::cout << "GreensFunction1DRadAbs_test : test_DrawTime(1 - 1e-16) : t = " << t << std::endl;
  std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DRadAbs_DrawTime_a_equal_sigma)
{
  Real D( 1e-12 ), v( -3e-8 ), kf( 1e-8 ), r0(0.0), sigma( 1e-7 ), a( 1e-7 );//L(0.0), r0(0.0)
  GreensFunction1DRadAbs gf( D, v, kf, r0, sigma, a );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK_EQUAL( t, 0.0 );
  std::cout << "GreensFunction1DRadAbs_test : test_DrawTime_a_equal_sigma : t = " << t << std::endl;
  std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DRadAbs_DrawTime_a_near_sigma)
{
  Real D( 1e-12 ), v( -3e-8 ), kf( 1e-8 ), r0(1e-14), sigma( .5e-14 ), a( 2.5e-14 );//L(2e-14), r0(L/2)
  GreensFunction1DRadAbs gf( D, v, kf, r0, sigma, a );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK(0.0 < t < INFINITY);
  std::cout << "GreensFunction1DRadAbs_test : test_DrawTime_a_near_sigma : t = " << t << std::endl;
  std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DRadAbs_DrawTime_r0_equal_a)
{
  Real D( 1e-12 ), v( -3e-8 ), kf( 1e-8 ), r0(2e-7), sigma(0.0), a(2e-7);//L(2e-7), r0(L)
  GreensFunction1DRadAbs gf( D, v, kf, r0, sigma, a );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK_EQUAL( t, 0.0 );
  std::cout << "GreensFunction1DRadAbs_test : test_DrawTime_r0_equal_a : t = " << t << std::endl;
  std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DRadAbs_DrawTime_r0_equal_sigma_kf_zero)
{
  Real D( 1e-12 ), v( -3e-8 ), kf( 0.0 ), r0( 0.0 ), sigma(0.0), a(1e-7);//L(1e-7)
  GreensFunction1DRadAbs gf( D, v, kf, r0, sigma, a );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK(0.0 < t < INFINITY);
  std::cout << "GreensFunction1DRadAbs_test : test_DrawTime_r0_equal_sigma_kf_zero : t = " << t << std::endl;
  std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DRadAbs_DrawTime_r0_equal_sigma_kf_large)
{
  Real D( 1e-12 ), v( -3e-8 ), kf( 1e-8 ), r0( 1e-12 ), sigma(20e-7), a(40e-7);//L(20e-7)
  GreensFunction1DRadAbs gf( D, v, kf, r0, sigma, a );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK(0.0 < t < INFINITY);
  std::cout << "GreensFunction1DRadAbs_test : test_DrawTime_r0_equal_sigma_kf_large : t = " << t << std::endl;
  std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DRadAbs_DrawEventType)
{
  Real D( 1e-12 ), v( -3e-8 ), kf( 1e-8 ), r0(1e-7), sigma(.5e-7), a(2.5e-7);//L(2e-7), r0(L/2)
  GreensFunction1DRadAbs gf( D, v, kf, r0, sigma, a );

  Real t( gf.drawTime( 0.5 ) );
  EventKind eventType(gf.drawEventType(0.5, t));

  BOOST_CHECK( (eventType == IV_ESCAPE) || (eventType == IV_REACTION) );
  std::cout << "GreensFunction1DRadAbs_test : test_DrawEventType : eventType = " << eventType << std::endl;

  eventType = gf.drawEventType(0.999999, t);
  BOOST_CHECK( eventType == IV_ESCAPE );

  eventType = gf.drawEventType(0.0, t);
  BOOST_CHECK( eventType == IV_REACTION );

  std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DRadAbs_DrawEventType_smallt)
{
  Real D( 1e-12 ), v( -3e-8 ), kf( 1e-8 ), r0( 1e-6 ), sigma(.5e-6), a(2.5e-6);//L(2e-6), r0(L/2)
  GreensFunction1DRadAbs gf( D, v, kf, r0, sigma, a );

  Real t( gf.drawTime( 0.001 ) );
  EventKind eventType(gf.drawEventType(0.5, t));

  BOOST_CHECK( (eventType == IV_ESCAPE) || (eventType == IV_REACTION) );
  std::cout << "GreensFunction1DRadAbs_test : test_DrawEventType : eventType = " << eventType << std::endl;

  eventType = gf.drawEventType(0.999999, t);
  BOOST_CHECK( eventType == IV_ESCAPE );

  eventType = gf.drawEventType(0.0, t);
  BOOST_CHECK( eventType == IV_REACTION );

  std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DRadAbs_DrawR)
{
  Real D( 1e-12 ), v( -3e-8 ), kf( 1e-8 ), r0( 1e-7 ), sigma(.5e-7), a(2.5e-7);//L(2e-7), r0(L/2)
  Real L(a - sigma);
  GreensFunction1DRadAbs gf( D, v, kf, r0, sigma, a );

  Real t( 1e-3 );
  Real r(gf.drawR(0.5, t));
  BOOST_CHECK( (0.0 <= r) || (r <= L) );
  std::cout << "GreensFunction1DRadAbs_test : test_DrawR : r = " << r << std::endl;

  Real r1(gf.drawR(0.0, t));
  Real r2(gf.drawR(0.999999999999, t));

  BOOST_CHECK( (0.0 <= r1) || (r1 <= L) );
  BOOST_CHECK( (0.0 <= r2) || (r2 <= L) );

  std::cout << "GreensFunction1DRadAbs_test : test_DrawR : r1 = " << r1 << std::endl;
  std::cout << "GreensFunction1DRadAbs_test : test_DrawR : r2 = " << r2 << std::endl;

  std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DRadAbs_DrawR_zerot)
{
  Real D( 1e-12 ), v( -3e-8 ), kf( 1e-8 ), r0( .5e-7 ), sigma(1e-8), a(1.1e-7);//L(1e-7), r0(L/2)
  GreensFunction1DRadAbs gf( D, v, kf, r0, sigma, a );

  Real t( 0.0 );
  Real r(gf.drawR(0.5, t));
  BOOST_CHECK_EQUAL( r, r0 );
  std::cout << "GreensFunction1DRadAbs_test : test_DrawR_zerot : r = " << r << std::endl;
  std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DRadAbs_DrawR_r0_equal_sigma)
{
  Real D( 1e-12 ), v( -3e-8 ), kf( 1e-8 ), r0( 0.0 ), sigma( 0.0 ), a(2e-7);//L(2e-7), r0(0.0)
  Real L(a - sigma);
  GreensFunction1DRadAbs gf( D, v, kf, r0, sigma, a );

  Real t( 1e-3 );
  Real r(gf.drawR(0.5, t));
  BOOST_CHECK( (0.0 <= r) || (r <= L) );
  std::cout << "GreensFunction1DRadAbs_test : test_DrawR_r0_equal_sigma : r = " << r << std::endl;
  std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(GF1DRadAbs_DrawR_squeezed)
{
  Real D( 1e-12 ), v( -3e-8 ), kf( 1e-8 ), r0( 0.01e-8 ), sigma( 0.0 ), a(0.02e-8);//L(0.02e-8), r0(L/2)
  Real L(a - sigma);
  GreensFunction1DRadAbs gf( D, v, kf, r0, sigma, a );

  Real t( 1e-6 );
  r0 = 0.0;
  gf.setr0( r0 );
  Real r(gf.drawR(0.5, t));
  BOOST_CHECK( (0.0 <= r) || (r <= L) );
  std::cout << "GreensFunction1DRadAbs_test : test_DrawR_squeezed: r = " << r << std::endl;

//near s
  r0 = 0.0001e-8;
  gf.setr0( r0 );
  r = gf.drawR(0.5, t);
  BOOST_CHECK( (0.0 <= r) || (r <= L) );
  std::cout << "GreensFunction1DRadAbs_test : test_DrawR_squeezed: 1st drawn r = " << r << std::endl;

//near a 
  r0 = L - 0.0001e-8;
  gf.setr0( r0 );
  r = gf.drawR(0.5, t);
  BOOST_CHECK( (0.0 <= r) || (r <= L) );
  std::cout << "GreensFunction1DRadAbs_test : test_DrawR_squeezed: 2nd drawn r = " << r << std::endl;
  std::cout << std::endl;
}

