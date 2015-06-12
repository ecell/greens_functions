#define BOOST_TEST_MODULE "GreensFunction3D_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include "../GreensFunction.hpp"
#include "../GreensFunction3D.hpp"

using namespace greens_functions;

BOOST_AUTO_TEST_CASE(GF3D_constructor)
{
  Real D( 1e-12 ), r0( 1e-7 );
  GreensFunction3D gf( D, r0 );
}

BOOST_AUTO_TEST_CASE(GF3D_DrawR)
{
  Real D( 1e-12 ), r0( 2e-8 );
  GreensFunction3D gf( D, r0 );

  Real t( 1e-3 );

  Real r( gf.drawR(0.5, t) );
  BOOST_CHECK( 0.0 <= r );

  r = gf.drawR( 0.0 , t);
  BOOST_CHECK(0.0 <= r);

  r = gf.drawR( 1.0, t );
  BOOST_CHECK(0.0 <= r);
}

BOOST_AUTO_TEST_CASE(GF3D_DrawR_zerot_is_r0)
{
  Real D( 1e-12 ), r0( 2e-8 );
  GreensFunction3D gf( D, r0 );

  Real t( 0.0 );

  Real r( gf.drawR(0.5, t) );
  BOOST_CHECK_EQUAL( r, r0 );

  r = gf.drawR( 0.0 , t);
  BOOST_CHECK_EQUAL( r, r0 );

  r = gf.drawR( 1.0, t );
  BOOST_CHECK_EQUAL( r, r0 );
}

BOOST_AUTO_TEST_CASE(GF3D_DrawR_smallt)
{
  Real D( 1e-12 ), r0( 2e-8 );
  GreensFunction3D gf( D, r0 );

  Real t( 1e-4 ), r;

  while(t > 1e-60)
  {
    r = gf.drawR(0.5, t);
    BOOST_CHECK( 0.0 <= r );

    r = gf.drawR( 0.0 , t);
    BOOST_CHECK( 0.0 <= r );

    r = gf.drawR( 1.0, t );
    BOOST_CHECK( 0.0 <= r );

    t *= 1e-3;
  }
}

BOOST_AUTO_TEST_CASE(GF3D_DrawTheta)
{
  Real D( 1e-12 ), r0( 5e-8 );
  GreensFunction3D gf( D, r0 );

  Real t( 1e-4 ), r( r0 ), theta;

  theta = gf.drawTheta( 0.5, r, t );
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );

  theta = gf.drawTheta( 0.0, r, t);
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );

  theta = gf.drawTheta( 1.0, r, t );
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );
}

BOOST_AUTO_TEST_CASE(GF3D_ip_r_infinity_is_one)
{
  Real D( 1e-12 ), r0( 5e-8 );
  GreensFunction3D gf( D, r0 );

  Real t( 1e-5 ), r( 2.5e-8 ), ip;

  ip = gf.ip_r( INFINITY, t );
  BOOST_CHECK_EQUAL( 1.0 ,ip );
}

BOOST_AUTO_TEST_CASE(GF3D_p_r_is_ip_r)
{
/*
import scipy.integrate

        D = 1e-12
        
        t = 1e-5
        r0 = 5e-8
        r = 2.5e-8
        
        gf = mod.GreensFunction3D(D, r0)

        ip = gf.ip_r(0.0, t)
        self.assertEqual(0.0, ip)

        maxr = 1e-6

        resolution = 20
        for i in range(1, resolution):
            r = i * maxr / resolution 
            ip = gf.ip_r(r, t) 
            result = scipy.integrate.quad(gf.p_r, 0.0, r,
                                          args=(t, ))
            np = result[0]
            self.assertAlmostEqual(0.0, (np-ip)/ip)


*/
}

BOOST_AUTO_TEST_CASE(GF3D_int_p_theta_is_p_r)
{
/*
 import scipy.integrate

        D = 1e-12
        
        t = 1e-5
        r0 = 5e-8
        r = 2.5e-8
        
        gf = mod.GreensFunction3D(D, r0)

        ip = gf.ip_theta(numpy.pi, r, t)
        result = scipy.integrate.quad(gf.p_theta, 0.0, numpy.pi,
                                      args=(r, t))
        np = result[0]

        pr = gf.p_r(r, t) / (2 * numpy.pi * r * r)
        
        self.assertAlmostEqual(0.0, (pr-ip)/pr)
        self.assertAlmostEqual(0.0, (pr-np)/pr)

*/
}

BOOST_AUTO_TEST_CASE(GF3D_int_p_theta_is_ip_theta)
{

/*
        import scipy.integrate

        D = 1e-12
        
        t = 1e-3
        r0 = 5e-8
        r = 2.5e-8
        
        gf = mod.GreensFunction3D(D, r0)

        ip = gf.ip_theta(0.0, r, t)
        self.assertEqual(0.0, ip)

        resolution = 20
        for i in range(1, resolution):
            theta = i * numpy.pi / resolution 
            ip = gf.ip_theta(theta, r, t)
            result = scipy.integrate.quad(gf.p_theta, 0.0, theta,
                                          args=(r, t))
            np = result[0]
            self.assertAlmostEqual(0.0, (np-ip)/ip)
*/
}
