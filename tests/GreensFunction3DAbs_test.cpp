#define BOOST_TEST_MODULE "GreensFunction3DAbs_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include "../GreensFunction.hpp"
#include "../GreensFunction3DAbs.hpp"

using namespace greens_functions;

BOOST_AUTO_TEST_CASE(GF3DAbs_constructor)
{
  Real D( 1e-12 ), r0( 5e-8 ), sigma( 1e-8 ), a( 1e-7 );
  GreensFunction3DAbs gf( D, r0, a );
}

BOOST_AUTO_TEST_CASE(GF3DABs_DrawTime)
{
  Real D( 1e-12 ), a( 1e-7 ), r0( 5e-8 );
  GreensFunction3DAbs gf( D, r0, a );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK(0.0 < t && t < INFINITY);

  t = gf.drawTime( 0.0 );
  BOOST_CHECK(0.0 < t && t < INFINITY);

  t = gf.drawTime( 1.0 ) ;
  BOOST_CHECK(0.0 < t && t < INFINITY);
}

BOOST_AUTO_TEST_CASE(GF3DABs_DrawTime_r0_equal_a)
{
  Real D( 1e-12 ), a( 1e-7 ), r0( a );
  GreensFunction3DAbs gf( D, r0, a );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK_EQUAL(0.0 , t );
}

BOOST_AUTO_TEST_CASE(GF3DABs_DrawR)
{
  Real D( 1e-12 ), a( 1e-7 ), r0( 2e-8 );
  GreensFunction3DAbs gf( D, r0, a );

  Real t( 1e-3 );

  Real r( gf.drawR(0.5, t) );
  BOOST_CHECK( 0.0 <= r && r <= a);

  Real r1( gf.drawR( 0.0, t ) );
  Real r2( gf.drawR( 1.0, t ) );
  BOOST_CHECK( 0.0 <= r1 && r1 <= a);
  BOOST_CHECK( 0.0 <= r2 && r2 <= a);

  BOOST_CHECK( abs(r1) <= 1e-15 );
  BOOST_CHECK_CLOSE( abs(r2 - a), 0.0, 0.0001);
}

BOOST_AUTO_TEST_CASE(GF3DABs_DrawR_zerot)
{
  Real D( 1e-12 ), a( 1e-7 ), r0( 2e-8 );
  GreensFunction3DAbs gf( D, r0, a );

  Real t( 0.0 );

  Real r( gf.drawR(0.5, t) );
  BOOST_CHECK_EQUAL( r0, r );
}

BOOST_AUTO_TEST_CASE(GF3DABs_DrawR_squeezed)
{
  Real D( 1e-12 ), a( 1.01e-8 ), r0( 1.005e-8 );
  GreensFunction3DAbs gf( D, r0, a );

  Real t( 1e-6 );

  Real r( gf.drawR(0.5, t) );
  BOOST_CHECK( 0.0 <= r && r <= a );

//near 0
  r0 = 0.0001e-8;
  GreensFunction3DAbs gf1( D, r0, a );
  r = gf1.drawR(0.5, t);
  BOOST_CHECK( 0.0 <= r && r <= a );

//near a
  r0 = 1.0099e-8;
  GreensFunction3DAbs gf2( D, r0, a );
  r = gf2.drawR(0.5, t);
  BOOST_CHECK( 0.0 <= r && r <= a );
}

BOOST_AUTO_TEST_CASE(GF3DABs_DrawTheta)
{
  Real D( 1e-12 ), a( 1e-7 ), r0( 5e-8 );
  GreensFunction3DAbs gf( D, r0, a );

  Real t( gf.drawTime(0.5) );
  Real r( gf.drawR(0.5, t) );
  Real theta;

  theta = gf.drawTheta(0.5, r, t);
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );

  theta = gf.drawTheta(0.0, r, t);
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );

  theta = gf.drawTheta(1.0, r, t);
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );
}

BOOST_AUTO_TEST_CASE(GF3DABs_DrawTheta_zerot)
{
  Real D( 1e-12 ), a( 1e-7 ), r0( 5e-8 );
  GreensFunction3DAbs gf( D, r0, a );

  Real t( 0.0 );
  Real theta(gf.drawTheta(0.5, r0, t));
  BOOST_CHECK_EQUAL( 0.0, theta );
}

BOOST_AUTO_TEST_CASE(GF3DABs_DrawTheta_smallt)
{
  Real D( 1e-12 ), a( 1e-7 ), r0( 5e-8 );
  GreensFunction3DAbs gf( D, r0, a );

  Real r( 5e-8 );
  Real t( 1e-4 );//it is not very small though
  Real theta(gf.drawTheta(0.5, r, t));
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );
}

BOOST_AUTO_TEST_CASE(GF3DABs_DrawTheta_larget)
{
  Real D( 1e-12 ), a( 1e-7 ), r0( 5e-8 );
  GreensFunction3DAbs gf( D, r0, a );

  Real r( 5e-8 );
  Real t( 1e5 );
  Real theta(gf.drawTheta(0.5, r, t));
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );
}

BOOST_AUTO_TEST_CASE(GF3DABs_DrawTheta_near_a)
{
  Real D( 1e-12 ), a( 1.01e-8 ), r0( 1.009999e-8 ), r( 1.009999e-08 );
  GreensFunction3DAbs gf( D, r0, a );

  Real t( 1e-5 );
  Real theta(gf.drawTheta(0.5, r, t));
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );
}

BOOST_AUTO_TEST_CASE(GF3DABs_DrawTheta_r_equal_a)
{
  Real D( 1e-12 ), a( 1e-7 ), r0( 9e-8 ), r( a );
  GreensFunction3DAbs gf( D, r0, a );

  Real t( 1e-4 );
  Real theta(gf.drawTheta(0.5, r, t));
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );
}

BOOST_AUTO_TEST_CASE(GF3DABs_DrawTheta_r0_near_a_r_equal_a)
{
  Real D( 1e-12 ), a( 1e-7 ), r0( a - 1e-9 ), r( a );
  GreensFunction3DAbs gf( D, r0, a );

  Real t( 1e-2 );
  Real theta(gf.drawTheta(0.5, r, t));
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );
}

BOOST_AUTO_TEST_CASE(GF3DAbs_p_int_r)
{
  Real D( 1e-12 ), a( 1e-7 ), r0( 5e-8 ), r( r0 );
  GreensFunction3DAbs gf( D, r0, a );

  Real t( 1e-3 );
  Real pintr( gf.p_int_r( r, t ) );
  BOOST_CHECK( 0.0 <= pintr && pintr <= 1.0 );
}

BOOST_AUTO_TEST_CASE(GF3DAbs_p_int_r_at_a_is_p_survival)
{
  Real D( 1e-12 ), a( 1e-7 ), r0( 5e-8 ), r( a );//python code says r = r0 and p_int_r( a,t )
  GreensFunction3DAbs gf( D, r0, a );

  Real t( 1e-3 );
  Real pintr( gf.p_int_r( r, t ) );
  Real psurv( gf.p_survival( t ) );

  BOOST_CHECK_CLOSE( pintr, psurv, 0.0001 );
}

BOOST_AUTO_TEST_CASE(GF3DAbs_p_int_r_at_zero_is_zero)
{
  Real D( 1e-12 ), a( 1e-7 ), r0( 5e-8 ), r( r0 );
  GreensFunction3DAbs gf( D, r0, a );

  Real t( 1e-3 );
  Real pintr( gf.p_int_r( 0.0, t ) );

  BOOST_CHECK_EQUAL( pintr, 0.0 );
}

/*
    def test_ip_theta_is_int_p_theta(self):

        import scipy.integrate

        D = 1e-11

        t = 1e-3
        r0 = 5e-8

        a = 1e-7
        
        gf = mod.GreensFunction3DAbs(D, r0, a)
        r = r0

        ip = gf.ip_theta(0.0, r, t)
        self.assertEqual(0.0, ip)
        
        resolution = 10
        for i in range(1, resolution):
            theta = i * numpy.pi / resolution 
            ip = gf.ip_theta(theta, r, t)
            result = scipy.integrate.quad(gf.p_theta, 0.0, theta,
                                          args=(r, t))
            np = result[0]
            #print theta, np, ip
            self.assertAlmostEqual(0.0, (np-ip)/ip)

    '''
    def test_ip_theta_pi_is_p_0(self):
        D = 1e-12
        sigma = 1e-8
        kf = 1e-8
        t = 1e-5
        r0 = 5e-8
        r = r0
        a = 1e-7
        
        gf = mod.GreensFunction3DAbs(D, kf, sigma, a)
        ip = gf.ip_theta(numpy.pi, r, r0, t)
        p0 = gf.p_0(t, r, r0) * 2
        self.assertNotEqual(0.0, ip)
        self.assertAlmostEqual(1.0, ip/p0)
'''
*/

BOOST_AUTO_TEST_CASE(GF3DAbs_p_theta_never_negative)
{
  Real D( 1e-12 ), a( 1e-7 ), r0( 5e-8 ), r( r0 );
  GreensFunction3DAbs gf( D, r0, a );

  Real t( 1e-3 );
  Real pint( gf.ip_theta( M_PI, r, t ) );
  Real pmin(0.0), p, theta;
  int resolution(50);

  for(int i(0); i<resolution; i++)
  {
    theta = i * M_PI / resolution;
    p = gf.p_theta(theta, r, t) / pint / resolution;
    pmin = std::min( pmin, p );
  }

  BOOST_CHECK( 0.0 <= pmin );
}

BOOST_AUTO_TEST_CASE(GF3DAbs_p_theta_never_decrease)
{
  Real D( 1e-12 ), a( 1e-7 ), r0( 5e-8 ), r( r0 );
  GreensFunction3DAbs gf( D, r0, a );

  Real t( 1e-3 );
  Real pint_prev( 0.0 ), pint, theta;
  int resolution(50);

  for(int i(0); i<resolution; i++)
  {
    theta = i * M_PI / resolution;
    pint = gf.ip_theta(theta, r, t);
    BOOST_CHECK( pint >= pint_prev );
    pint_prev = pint;
  }
}

BOOST_AUTO_TEST_CASE(GF3DAbs_idp_theta_at_a_is_dp_survival)
{
  Real D( 1e-12 ), a( 1e-7 ), r0( 9e-8 );
  GreensFunction3DAbs gf( D, r0, a );

  Real t( 1e-3 );
  Real dp( gf.dp_survival(t) );
  Real iptheta( gf.idp_theta( M_PI, a, t ) * M_PI * a * a * 2);

  BOOST_CHECK_CLOSE(dp, iptheta, 0.0001);
}

