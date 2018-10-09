#define BOOST_TEST_MODULE "GreensFunction3DRadInf_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include "../freeFunctions.hpp"
#include "../GreensFunction.hpp"
#include "../GreensFunction3DRadInf.hpp"

using namespace greens_functions;
using namespace std;
using namespace boost::test_tools;

BOOST_AUTO_TEST_CASE(GF3DRadInf_constructor)
{
  Real D( 1e-12 ), kf( 1e-8 ), r0( 5e-8 ), sigma( 1e-8 );
  GreensFunction3DRadInf gf( D, kf, r0, sigma);
}

BOOST_AUTO_TEST_CASE(GF3DRadInf_DrawTime)
{
  Real D( 1e-12 ), kf( 1e-8 ), r0( 5e-8 ), sigma( 1e-8 );
  GreensFunction3DRadInf gf( D, kf, r0, sigma );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK( 0.0 < t );

  t = gf.drawTime( 0.0 );
  BOOST_CHECK( 0.0 <= t );

  t = gf.drawTime( 0.9999999 ) ;
  BOOST_CHECK( 0.0 < t );
}

BOOST_AUTO_TEST_CASE(GF3DRadInf_DrawTime_r0_equal_sigma_kf_zero)
{
  Real D( 1e-12 ), kf( 0.0 ), sigma( 1e-8 ), r0( sigma );
  GreensFunction3DRadInf gf( D, kf, r0, sigma );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK( 0.0 <= t );
}

BOOST_AUTO_TEST_CASE(GF3DRadInf_DrawR)
{
  Real D( 1e-12 ), kf( 1e-8 ), r0( 2e-8 ), sigma( 1e-8 );
  GreensFunction3DRadInf gf( D, kf, r0, sigma );

  Real t( 1e-3 );

  Real r( gf.drawR( 0.5, t ) );
  BOOST_CHECK( sigma <= r );

  Real r1( gf.drawR( 0.0, t ) );
  Real r2( gf.drawR( 0.9999999, t ) );
  BOOST_CHECK( sigma <= r1 );
  BOOST_CHECK( sigma <= r2 );

  BOOST_CHECK( fabs(r1 - sigma) <= 1e-15 );
}

BOOST_AUTO_TEST_CASE(GF3DRadInf_DrawR_zerot)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), r0( 2e-8 );
  GreensFunction3DRadInf gf( D, kf, r0, sigma );

  Real t( 0.0 );

  Real r( gf.drawR( 0.5, t ) );
  BOOST_CHECK_EQUAL( r0, r );
}

BOOST_AUTO_TEST_CASE(GF3DRadInf_DrawR_smallt)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), r0( 1.000001e-8 );
  GreensFunction3DRadInf gf( D, kf, r0, sigma );

  Real t( 1e-12 );

  Real r( gf.drawR( 0.5, t ) );
  BOOST_CHECK( sigma <= r );
}

BOOST_AUTO_TEST_CASE(GF3DRadInf_DrawR_r0_equal_sigma)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), r0( sigma );
  GreensFunction3DRadInf gf( D, kf, r0, sigma );

  Real t( 1e-5 );

  Real r( gf.drawR( 0.5, t ) );
  BOOST_CHECK( sigma <= r );
}

BOOST_AUTO_TEST_CASE(GF3DRadInf_DrawTheta)
{
  Real D( 1e-12 ), kf( 1e-8 ), r0( 2e-8 ), sigma( 1e-8 );
  GreensFunction3DRadInf gf( D, kf, r0, sigma );

  Real t( 1e-3 ), r( 2.1e-8 );

  Real theta( gf.drawTheta( 0.5, r, t ) );
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );

  theta = gf.drawTheta( 0.0, r, t );
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );

  theta = gf.drawTheta( 0.9999999, r, t );
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );
}

/*
    '''
    def test_draw_theta2(self):
        D = 2e-12
        kf = 0
        sigma = 5e-9
        r0 = 5.064e-9
        r = 5.05e-9
        t = 1e-9
        gf = mod.GreensFunction3DRadInf(D, kf, sigma)
        theta = gf.drawTheta(0.5, r, r0, t)
        self.failIf(theta < 0.0 or theta > numpy.pi)
        #theta = gf.drawTheta(0.0, r, r0, t)
        #self.failIf(theta < 0.0 or theta > numpy.pi)
        #theta = gf.drawTheta(0.9999999, r, r0, t)
        #self.failIf(theta < 0.0 or theta > numpy.pi)
'''
*/

BOOST_AUTO_TEST_CASE(GF3DRadInf_DrawTheta_zerot)
{
  Real D( 1e-12 ), kf( 1e-8 ), r0( 5e-8 ), sigma( 1e-8 );
  GreensFunction3DRadInf gf( D, kf, r0, sigma );

  Real t( 0.0 );
  // Real t( 0.0 ), r( 5e-8 );

  Real theta( gf.drawTheta( 0.5, r0, t ) );
  BOOST_CHECK_EQUAL( 0.0, theta );
}

BOOST_AUTO_TEST_CASE(GF3DRadInf_DrawTheta_smallt)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 );
  Real t( 1e-11 ), r, r0;

  Real disp( 3 * sqrt(6 * D * t) );
  r = sigma + disp + disp;
  r0 = sigma + disp;

  GreensFunction3DRadInf gf( D, kf, r0, sigma );

  Real theta( gf.drawTheta( 0.5, r, t ) );
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );
}

BOOST_AUTO_TEST_CASE(GF3DRadInf_DrawTheta_r0_equal_sigma)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), r0( sigma );
  GreensFunction3DRadInf gf( D, kf, r0, sigma );

  Real t( 1e-3 ), r( r0 );

  Real theta( gf.drawTheta( 0.5, r, t ) );
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );
}

BOOST_AUTO_TEST_CASE(GF3DRadInf_p_int_r_at_s_is_zero)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), r0( 2e-8 );
  GreensFunction3DRadInf gf( D, kf, r0, sigma );

  Real t( 1e-3 );

  Real pintr( gf.p_int_r( sigma, t ) );
  BOOST_CHECK_EQUAL( 0.0, pintr );
}

BOOST_AUTO_TEST_CASE(GF3DRadInf_p_int_r0_at_s_zerot_is_zero)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), r0( 2e-8 );
  GreensFunction3DRadInf gf( D, kf, r0, sigma );

  Real t( 0.0 );

  Real pintr( gf.p_int_r( sigma, t ) );
  BOOST_CHECK_EQUAL( 0.0, pintr );
}

BOOST_AUTO_TEST_CASE(GF3DRadInf_p_int_r_large_is_p_survival)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), r0( 2e-8 );
  GreensFunction3DRadInf gf( D, kf, r0, sigma );

  Real t( 1e-3 );

  Real pintr( gf.p_int_r( sigma * 1e8, t ) );
  Real psurv( gf.p_survival(t) );
  BOOST_CHECK_CLOSE( psurv, pintr, 0.0001 );//?
}

/*
BOOST_AUTO_TEST_CASE(GF3DRadInf_ip_theta_is_int_p_theta)
{
        import scipy.integrate

        D = 1e-12
        sigma = 1e-8
        kf = 1e-9

        t = 1e-4
        r0 = sigma*1.1

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)
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
            #print 'theta, np, ip', theta, np, ip
            self.assertAlmostEqual(0.0, (np-ip)/ip)
}
*/

/*
BOOST_AUTO_TEST_CASE(GF3DRadInf_ip_theta_r0_is_sigma)
{
        import scipy.integrate

        D = 1e-12
        sigma = 1e-8
        kf = 1e-8

        t = 1e-3
        r0 = sigma

        gf = mod.GreensFunction3DRadInf(D, kf, r0, sigma)
        r = 1.1 * r0

        ip = gf.ip_theta(0.0, r, t)
        self.assertEqual(0.0, ip)

        ip = gf.ip_theta(numpy.pi, r, t)
        #print 'ip', ip
        #self.assertEqual(0.0, ip)
}
*/

BOOST_AUTO_TEST_CASE(GF3DRadInf_ip_theta_pi_is_p_irr)
{
  Real D( 1e-12 ), kf( 0.0 ), sigma( 1e-8 ), r0( 1.1e-8 );
  GreensFunction3DRadInf gf( D, kf, r0, sigma );

  Real t( 1e-3 ), r( r0 );

  Real ip( gf.ip_theta( M_PI, r, t ) * ( 2 * M_PI * r * r ) );
  Real pirr( p_irr(r, t, r0, kf, D, sigma) );
  gf.p_corr(M_PI, r, t);
  // Real pcorr( gf.p_corr(M_PI, r, t) * ( 2 * M_PI * r * r ) );
  Real pfree( gf.p_free(M_PI, r, t) * ( 2 * M_PI * r * r ) );

  Real dif( fabs(pirr - pfree) );
  BOOST_CHECK( dif > 1e-7 );//NotAlmostEqual
  BOOST_CHECK( 0.0 != ip );
  BOOST_CHECK_CLOSE( ip/pirr, 1, 0.0001 );
}

BOOST_AUTO_TEST_CASE(GF3DRadInf_p_theta_never_negative)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), r0( 5e-8 );
  GreensFunction3DRadInf gf( D, kf, r0, sigma );

  Real t( 1e-3 ), r( r0 ), theta;

  Real pint( gf.ip_theta( M_PI, r, t ) );
  Real pmin(0.0), p;
  int resolution(50);

  for(int i=0; i<resolution; i++)
  {
    theta = i * M_PI / resolution;
    p = gf.p_theta(theta, r, t) / pint / resolution;
    pmin = min(pmin, p);
  }

  BOOST_CHECK( 0.0 <= pmin );
/*python code: 
        resolution = 50
        for i in range(resolution):
            theta = i * numpy.pi / resolution
            p = gf.p_theta(theta, r, t) / pint / resolution 
            pmin = min(pmin, p)
            #print 'theta: ', theta, '\tp: ', p
            
        self.failIf(pmin < 0.0, 'Negative p_theta; t= %.16g, %s'
                    % (t, gf.dump()))
*/
}

BOOST_AUTO_TEST_CASE(GF3DRadInf_ip_theta_never_decrese)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), r0( 5e-8 );
  GreensFunction3DRadInf gf( D, kf, r0, sigma );

  Real t( 1e-3 ), r( r0 ), theta;

  Real pint, pint_prev( 0.0 );
  int resolution(50);

  for(int i=0; i<resolution; i++)
  {
    theta = i * M_PI / resolution;
    pint = gf.ip_theta(theta, r, t);
    BOOST_CHECK( pint_prev <= pint );
    pint_prev = pint;
  }
/*python code:
        resolution = 50
        for i in range(resolution):
            theta = i * numpy.pi / resolution
            pint = gf.ip_theta(theta, r, t)
            #print pint
            self.failIf(pint < pint_prev)
            pint_prev = pint
*/
}

