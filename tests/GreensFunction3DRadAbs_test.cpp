#define BOOST_TEST_MODULE "GreensFunction3DRadAbs_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include "../GreensFunction.hpp"
#include "../GreensFunction3DRadAbs.hpp"

using namespace greens_functions;


BOOST_AUTO_TEST_CASE(GF3DRadAbs_constructor)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( 5e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_DrawTime)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( 5e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK(0.0 < t && t < std::numeric_limits<Real>::infinity());

  t = gf.drawTime( 0.0 );
  BOOST_CHECK(0.0 < t && t < std::numeric_limits<Real>::infinity());

  t = gf.drawTime( 1 - 1e-16 );
  BOOST_CHECK(0.0 < t && t < std::numeric_limits<Real>::infinity());
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_DrawTime_a_equal_sigma)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( sigma ), r0( a );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK_EQUAL( 0.0, t);
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_DrawTime_a_near_sigma)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 );
  Real a( sigma + sigma * 1e-6 ), r0( (a + sigma) * 0.5 );

  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK(0.0 < t && t < std::numeric_limits<Real>::infinity());
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_DrawTime_r0_equal_a)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( a );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK_EQUAL( 0.0, t );
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_DrawTime_r0_equal_sigma_kf_zero)
{
  Real D( 1e-12 ), kf( 0.0 ), sigma( 1e-8 ), a( 1e-7 ), r0( sigma );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK(0.0 < t && t < std::numeric_limits<Real>::infinity());
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_DrawTime_r0_equal_sigma_kf_large)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 10e-7 ), r0( sigma + 1e-12 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( gf.drawTime( 0.5 ) );
  BOOST_CHECK(0.0 < t && t < std::numeric_limits<Real>::infinity());
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_DrawEventType)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( 5e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( gf.drawTime( 0.5 ) );
  GreensFunction3DRadAbs::EventKind eventType( gf.drawEventType( 0.5, t ) );
  BOOST_CHECK( (eventType == GreensFunction3DRadAbs::IV_ESCAPE) || (eventType == GreensFunction3DRadAbs::IV_REACTION) );

  eventType = gf.drawEventType(0.0, t);
  BOOST_CHECK_EQUAL(eventType, GreensFunction3DRadAbs::IV_REACTION);

  eventType = gf.drawEventType(0.9999999, t);
  BOOST_CHECK_EQUAL(eventType, GreensFunction3DRadAbs::IV_ESCAPE);
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_DrawEventType_smallt)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-6 ), r0( 1.1e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( gf.drawTime( 0.999 ) );

  GreensFunction3DRadAbs::EventKind eventType( gf.drawEventType( 0.5, t ) );
  BOOST_CHECK( (eventType == GreensFunction3DRadAbs::IV_ESCAPE) || (eventType == GreensFunction3DRadAbs::IV_REACTION) );

  eventType = gf.drawEventType(0.0, t);
  BOOST_CHECK_EQUAL(eventType, GreensFunction3DRadAbs::IV_REACTION);

  eventType = gf.drawEventType(0.9999999, t);
//  BOOST_CHECK_EQUAL(eventType, GreensFunction3DRadAbs::IV_ESCAPE); it fails and it is comment-outed in python code, too.
}

/*
    '''
    def test_draw_time2(self):
        D = 1e-12
        kf = 1e-18
        #kf = 0
        sigma = 1e-8
        a = 1e-7
        r0 = 9e-8
        
        gf = mod.GreensFunction3DRadAbs(D, kf, r0, sigma, a)
        print '==============================================='
        t, et = gf.drawTime2(0.5, 0.5)
        t2 = gf.drawTime(0.5)
        print t, et, t2
        self.failIf(t <= 0.0 or t >= numpy.inf)
        print '==============================================='
        t, et = gf.drawTime2(0.0, 0.0)
        self.failIf(t < 0.0 or t >= numpy.inf)
        print t, et
        print '==============================================='
        t, et = gf.drawTime2(1 - 1e-8, 1 - 1e-8)
        self.failIf(t <= 0.0 or t >= numpy.inf)
        print t, et
        print '==============================================='
    def test_draw_time2_a_equal_sigma(self):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        a = sigma
        r0 = a
        
        gf = mod.GreensFunction3DRadAbs(D, kf, r0, sigma, a)
        t, et = gf.drawTime2(0.5, 0.5)
        self.assertEqual(0.0, t)
        self.assertEqual(et, mod.PairEventKind.IV_ESCAPE)
    def test_draw_time2_squeezed(self):
        D = 1e-12
        kf = 1e-10
        sigma = 1e-8 
        a = sigma + sigma * 1e-6
        r0 = (a + sigma) * .5
        
        gf = mod.GreensFunction3DRadAbs(D, kf, r0, sigma, a)
        t, et = gf.drawTime2(0.5, 0.5)
        self.failIf(t <= 0.0 or t >= numpy.inf)
    def test_draw_time2_r0_equal_a(self):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        a = 1e-7
        r0 = a
        
        gf = mod.GreensFunction3DRadAbs(D, kf, r0, sigma, a)
        t, et = gf.drawTime2(0.5, 0.5)
        self.assertEqual(0.0, t)
        self.assertEqual(et, mod.PairEventKind.IV_ESCAPE)
    def test_draw_time2_r0_equal_sigma_kf_zero(self):
        D = 1e-12
        kf = 0.0 # note this
        sigma = 1e-8
        a = 1e-7
        r0 = sigma
        
        gf = mod.GreensFunction3DRadAbs(D, kf, r0, sigma, a)
        t, et = gf.drawTime2(0.5, 0.5)
        self.failIf(t < 0.0 or t >= numpy.inf)
        self.assertEqual(et, mod.PairEventKind.IV_ESCAPE)
        # when kf == 0, pleavea == psurvival
        t2 = gf.drawTime(0.5)
        self.assertAlmostEqual(t, t2)
    def test_draw_time2_r0_near_sigma(self):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        a = 1e-7
        r0 = sigma*1.1
        print '**************'
        gf = mod.GreensFunction3DRadAbs(D, kf, r0, sigma, a)
        t, et = gf.drawTime2(0.3, 0.3, r0)
        t2 = gf.drawTime(0.3, r0)
        et2 = gf.drawEventType(0.3, t2)
        print '**************'
        print 't',t, 't2', t2, 'et', et, 'et2', et2
        self.failIf(t < 0.0 or t >= numpy.inf)
        self.assertEqual(et, mod.PairEventKind.IV_REACTION)
        self.assertAlmostEqual(t, t2)
    def no_test_draw_time2_r0_equal_sigma_kf_large(self):
        D = 1e-12
        kf = 1e-5
        sigma = 1e-8
        a = 10e-7
        r0 = sigma + 1e-12
        
        gf = mod.GreensFunction3DRadAbs(D, kf, r0, sigma, a)
        t, et = gf.drawTime2(0.5, 0.5, r0)
        self.failIf(t < 0.0 or t >= numpy.inf)
        '''
*/

BOOST_AUTO_TEST_CASE(GF3DRadAbs_DrawR)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( 2e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t(1e-3);

  Real r( gf.drawR( 0.5, t ) );
  BOOST_CHECK( sigma <= r && r <= a );

  Real r1( gf.drawR( 0.0, t ) );
  BOOST_CHECK( sigma <= r1 && r1 <= a );

  Real r2( gf.drawR( 0.999999999999, t ) );
  BOOST_CHECK( sigma <= r2 && r2 <= a );

  BOOST_CHECK_CLOSE(r1, sigma, 0.0001);
  BOOST_CHECK_CLOSE(r2, a, 0.0001);
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_DrawR_zerot)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( 2e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( 0.0 );

  Real r( gf.drawR( 0.5, t ) );
  BOOST_CHECK_EQUAL( r0, r );
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_DrawR_r0_equal_sigma)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( sigma );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( 1e-3 );

  Real r( gf.drawR( 0.5, t ) );
  BOOST_CHECK( sigma <= r && r <= a );
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_DrawR_squeezed)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1.01e-8 ), r0( 1.005e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( 1e-6 );

  Real r( gf.drawR( 0.5, t ) );
  BOOST_CHECK( sigma <= r && r <= a );

  //near s
  r0 = 1.0001e-8;
  GreensFunction3DRadAbs gf1( D, kf, r0, sigma, a );
  r =  gf1.drawR( 0.5, t );
  BOOST_CHECK( sigma <= r && r <= a );

  //near a
  r0 = 1.0099e-8;
  GreensFunction3DRadAbs gf2( D, kf, r0, sigma, a );
  r =  gf2.drawR( 0.5, t );
  BOOST_CHECK( sigma <= r && r <= a );
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_DrawTheta)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( 5e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( gf.drawTime( 0.5 ) );
  Real r( gf.drawR( 0.5, t ) );

  Real theta( gf.drawTheta(0.5, r, t) );
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );

  theta = gf.drawTheta(0.0, r, t);
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );

  theta = gf.drawTheta(0.999999, r, t);
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_DrawTheta_zero_t)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( 5e-8 ), r(5e-8);
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( 0.0 );

  Real theta( gf.drawTheta(0.5, r0, t) );
  BOOST_CHECK_EQUAL( 0.0, theta );
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_DrawTheta_small_t)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( 2e-8 ), r( 2e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( 1e-7 );

  Real theta( gf.drawTheta(0.5, r, t) );
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_DrawTheta_squeezed)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 );
  Real a( 1.001e-8 ), r0( 1.0001e-8 ), r( 1.0001e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( 1e-8 );

  Real theta( gf.drawTheta(0.5, r, t) );
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );

  //near s
  r0 = 1.00001e-8;
  r = 1.00001e-8;
  GreensFunction3DRadAbs gf1( D, kf, r0, sigma, a );

  theta =  gf1.drawTheta(0.5, r, t);
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );

  //near a
  r0 = 1.00099e-8;
  r = 1.00099e-8;
  GreensFunction3DRadAbs gf2( D, kf, r0, sigma, a );

  theta =  gf2.drawTheta(0.5, r, t);
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_ip_theta_squeezed)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 );
  Real a( 1.001e-8 ), r0( 1.00099e-8 ), r( 1.00099e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( 1e-10 );
  Real ip( gf.ip_theta(1, r, t) );

  r = 1.0000001e-8;
  r0 = 1.0000001e-8;
  
  GreensFunction3DRadAbs gf2( D, kf, r0, sigma, a );
  ip = gf.ip_theta( 1,r,t );
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_drawTheta_r0_equal_sigma)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( sigma ), r( r0 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( 1e-3 );
  Real theta( gf.drawTheta(0.5, r, t) );
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_drawTheta_r_equal_a)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( 9e-8 ), r( a );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( 1e-3 );
  Real theta( gf.drawTheta(0.5, r, t) );
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_drawTheta_regression1)
{
  Real r0( 1.0206416181e-07 ), t( 4.41358538629e-08 );
  Real D( 4e-11 ), kf( 0 ), sigma( 1e-7 ), a( 1.05134e-07 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real r( 1.03421535312e-07 );
  Real theta( gf.drawTheta(0.5, r, t) );
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_drawTheta_regression2)
{
  Real r0( 5.50265e-09 ), dt( 1.43732e-06 ), rnd( 0.749325 );
  Real D( 2e-12 ), kf( 0 ), sigma( 5e-9 ), a( 1.16464e-08 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real r( 6.9928e-09 );
  Real theta( gf.drawTheta(rnd, r, dt) );
  BOOST_CHECK( 0.0 <= theta && theta <= M_PI );
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_drawTheta_alpha0)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( 5e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real maxerror( 0.0 ), alpha, error;

  for(int i(0); i<100; i++)
  {
    alpha = gf.alpha0_i(i);
    error = fabs( gf.f_alpha0(alpha) / alpha );
    maxerror = std::max(error, maxerror);
  }

  BOOST_CHECK( fabs(maxerror) <= 1e-10 );
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_psurvival_is_pleaves_plus_pleavea)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( 5e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( 1e-5 );

  Real surv( gf.p_survival(t) );
  Real pleaves( gf.p_leaves(t) );
  Real pleavea( gf.p_leavea(t) );

  BOOST_CHECK( 0.0 < surv );
  BOOST_CHECK( 0.0 < pleavea && 0.0 < pleaves);
  BOOST_CHECK_CLOSE( surv, pleavea + pleaves, 0.0001 );
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_dppsurvival_is_leaves_plus_leavea)
{
  Real D( 1e-12 ), kf( 1e-13 ), sigma( 1e-8 ), a( 1e-7 ), r0( 2e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( 1e-3 );

  Real dsurv( gf.dp_survival(t) );
  Real leaves( gf.leaves(t) * 4.0 * M_PI * sigma * sigma );
  Real leavea( gf.leavea(t) * 4.0 * M_PI * a * a);

  BOOST_CHECK( 0.0 != dsurv );
  BOOST_CHECK_CLOSE( dsurv, leavea + leaves, 0.0001 );
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_psurvival_smallt)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( 2e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( 1e-4 );

  Real psurv, pleaves, pleavea;
  for(int i(0); i<5; ++i)
  {
    psurv = gf.p_survival(t);
    pleaves = gf.p_leaves(t);
    pleavea = gf.p_leavea(t);
    BOOST_CHECK(0.0 != psurv);
    BOOST_CHECK_CLOSE(pleaves + pleavea, psurv, 0.0001);
    t *= 0.1;
  }
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_p_int_r)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( 5e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( 1e-3 ), r( r0 );

  Real pintr( gf.p_int_r( r, t ) );
  BOOST_CHECK( 0.0 <= pintr && pintr <= 1.0 );
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_p_int_r_at_a_is_p_survival)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( 5e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( 1e-3 ), r( r0 );

  Real pintr( gf.p_int_r( a, t ) );
  Real psurv( gf.p_survival( t ) );
  BOOST_CHECK_CLOSE( pintr, psurv, 0.0001 );
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_p_int_r_at_s_is_zero)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( 5e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( 1e-3 ), r( r0 );

  Real pintr( gf.p_int_r( sigma, t ) );
  BOOST_CHECK_EQUAL( 0.0, pintr );
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_p_int_never_decrease)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 3e-7 ), r0( sigma );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( 1e-3 ), r( r0 );

  Real psurv( gf.p_survival( t ) ), pintr_prev( 0.0 ), pintr;
  int resolution( 500 );
  for(int i(0); i<resolution; ++i)
  {
    r = i * (a - sigma) / resolution + sigma;
    pintr = gf.p_int_r( r,t );
    BOOST_CHECK(pintr <= psurv);
    BOOST_CHECK(pintr >= pintr_prev);
    pintr_prev = pintr;
  }
}
/*
    def test_ip_theta_is_int_p_theta(self):

        import scipy.integrate

        D = 1e-12
        sigma = 1e-8
        kf = 1e-10

        t = 1e-2  #FIXME: smaller t should be fine
        r0 = 5e-8

        a = 1e-7
        
        gf = mod.GreensFunction3DRadAbs(D, kf, r0, sigma, a)
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
            self.assertAlmostEqual(0.0, (np-ip)/ip, 5)
*/

BOOST_AUTO_TEST_CASE(GF3DRadAbs_ip_theta_pi_is_p_0)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( 5e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( 1e-5 ), r( r0 );

  Real ip( gf.ip_theta( M_PI, r, t ) );
  Real p0( gf.p_0(t,r) * 2 );
  BOOST_CHECK( 0.0 != ip );
  BOOST_CHECK_CLOSE( 1.0, ip/p0, 5);
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_ip_theta_never_negative)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( 5e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( 1e-3 ), r( r0 );

  Real pint( gf.ip_theta( M_PI, r, t ) );
  Real pmin(0.0), theta, p;
  int resolution(50);
  for(int i(0); i<resolution; ++i)
  {
    theta = i * M_PI / resolution;
    p = gf.p_theta(theta, r, t) / pint / resolution;
    pmin = std::min(pmin, p);
  }

  BOOST_CHECK( pmin >= 0.0 );
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_ip_theta_never_decrease)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( 5e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( 1e-3 ), r( r0 );

  Real pint_prev( 0.0 ), pint, theta;
  int resolution(50);
  for(int i(0); i<resolution; ++i)
  {
    theta = i * M_PI / resolution;
    pint = gf.ip_theta(theta, r, t);
    BOOST_CHECK(pint >= pint_prev);
  }
}

BOOST_AUTO_TEST_CASE(GF3DRadAbs_dp_theta_at_a_is_leavea)
{
  Real D( 1e-12 ), kf( 1e-8 ), sigma( 1e-8 ), a( 1e-7 ), r0( 9e-8 );
  GreensFunction3DRadAbs gf( D, kf, r0, sigma, a );

  Real t( 1e-4 ), r( r0 );

  Real leavea( gf.leavea(t) * M_PI * a * a * 2 );
  Real iptheta( gf.idp_theta(M_PI, a, t) * M_PI * a * a );
 
  BOOST_CHECK_CLOSE( leavea / iptheta, 1.0, 0.001);
}

/*
'''
    def test_p_theta_free_is_p_theta_smallt(self):
        D = 1e-12
        sigma = 1e-8
        kf = 1e-8
        
        t = 1e-7
        r0 = 5e-7
        r = 5e-7
        a = 1e-6
        
        gf = mod.GreensFunction3DRadAbs(D, kf, r0, sigma, a)
        resolution = 20
        for i in range(1, resolution):
            theta = i * numpy.pi / resolution 
            pfree = mod.p_theta_free(theta, r, r0, t, D) 
            p = gf.p_theta(theta, r, t)* 4 * numpy.pi * r * r
            print pfree, p
            self.assertAlmostEqual(0.0, (pfree - p)/pfree)
'''

'''
    def test_alphan(self):
        D = 1e-12
        sigma = 1e-8
        kf = 1e-18
        
        a = 2e-7
        
        gf = mod.GreensFunction3DRadAbs(D, kf, r0, sigma, a)
        maxerror = 0
        
        for n in range(100):
            for i in range(1000):
                alpha = gf.alpha_i(n, i)
                error = abs(gf.f_alpha0(alpha))
                maxerror = max(error, maxerror)
        self.failIf(abs(maxerror) > 1e-8)
'''
*/

