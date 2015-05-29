#define BOOST_TEST_MODULE "GreensFunction3DSym_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include "../GreensFunction3DSym.hpp"

using namespace greens_functions;

BOOST_AUTO_TEST_CASE(GF3DSym_test_constructor)
{
  Real D(1e-12);
  GreensFunction3DSym gf(D);

}

BOOST_AUTO_TEST_CASE(GF3DSym_drawR)
{
  Real D(1e-12);
  GreensFunction3DSym gf(D);

  Real r, t(1e-3);

  r = gf.drawR(0.5, t);
  BOOST_CHECK((r >= 0.0));

  r = gf.drawR(0.0, t);
  BOOST_CHECK((r >= 0.0));

  r = gf.drawR(1.0, t);
  BOOST_CHECK((r >= 0.0));
}

BOOST_AUTO_TEST_CASE(GF3DSym_drawR_zero_t_is_zero)
{
  Real D(1e-12);
  GreensFunction3DSym gf(D);

  Real r, t(0.0);

  r = gf.drawR(0.5, t);
  BOOST_CHECK_EQUAL(r, 0.0);

  r = gf.drawR(0.0, t);
  BOOST_CHECK_EQUAL(r, 0.0);

  r = gf.drawR(1.0, t);
  BOOST_CHECK_EQUAL(r, 0.0);
}

BOOST_AUTO_TEST_CASE(GF3DSym_no_test_ip_r_infinity_is_one)
{
  Real D(1e-12);
  GreensFunction3DSym gf(D);

  Real r(2.5e-8), t(1e-5), ip;

  ip = gf.ip_r(INFINITY, t);
  BOOST_CHECK((ip == 1.0));
}
/*
def test_int_p_r_is_ip_r(self):

        import scipy.integrate

        D = 1e-12
        t = 1e-5
        
        gf = mod.GreensFunction3DSym(D)

        ip = gf.ip_r(0.0, t)
        self.assertEqual(0.0, ip)

        maxr = 5e-8

        resolution = 20
        for i in range(1, resolution):
            r = i * maxr / resolution 
            ip = gf.ip_r(r, t)
            result = scipy.integrate.quad(gf.p_r, 0.0, r,
                                          args=(t, ))
            np = result[0]
            self.assertAlmostEqual(0.0, (np-ip)/ip)
*/

