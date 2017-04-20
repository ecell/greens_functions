#include "compat.h"

#include <algorithm>
#include <stdexcept>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

//#include "Logger.hpp"
#include "freeFunctions.hpp"

namespace greens_functions
{

/**
   Calculates std::exp(x^2) * erfc(x)

   See asymptotic expansion here:
   http://en.wikipedia.org/wiki/Error_function
*/  
Real expxsq_erfc(Real x)
{
    Real result;

    const Real xsq(x * x);
    if(x > 26.0)
    {
        const Real M_1_SQRTPI(M_2_SQRTPI * 0.5); 

        const Real x2sq_r(1.0 / (2.0 * xsq));  // 2 / (2 x)^2

        /*
          up to second term in the expansion.
          abs err ~= 9e-8 at x == 20, 3e-8 at x == 25

          the third term 
          - (8 / (x2sq * x2sq * x2sq))       
          and beyond doesn't have a major contribution for large x.
        */

        result = (M_1_SQRTPI / x) * 
            (1.0 - x2sq_r +      // term 1
              x2sq_r * x2sq_r);  // term 2
    }
    else
    {
        result = std::exp(xsq) * erfc(x);
    }

    return result;
}


/**
   W(a, b) := std::exp(2 a b + b^2) erfc(a + b)
*/
Real W(Real a, Real b)
{
    // std::exp(2 a b + b^2) erfc(a + b) == 
    //               std::exp(- a^2) std::exp((a + b)^2) erfc(a + b)
    return std::exp(- a * a) * expxsq_erfc(a + b);
}


/* List of 1D functions used as approximations for the greensfunctions
   of the full domain. For small t, a particle only scans part of the domain
   and thus not all the boundary conditions have to be included. Because the
   greensfunctions converge bad for small t, we use these approximations.

   note: drift not included yet!

   *** Naming: 
   All 1D function names start with an X, followed by:

   P for position density function.
   S for survival probability function.
   I for culumative distribution function of position.

   ** Last two numer give boundaries 
   (By J.V. Beck - Heat Conduction in Green's Functions)

   0 - No boundary.
   1 - Absorbing boundary.
   2 - Reflective boundary.
   3 - Radiation boundary with rate ka.

   The approximation for the greensfunction with a sink has 3 boundaries: 030
   left absorbing boundary, sink, right absorbing boundary.
 */

Real XP00( Real r, Real t, Real r0, Real D, Real v )
{
    const Real fourDt( 4 * D * t );
    const Real rminr0minv2( gsl_pow_2( r - r0 - v * t) );

    return 1. / sqrt( fourDt * M_PI ) * exp( - rminr0minv2 / fourDt );
}


Real XI00( Real r, Real t, Real r0, Real D, Real v )
{
    const Real sqrt4Dt( sqrt( 4 * D * t ) );
    const Real rminr0minv( r - r0 - v * t);

    return 0.5 * ( 1.0 + erf( rminr0minv / sqrt4Dt ) );
}


Real XS00( Real t, Real r0, Real D, Real v )
{
    return 1.0;
}


Real XP10( Real r, Real t, Real r0, Real D, Real v )
{
    const Real fourDt( 4 * D * t );
    const Real rminr02( gsl_pow_2( r - r0 ) );
    const Real rplusr02( gsl_pow_2( r + r0 ) );
    const Real v_2( v /2.0 );

    Real drift_prefac;
    if( v == 0.0 )
        drift_prefac = 1.0;
    else
        drift_prefac = exp( v_2 / D * ( r - r0 - v_2 * t ) );   

    return drift_prefac / sqrt( fourDt * M_PI ) * 
        ( exp( - rminr02 / fourDt ) - exp( - rplusr02 / fourDt ) );
}


Real XI10( Real r, Real t, Real r0, Real D, Real v )
{
    const Real sqrt4Dt( sqrt( 4 * D * t ) );

    if( v == 0.0 )
    {
        const Real rminr0( r - r0 );
        const Real rplusr0( r + r0 );

        const Real temp( 2.0 * erf( r0 / sqrt4Dt ) 
                         + erf( rminr0/sqrt4Dt ) - erf( rplusr0/sqrt4Dt ) );

        return 0.5 * temp;
    }
    else
    {
        const Real r0minvt( r0 - v*t);
        const Real r0plusvt( r0 + v*t );

        const Real term1( erfc( (r0minvt + r)/sqrt4Dt ) - erfc( r0minvt/sqrt4Dt ) );
        const Real term2( erf( r0plusvt/sqrt4Dt ) - erf( (r0plusvt - r)/sqrt4Dt ) );

        return 0.5 * ( exp(- v * r0 / D ) * term1 + term2 );
    }
}


Real XS10( Real t, Real r0, Real D, Real v )
{
    const Real sqrt4Dt( sqrt( 4 * D * t ) );
    const Real vt( v * t );

    if( v == 0.0 )
        return erf( r0 / sqrt4Dt );
    else
    {
        const Real corr( exp((6*v*r0 - v*vt ) / (4 * D)) ); 
        //erfc( -(r0 + vt) / sqrt4Dt ) - exp(- v * r0 / D ) * erfc( (r0 - vt) / sqrt4Dt) );
        return 0.5 * corr * W(r0/sqrt4Dt, -vt/sqrt4Dt);
    }
}


Real XP20( Real r, Real t, Real r0, Real D, Real v )
{
    const Real fourDt( 4 * D * t );
    const Real sqrt4Dt( fourDt );
    const Real rminr02( gsl_pow_2( r - r0 ) );
    const Real rplusr02( gsl_pow_2( r + r0 ) );

    const Real XP20_nov( 1. / sqrt( fourDt * M_PI ) * 
                         ( exp( - rminr02 / fourDt ) + 
                           exp( - rplusr02 / fourDt ) ) );

    if( v == 0.0 )
    {
        return XP20_nov;
    }
    else
    {
        const Real v_2( v / 2.0 );
        const Real drift_prefac( exp( v_2 / D * ( r - r0 - v_2 * t ) ) );

        const Real XP20_v( exp(v_2 / D * (r + r0) + v_2 * v_2 / D * t) 
                           * erfc( (r + r0 + v * t)/sqrt4Dt ) );

        return drift_prefac * ( v_2 / D * XP20_v + XP20_nov );
    }
}


Real XI20( Real r, Real t, Real r0, Real D, Real v )
{
    const Real sqrt4Dt( sqrt( 4 * D * t ) );

    if( v == 0.0 )
    {
        const Real temp( erf( (r - r0)/sqrt4Dt ) + erf( (r + r0)/sqrt4Dt ) );
        return 0.5 * temp;
    }
    else
    {
        //TODO: c.d.f. with drift. It exist!
        return Real();
    }
}


Real XS20( Real t, Real r0, Real D, Real v )
{
    return 1.0;
}


Real XP30term_nov( Real r, Real t, Real r0, Real ka, Real D )
{
    r0 = fabs( r0 );
    const Real fourDt( 4 * D * t );
    const Real k_D( ka / D );
    const Real rplusr02( gsl_pow_2( r + r0 ) );
    const Real arg( (r + r0)/sqrt(fourDt) + ka * sqrt( t / D ) );

    return -k_D * exp( - rplusr02 / fourDt ) * expxsq_erfc( arg );
}


Real XP30term_v( Real r, Real t, Real r0, Real ka, Real D, Real v )
{
    r0 = fabs( r0 );
    const Real sqrt4Dt( sqrt( 4 * D * t ) );
    const Real v_2( v / 2.0 );
    const Real kplusv2( ka + v_2 );

    const Real erfc_arg( (r + r0 + 2 * kplusv2 * t ) / sqrt4Dt);

    const Real exp_arg( 1.0 / D * ( kplusv2 * kplusv2 * t + 
                                    kplusv2 * (r + r0) ) );

    return -exp( exp_arg ) * erfc( erfc_arg );
}


Real XP30( Real r, Real t, Real r0, Real ka, Real D, Real v )
{
    r0 = fabs( r0 );
    const Real XP20temp( XP20(r, t, r0, D, 0.0) );

    if( v == 0.0 )
        return XP20temp + XP30term_nov(r, t, r0, ka, D);
    else
    {
        const Real v_2( v / 2.0 );
        const Real drift_prefac( exp( v_2 / D * ( r - r0 - v_2 * t ) ) );
        
        return drift_prefac * ( XP20temp + 
                                1/D * (ka + v/2) * XP30term_v(r, t, r0, ka, D, v) );
    }
}


Real XI30term_nov( Real r, Real t, Real r0, Real ka, Real D )
{
    const Real sqrt4Dt( sqrt(4 * D * t) );
    r0 = fabs( r0 );
        
    const Real term1( erf( r0/sqrt4Dt ) - erf( (r+r0)/sqrt4Dt ) );

    const Real term2( W( r0/sqrt4Dt, 2*ka*t/sqrt4Dt ) );

    //exp( k_D * (ka * t + r0) ) * erfc( (2 * ka * t + r0)/sqrt4Dt ) );

    const Real term3( W( (r + r0)/sqrt4Dt, 2*ka*t/sqrt4Dt ) ); 

    //exp( k_D * (ka * t + r0 + r) ) * erfc( (2 * ka * t + r0 + r)/sqrt4Dt ) );

    return term1 + term2 - term3;
}


Real XI30( Real r, Real t, Real r0, Real ka, Real D, Real v )
{
    r0 = fabs( r0 );
    if( v == 0.0 )
        return XI20(r, t, r0, D, 0.0) + XI30term_nov(r, t, r0, ka, D);
    else
        //TODO: c.d.f with drift not included yet. Does exist!
        return Real();
}


Real XS30( Real t, Real r0, Real ka, Real D, Real v )
{
    r0 = fabs( r0 );
    const Real sqrt4Dt( sqrt( 4 * D * t ) );
    const Real k_D( ka / D );

    if( v == 0.0 )
        return erf( r0 / sqrt4Dt ) + W( r0 / sqrt4Dt, 2 * ka * t / sqrt4Dt );
    // exp( k_D * ka * t + k_D * r0 ) * erfc( (2 * ka * t + r0) / sqrt4Dt );
    else
    {
        const Real v_2( v / 2.0 );
        const Real kplusv2( ka + v_2 );
        const Real r0plusvt( r0 + v * t );
        const Real r0minvt( r0 - v * t );

        const Real erfc_arg( ( r0 + 2 * kplusv2 * t ) / sqrt4Dt );

        const Real exp_arg( k_D * ( r0 + (ka + v) * t ) );

        const Real term2( erfc( -r0plusvt / sqrt4Dt ) - 
                          ka / (ka + v) * exp(- v / D * r0) * erfc( r0minvt/sqrt4Dt ) );
        
        return (ka + v_2)/(ka + v) * exp( exp_arg ) * erfc( erfc_arg ) + 0.5 * term2;
    }
}


Real XP030( Real r, Real t, Real r0, Real ka, Real D )
{
    r0 = fabs( r0 );

    return XP00(r, t, r0, D, 0.0) + .5 * XP30term_nov(fabs(r), t, r0, .5 * ka, D);
}


Real XI030( Real r, Real t, Real r0, Real ka, Real D )
{
    ka *= .5;
    r0 = fabs( r0 );
    Real sign( 1 );

    if( r < 0 )
        sign *= -1;

    return XI00(r, t, r0, D, 0.0) + .5 * ( ( XS30(t, r0, ka, D, 0.0) - 1 ) + 
                                           sign * XI30term_nov( fabs(r), t, r0, ka, D ) );
}


Real XS030( Real t, Real r0, Real ka, Real D )
{
    return XS30( t, fabs( r0 ), 0.5 * ka, D, 0.0 );
}


Real __p_irr(Real r, Real t, Real r0, Real kf, Real D, Real sigma, Real alpha)
{
    //  printf("irrp %.16g %.16g %.16g\n",r,r0,t);
    const Real sqrtD(std::sqrt(D));

    const Real Dt4(4.0 * D * t);
    const Real r_plus_r0_minus_2sigma(r + r0 - 2.0 * sigma);

    const Real num1(std::exp(- gsl_pow_2(r - r0) / Dt4));
    const Real num2(std::exp(- gsl_pow_2(r_plus_r0_minus_2sigma) / Dt4));
    const Real num3(W(r_plus_r0_minus_2sigma / std::sqrt(Dt4), 
                        alpha * std::sqrt(t)));

    const Real num((num1 + num2) / std::sqrt(4.0 * M_PI * t) -  alpha * num3);

    const Real den(4.0 * M_PI * r * r0 * sqrtD);

    const Real result(num / den);

    const Real jacobian(4.0 * M_PI * r * r);

    return result * jacobian;
}

Real p_irr(Real r, Real t, Real r0, Real kf, Real D, Real sigma)
{
    const Real kD(4.0 * M_PI * sigma * D);
    const Real alpha((1.0 + (kf / kD)) * (std::sqrt(D) / sigma));

    const Real p(__p_irr(r, t, r0, kf, D, sigma, alpha));

    return p;
}


Real p_survival_irr(Real t, Real r0, Real kf, Real D, Real sigma)
{
    const Real kD(4.0 * M_PI * sigma * D);
    const Real alpha((1.0 + (kf / kD)) * (std::sqrt(D) / sigma));

    const Real p(__p_reaction_irr(t, r0, kf, D, sigma, alpha, kD));

    return 1.0 - p;
}

Real __p_reaction_irr(Real t, Real r0, Real kf, Real D, Real sigma,
                       Real alpha, Real kD)
{
    const Real sqrtt(std::sqrt(t));
    const Real sqrtD(std::sqrt(D));

    const Real r0_m_sigma_over_sqrt4D_t((r0 - sigma) 
                                         / ((sqrtD + sqrtD) * sqrtt));

    const Real Wf(W(r0_m_sigma_over_sqrt4D_t, alpha * sqrtt));
    const Real factor(sigma * kf / (r0 * (kf + kD)));

    return factor * (erfc(r0_m_sigma_over_sqrt4D_t) - Wf);
}


Real 
__p_reaction_irr_t_inf(Real r0, Real kf, Real sigma, Real kD)
{
    const Real kf_kD_r0((kf + kD) * r0);
    return 1 - (kf_kD_r0 - kf * sigma) / kf_kD_r0;
}


Real p_survival_nocollision(Real t, Real r0, Real D, Real a)
{
    const Real Dt(D * t);
    const Real asq(a * a);
    const Real a_r(1.0 / a);
    const Real asq_r(a_r * a_r);

    const Real PIr0(M_PI * r0);

    const Real angle_factor(PIr0 * a_r);
    const Real exp_factor(- Dt * M_PI * M_PI * asq_r);

    const Real TOLERANCE(1e-8);

    const unsigned int i_max(
        std::max(static_cast<unsigned int>(
                      std::ceil(std::sqrt(M_PI * M_PI 
                                  + asq * std::log(1.0 / TOLERANCE) / Dt) *
                            M_1_PI)), 2u));

    Real p(0.0);
    Real sign(1.0);
    unsigned int i(1);
    while(true)
    {
        const Real term(sign * 
                         std::exp(exp_factor * i * i) * 
                         std::sin(angle_factor * i) / i);
        
        p += term;

        if(i >= i_max)
        {
            break;
        }

        sign = -sign;
        ++i;
    }

    const Real factor((a + a) / PIr0);

    return p * factor;
}

Real dp_survival_nocollision(Real t, Real r0, Real D, Real a)
{
    const Real Dt(D * t);
    const Real asq(a * a);
    const Real a_r(1.0 / a);
    const Real asq_r(a_r * a_r);

    const Real PIr0(M_PI * r0);

    const Real angle_factor(PIr0 * a_r);
    const Real exp_factor(- Dt * M_PI * M_PI * asq_r);

    const Real TOLERANCE(1e-8);

    const unsigned int i_max(
        std::max(static_cast<unsigned int>(
                      std::ceil(std::sqrt(M_PI * M_PI 
                                  + asq * std::log(1.0 / TOLERANCE) / Dt) *
                            M_1_PI)), 2u));

    Real p(0.0);
    Real sign(- 1.0);
    unsigned int i(1);
    while(true)
    {
        const Real term(sign * 
                         std::exp(exp_factor * i * i) * 
                         std::sin(angle_factor * i) * i);
        
        p += term;

        if(i >= i_max)
        {
            break;
        }

        sign = -sign;
        ++i;
    }

    const Real factor(D * (M_PI + M_PI) / (a * r0));

    return p * factor;
}

Real p_theta_free(Real theta, Real r, Real r0, Real t, Real D)
{
    Real sin_theta;
    Real cos_theta;
    sincos(theta, &sin_theta, &cos_theta);

    const Real Dt4(4.0 * D * t);
    const Real Dt4Pi(Dt4 * M_PI);

    const Real term1(std::exp(- (r * r - 2.0 * cos_theta * r * r0 + r0 * r0) / 
                           Dt4));
    const Real term2(1.0 / std::sqrt(Dt4Pi * Dt4Pi * Dt4Pi));

    return term1 * term2 * sin_theta; // jacobian
}

Real ip_theta_free(Real theta, Real r, Real r0, Real t, Real D)
{
    const Real Dt(D * t);
    const Real Dt2(Dt + Dt);
    const Real rr0(r * r0);

    const Real rr0_over_2Dt(rr0 / Dt2);

    const Real rsqr0sq_over_4Dt((r * r + r0 * r0) / (Dt2 + Dt2));

    const Real term1(expm1(rr0_over_2Dt 
                             - rsqr0sq_over_4Dt));
    const Real term2(expm1(rr0_over_2Dt * cos(theta) 
                             - rsqr0sq_over_4Dt));

    const Real den(4.0 * std::sqrt(M_PI * M_PI * M_PI * Dt) * rr0);

    return (term1 - term2) / den;
}

Real g_bd_3D(Real r, Real sigma, Real t, Real D)
{
    const Real Dt4(4.0 * D * t);
    const Real mDt4_r(- 1.0 / Dt4);
    const Real sqrtDt4(std::sqrt(Dt4));
    const Real sqrtDt4_r(1.0 / sqrtDt4);
    const Real sqrtPi(std::sqrt(M_PI));

    const Real rps(r + sigma);
    const Real rms(r - sigma);

    const Real term1((std::exp(rps * rps * mDt4_r) - 
                        std::exp(rms * rms * mDt4_r)) * sqrtDt4 / 
                      (sqrtPi * r));
    const Real term2(erf(rps * sqrtDt4_r) - erf(rms * sqrtDt4_r));

    return 0.5 * (term1 + term2) * r * r;
}

Real g_bd_1D(Real r, Real sigma, Real t, Real D, Real v)
{
    const Real Dt4(4.0 * D * t);
    const Real sqrtDt4(std::sqrt(Dt4));
    const Real sqrtDt4_r(1.0 / sqrtDt4);
    const Real vt = v*t;

    const Real s_plus_r_plus_vt(sigma + r + vt);
    const Real s_min_r_min_vt(sigma - r - vt);

    const Real result(erfl(s_plus_r_plus_vt * sqrtDt4_r) + erfl(s_min_r_min_vt * sqrtDt4_r));

    return 0.5 * result;
}
    
Real I_bd_3D(Real sigma, Real t, Real D)
{
    const Real sqrtPi(std::sqrt(M_PI));

    const Real Dt(D * t);
    const Real Dt2(Dt + Dt);
    const Real sqrtDt(std::sqrt(Dt));
    const Real sigmasq(sigma * sigma);

    const Real term1(1.0 / (3.0 * sqrtPi));
    const Real term2(sigmasq - Dt2);
    const Real term3(Dt2 - 3.0 * sigmasq);
    const Real term4(sqrtPi * sigmasq * sigma * erfc(sigma / sqrtDt));

    const Real result(term1 * (- sqrtDt *
                                 (term2 * std::exp(- sigmasq / Dt) + term3)
                                 + term4));
    
    return result;
}

Real I_bd_1D(Real sigma, Real t, Real D, Real v)
{
    if(D == 0)
        return 0;

    const Real sqrtPi(std::sqrt(M_PI));

    const Real Dt4(4 * D * t);
    const Real sqrt4Dt(std::sqrt(Dt4));
    Real vt = v*t;

    const Real arg1(-(2*sigma + vt)*(2*sigma + vt)/Dt4);
    const Real term1(expl( -vt*vt/Dt4 ) - expl( arg1 ));
    const Real term2(vt*erfl( vt/sqrt4Dt ) - (2*sigma + vt)*erfl( (2*sigma + vt)/sqrt4Dt ));
    const Real result(1./2*(sqrt4Dt/sqrtPi*term1 + term2 + 2*sigma));

    return result;
    
}

Real I_bd_r_3D(Real r, Real sigma, Real t, Real D)
{
    const Real sqrtPi(std::sqrt(M_PI));

    const Real Dt(D * t);
    const Real Dt2(Dt + Dt);
    const Real Dt4(Dt2 + Dt2);
    const Real sqrtDt(std::sqrt(Dt));
    const Real sqrtDt4(std::sqrt(Dt4));
    const Real sigmasq(sigma * sigma);

    const Real sigmacb(sigmasq * sigma);
    const Real rcb(gsl_pow_3(r));

    const Real rsigma(r * sigma);

    const Real rps_sq(gsl_pow_2(r + sigma));
    const Real rms_sq(gsl_pow_2(r - sigma));

    const Real term1(- 2.0 * sqrtDt / sqrtPi);
    const Real term2(std::exp(- sigmasq / Dt) * (sigmasq - Dt2));
    const Real term3(- std::exp(- rps_sq / Dt4) * (rms_sq + rsigma - Dt2));
    const Real term4(std::exp(- rms_sq / Dt4) * (rps_sq - rsigma - Dt2));
    const Real term5(- sigmasq * 3.0 + Dt2);

    const Real term6((sigmacb - rcb) * erf((r - sigma) / sqrtDt4));
    const Real term7(- (sigmacb + sigmacb) * erf(sigma / sqrtDt));
    const Real term8((sigmacb + rcb) * erf((r + sigma) / sqrtDt4));

    const Real result((term1 * (term2 + term3 + term4 + term5)
                         // + sigmasq + rsigma + rsigma - Dt2)//expm1
                         + term6 + term7 + term8) / 6.0);
    
    return result;
}

Real I_bd_r_1D(Real r, Real sigma, Real t, Real D, Real v)
{
    if(D == 0)
        return 0;

    const Real sqrtPi(std::sqrt(M_PI));

    const Real Dt(D * t);
    const Real Dt2(Dt + Dt);
    const Real Dt4(Dt2 + Dt2);
    const Real sqrt4Dt(std::sqrt(Dt4));
    const Real vt = v*t;

    const Real smrmvt_sq(gsl_pow_2(sigma - r - vt));
    const Real sprpvt_sq(gsl_pow_2(sigma + r + vt));
    const Real twospvt_sq(gsl_pow_2(2*sigma + vt));

    const Real temp1(-expl( -smrmvt_sq/Dt4 ) + expl( -sprpvt_sq/Dt4 ));
    const Real temp2(expl( -vt*vt/Dt4 ) - expl( -twospvt_sq/Dt4 ));
    const Real term1(sqrt4Dt/sqrtPi*( temp1 + temp2 ));

    const Real term2(vt*erfl(vt/sqrt4Dt) - (2*sigma + vt)*erfl( (2*sigma + vt)/sqrt4Dt ));
    const Real term3((r - sigma + vt)*erfl( (sigma - r - vt)/sqrt4Dt ));
    const Real term4((r + sigma + vt)*erfl( (r + sigma + vt)/sqrt4Dt ));
    const Real result(1./2*(term1 + term2 + term3 + term4));
    
    return result;
}

struct g_bd_params
{ 
    const Real sigma;
    const Real t;
    const Real D;
    const Real target;
    const Real v;
};


static Real I_gbd_r_3D_F(Real r, const g_bd_params* params)
{
    const Real sigma(params->sigma);
    const Real t(params->t);
    const Real D(params->D);
    const Real target(params->target);

    return I_bd_r_3D(r, sigma, t, D) - target;
}

static Real I_gbd_r_1D_F(Real r, const g_bd_params* params)
{
    const Real sigma(params->sigma);
    const Real t(params->t);
    const Real D(params->D);
    const Real target(params->target);
    const Real v(params->v);

    return I_bd_r_1D(r, sigma, t, D, v) - target;
}


Real drawR_gbd_3D(Real rnd, Real sigma, Real t, Real D)
{
    const Real I(I_bd_3D(sigma, t, D));

    g_bd_params params = { sigma, t, D, rnd * I, 0 };

    gsl_function F =
    {
        reinterpret_cast<typeof(F.function)>(&I_gbd_r_3D_F),
        &params
    };

    Real low(sigma);
    Real high(sigma + 10.0 * std::sqrt (6.0 * D * t));

    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));
    gsl_root_fsolver_set(solver, &F, low, high);

    const unsigned int maxIter(100);

    unsigned int i(0);
    while(true)
    {
        gsl_root_fsolver_iterate(solver);

        low = gsl_root_fsolver_x_lower(solver);
        high = gsl_root_fsolver_x_upper(solver);
        int status(gsl_root_test_interval(low, high, 1e-18, 1e-12));

        if(status == GSL_CONTINUE)
        {
            if(i >= maxIter)
            {
                gsl_root_fsolver_free(solver);
                throw std::runtime_error("drawR_gbd: failed to converge");
            }
        }
        else
        {
            break;
        }

        ++i;
    }
  
    gsl_root_fsolver_free(solver);

    return low;
}

Real drawR_gbd_1D(Real rnd, Real sigma, Real t, Real D, Real v)
{
    const Real I(I_bd_1D(sigma, t, D, v));

    g_bd_params params = { sigma, t, D, rnd * I, v };

    gsl_function F =
    {
        reinterpret_cast<typeof(F.function)>(&I_gbd_r_1D_F),
        &params
    };

    Real low(sigma);
    Real high(sigma + 100.0 * std::sqrt (2.0 * D * t));

    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));
    gsl_root_fsolver_set(solver, &F, low, high);

    const unsigned int maxIter(100);

    unsigned int i(0);
    while(true)
    {
        gsl_root_fsolver_iterate(solver);

        low = gsl_root_fsolver_x_lower(solver);
        high = gsl_root_fsolver_x_upper(solver);
        int status(gsl_root_test_interval(low, high, 1e-18, 1e-12));

        if(status == GSL_CONTINUE)
        {
            if(i >= maxIter)
            {
                gsl_root_fsolver_free(solver);
                throw std::runtime_error("drawR_gbd: failed to converge");
            }
        }
        else
        {
            break;
        }

        ++i;
    }
  
    gsl_root_fsolver_free(solver);

    return low;
}

}
