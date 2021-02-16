#ifndef GREENS_FUNCTIONS_3D_SYM_HPP
#define GREENS_FUNCTIONS_3D_SYM_HPP

#include "compat.h"

#include <gsl/gsl_integration.h>

//added
#include <string>
//#include "Logger.hpp"

#include "GreensFunction.hpp"

/**
  Green's Function for a free diffusion particle.
*/

namespace greens_functions{

class GreensFunction3DSym
    :
    public GreensFunction
{

private:

    static const Real TOLERANCE;

    static const Real H;

public:

    GreensFunction3DSym( const Real D )
        :
        GreensFunction( D )
    {
        ; // do nothing
    }


    ~GreensFunction3DSym()
    {
        ; // do nothing
    }

    Real drawTime( const Real ) const
    {
        return std::numeric_limits<Real>::infinity();
    }

    Real drawR( const Real rnd, const Real t ) const;

    Real p_r( const Real r, const Real t ) const;

    Real ip_r( const Real r, const Real t ) const;


    std::string dump() const;

    const char* getName() const
    {
        return "GreensFunction3DSym";
    }

private:

//    static Logger& log_;
};


}
#endif // __FREEGREENSFUNCTION
