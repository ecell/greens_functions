#if !defined( __GREENSFUNCTION_HPP )
#define __GREENSFUNCTION_HPP

#include "Defs.hpp"
namespace greens_functions{

enum EventKind
    {
        IV_ESCAPE,
        IV_REACTION
    };

class GreensFunction
{
public:
    GreensFunction( const Real D )
      : D( D ) {}
  
    ~GreensFunction() {}
  
    Real getD() const
    {
        return this->D;
    }

protected:
    const Real D;
};

}
#endif // __GREENSFUNCTION_HPP
