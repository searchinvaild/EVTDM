#ifndef easyTime_H
#define easyTime_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

namespace Foam
{
namespace viscosityModels
{

class easyTime
:
    public viscosityModel
{
    // Private data
    dictionary easyTimeCoeffs_;
    dimensionedScalar k_;
    dimensionedScalar n_;
    volScalarField nu_;
    
    // 如果需要调试场
    mutable volScalarField nuDebug_;

    // Private Member Functions
    tmp<volScalarField> calcNu() const;

public:

    //- Runtime type information
    TypeName("easyTime");

    // Constructors
    easyTime
    (
        const word& name,
        const dictionary& viscosityProperties,
        const volVectorField& U,
        const surfaceScalarField& phi
    );

    //- Destructor
    virtual ~easyTime()
    {}

    // Member Functions

    //- Return the laminar viscosity
    virtual tmp<volScalarField> nu() const
    {
        return nu_;
    }

    //- Return the laminar viscosity for patch
    virtual tmp<scalarField> nu(const label patchi) const
    {
        return nu_.boundaryField()[patchi];
    }

    //- Correct the laminar viscosity
    virtual void correct();

    //- Read transportProperties dictionary
    virtual bool read(const dictionary& viscosityProperties);
};

} // End namespace viscosityModels
} // End namespace Foam

#endif
