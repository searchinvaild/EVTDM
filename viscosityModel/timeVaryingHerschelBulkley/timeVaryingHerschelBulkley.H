/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::viscosityModels::timeVaryingHerschelBulkley

Description
    Time-varying Herschel-Bulkley non-Newtonian viscosity model.
    
    The consistency coefficient k varies with time according to:
    - Power law: k = A * t^B
    - Exponential: k = A * exp(B*t)

SourceFiles
    timeVaryingHerschelBulkley.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingHerschelBulkley_H
#define timeVaryingHerschelBulkley_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                Class timeVaryingHerschelBulkley Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingHerschelBulkley
:
    public viscosityModel
{
    // Private data

        //- Dictionary
        dictionary timeVaryingHerschelBulkleyCoeffs_;

        //- Initial consistency coefficient
        dimensionedScalar k0_;

        //- Flow behaviour index
        dimensionedScalar n_;

        //- Yield stress
        dimensionedScalar tau0_;

        //- Minimum viscosity
        dimensionedScalar nuMin_;

        //- Maximum viscosity
        dimensionedScalar nuMax_;

        //- Time variation coefficient A
        dimensionedScalar A_;

        //- Time variation coefficient B
        dimensionedScalar B_;

        //- Time variation type (power or exponential)
        word timeVariationType_;

        //- Current viscosity field
        volScalarField nu_;


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;

        //- Calculate time-varying k
        dimensionedScalar calcK(const scalar t) const;


public:

    //- Runtime type information
    TypeName("timeVaryingHerschelBulkley");


    // Constructors

        //- Construct from components
        timeVaryingHerschelBulkley
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~timeVaryingHerschelBulkley()
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
        virtual void correct()
        {
            nu_ = calcNu();
        }

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
