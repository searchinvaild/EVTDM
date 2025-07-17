/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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
    Foam::viscosityModels::timeSlurry

Description
     Herschel-Bulkley non-Newtonian viscosity model.

SourceFiles
    timeSlurry.C

\*---------------------------------------------------------------------------*/

#ifndef timeSlurry_H
#define timeSlurry_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class timeSlurry Declaration
\*---------------------------------------------------------------------------*/

class timeSlurry
:
    public viscosityModel
{
    // Private data

        dictionary timeSlurryCoeffs_;

        dimensionedScalar k_;
        dimensionedScalar n_;
        dimensionedScalar tau0_;
        dimensionedScalar nuMax_;
        dimensionedScalar timeCoeff_;
        volScalarField nu_;
        tmp<volScalarField> calcNu() const;

    // Debug fields
        mutable volScalarField Debug1_;
        mutable volScalarField Debug2_;
        mutable bool needsUpdate_;  // 添加更新标志
protected:
    mutable volScalarField nu_;  // 添加 mutable

    

public:

    //- Runtime type information
    TypeName("timeSlurry");


    // Constructors

        //- Construct from components
        timeSlurry
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~timeSlurry() = default;


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
        const_cast<timeSlurry*>(this)->correct();
        return nu_;
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Correct the laminar viscosity
        // virtual void correct()
        // {
        //     nu_ = calcNu();
        // }
        virtual void correct();
        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
