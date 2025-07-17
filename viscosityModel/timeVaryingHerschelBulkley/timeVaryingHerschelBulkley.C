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

\*---------------------------------------------------------------------------*/

#include "timeVaryingHerschelBulkley.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(timeVaryingHerschelBulkley, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        timeVaryingHerschelBulkley,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::dimensionedScalar
Foam::viscosityModels::timeVaryingHerschelBulkley::calcK(const scalar t) const
{
    if (timeVariationType_ == "power")
    {
        // k = A * t^B
        if (t > SMALL)
        {
            return A_ * pow(t, B_.value());
        }
        else
        {
            return k0_;
        }
    }
    else if (timeVariationType_ == "exponential")
    {
        // k = A * exp(B*t)
        return A_ * exp(B_.value() * t);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown time variation type: " << timeVariationType_
            << nl << "Valid types are: power, exponential"
            << exit(FatalError);
        
        return k0_;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::timeVaryingHerschelBulkley::calcNu() const
{
    // Get current time
    const scalar t = U_.time().timeOutputValue();
    
    // Calculate time-varying k
    dimensionedScalar k = calcK(t);
    
    Info<< "Time = " << t << " s, k = " << k.value() << " Pa.s^n" << endl;

    // Calculate strain rate magnitude
    tmp<volScalarField> tsr = strainRate();
    const volScalarField& sr = tsr();

    // Herschel-Bulkley model implementation
    return min
    (
        nuMax_,
        max
        (
            nuMin_,
            (tau0_ + k*pow(sr, n_))/(max(sr, dimensionedScalar("SMALL", dimless/dimTime, SMALL)))
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::timeVaryingHerschelBulkley::timeVaryingHerschelBulkley
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    timeVaryingHerschelBulkleyCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
    k0_("k0", dimViscosity*pow(dimTime, dimensionedScalar("n", dimless, 1)), timeVaryingHerschelBulkleyCoeffs_),
    n_("n", dimless, timeVaryingHerschelBulkleyCoeffs_),
    tau0_("tau0", dimViscosity/dimTime, timeVaryingHerschelBulkleyCoeffs_),
    nuMin_("nuMin", dimViscosity, timeVaryingHerschelBulkleyCoeffs_),
    nuMax_("nuMax", dimViscosity, timeVaryingHerschelBulkleyCoeffs_),
    A_("A", dimViscosity*pow(dimTime, dimensionedScalar("n", dimless, 1)), timeVaryingHerschelBulkleyCoeffs_),
    B_("B", dimless, timeVaryingHerschelBulkleyCoeffs_),
    timeVariationType_(timeVaryingHerschelBulkleyCoeffs_.lookup("timeVariationType")),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{
    // Validate time variation type
    if (timeVariationType_ != "power" && timeVariationType_ != "exponential")
    {
        FatalErrorInFunction
            << "Invalid time variation type: " << timeVariationType_
            << nl << "Valid types are: power, exponential"
            << exit(FatalError);
    }

    Info<< "Selecting time-varying Herschel-Bulkley model" << nl
        << "    Time variation type: " << timeVariationType_ << nl
        << "    Initial k0 = " << k0_.value() << nl
        << "    A = " << A_.value() << nl
        << "    B = " << B_.value() << nl
        << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::timeVaryingHerschelBulkley::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    timeVaryingHerschelBulkleyCoeffs_ = 
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    timeVaryingHerschelBulkleyCoeffs_.lookup("k0") >> k0_;
    timeVaryingHerschelBulkleyCoeffs_.lookup("n") >> n_;
    timeVaryingHerschelBulkleyCoeffs_.lookup("tau0") >> tau0_;
    timeVaryingHerschelBulkleyCoeffs_.lookup("nuMin") >> nuMin_;
    timeVaryingHerschelBulkleyCoeffs_.lookup("nuMax") >> nuMax_;
    timeVaryingHerschelBulkleyCoeffs_.lookup("A") >> A_;
    timeVaryingHerschelBulkleyCoeffs_.lookup("B") >> B_;
    timeVaryingHerschelBulkleyCoeffs_.lookup("timeVariationType") >> timeVariationType_;

    return true;
}


// ************************************************************************* //
