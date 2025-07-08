/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "timeSlurryPower.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
// 添加在#include之后


namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(timeSlurryPower, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        timeSlurryPower,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::timeSlurryPower::calcNu() const
{
    scalar timeIndex = readScalar(U_.time().timeName());

    dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);

    tmp<volScalarField> sr(strainRate());

    return
    (
        min
        (
            nuMax_,
            (
                tau0_
              + k_*pow(timeIndex,timeCoeff_)*rtone*pow(tone*sr(), n_)
            )
           /(max(sr(), dimensionedScalar ("VSMALL", dimless/dimTime, VSMALL)))
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::timeSlurryPower::timeSlurryPower
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    timeSlurryCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
    k_("k", dimViscosity, timeSlurryCoeffs_),
    n_("n", dimless, timeSlurryCoeffs_),
    tau0_("tau0", dimViscosity/dimTime, timeSlurryCoeffs_),
    nuMax_("nuMax", dimViscosity, timeSlurryCoeffs_),
    timeCoeff_("timeCoeff", dimless, timeSlurryCoeffs_),
    nu_
    (
        IOobject
        (
            "nu_sludge",  // 强制指定字段名称（原name参数）
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::timeSlurryPower::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    timeSlurryCoeffs_ =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    timeSlurryCoeffs_.readEntry("k", k_);
    timeSlurryCoeffs_.readEntry("n", n_);
    timeSlurryCoeffs_.readEntry("tau0", tau0_);
    timeSlurryCoeffs_.readEntry("nuMax", nuMax_);
    timeSlurryCoeffs_.readEntry("timeCoeff", timeCoeff_);

    return true;
}


// ************************************************************************* //
