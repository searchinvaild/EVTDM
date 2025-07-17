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

#include "timeSlurry.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(timeSlurry, 0);
    addToRunTimeSelectionTable
    (
        viscosityModel,
        timeSlurry,
        dictionary
    );
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::timeSlurry::calcNu() const
{
    scalar timeIndex = U_.time().value();

    dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);

    tmp<volScalarField> sr = strainRate();
    volScalarField srLimited = max(sr(), dimensionedScalar("VSMALL", dimless/dimTime, VSMALL));

    return
    min(
        nuMax_,
        (
            tau0_
            + k_*exp(timeCoeff_*timeIndex)*rtone*pow(tone*sr(), n_)
        )
        / srLimited
    );
}
void Foam::viscosityModels::timeSlurry::correct()
{
    // 获取剪切应变场
    tmp<volScalarField> sr = strainRate();
    volScalarField srLimited = max(sr(), dimensionedScalar("VSMALL", dimless/dimTime, VSMALL));

    // 安全获取 alpha.grout 场
    if (!U_.mesh().foundObject<volScalarField>("alpha.grout")) {
        FatalErrorInFunction
            << "Required field alpha.grout not found in mesh registry!"
            << exit(FatalError);
    }
    const volScalarField& alpha1 = U_.mesh().lookupObject<volScalarField>("alpha.grout");
    
    // 获取计算粘度场
    volScalarField nuCalc = calcNu(); 
    
    // 应用条件显示：仅在 alpha.grout > 0 的区域显示值
    Debug1_ = srLimited * pos(alpha1);  // pos(x) = (x > 0) ? 1 : 0
    Debug2_ = nuCalc * pos(alpha1);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::timeSlurry::timeSlurry
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
    Debug1_
    (
        IOobject
        (
            name + "_Debug1",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("zero", dimless/dimTime, 0.0)
    ),
    Debug2_
    (
        IOobject
        (
            name + "_Debug2",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("zero", dimViscosity, 0.0)
    ),
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
{} // 修复：移除了多余的括号


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



bool Foam::viscosityModels::timeSlurry::read
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
