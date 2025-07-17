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

#include "easyTime.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(easyTime, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        easyTime,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::easyTime::calcNu() const
{
    scalar timeIndex = U_.time().value();
    
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "nu",
                U_.time().timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            U_.mesh(),
            dimensionedScalar
            (
                "nu",
                dimViscosity,
                k_.value()*pow(timeIndex, n_.value())
            )
        )
    );
}

void Foam::viscosityModels::easyTime::correct()
{
    // 获取 alpha1 场
    const volScalarField& alpha1 = U_.mesh().lookupObject<volScalarField>("alpha.grout");
    
    nu_ = calcNu();
    nuDebug2_ = alpha1;  // 添加调试场2：alpha1场
  
    nuDebug_ = nu_;
    
    forAll(nuDebug_, celli)
    {
        nuDebug_[celli] = nu_[celli] * alpha1[celli];
    }
    
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::easyTime::easyTime
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    easyTimeCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
    k_("k", dimViscosity, easyTimeCoeffs_),
    n_("n", dimless, easyTimeCoeffs_),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE  // 改为 NO_WRITE不写入时间步
        ),
        calcNu()
    ),
    nuDebug_  // 添加调试场1：nu_field=alpha.grout*calcNu();
    (
        IOobject
        (
            name + "_field1",  // 例如 "nu1_field"
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        nu_
    ),
    nuDebug2_  // 修正：直接查找或使用占位符
    (
        IOobject
        (
            name + "alpha1",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh().lookupObject<volScalarField>("alpha.grout") // 方法2
    ) 
    
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



bool Foam::viscosityModels::easyTime::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    easyTimeCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    easyTimeCoeffs_.readEntry("k", k_);
    easyTimeCoeffs_.readEntry("n", n_);

    return true;
}


// ************************************************************************* //
