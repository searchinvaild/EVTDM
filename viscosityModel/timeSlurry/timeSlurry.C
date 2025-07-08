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
    defineTypeNameAndDebug(timeSlurry, 0); //set debug level to 0 for no debug output

    addToRunTimeSelectionTable //运行时选择表:菜单系统，告诉openFoam，我们的菜单有了新品：timeslurry牛排
    (
        viscosityModel, // 基类: 新品timeslurry是西餐（viscosity）
        timeSlurry,     // 新品名称: timeslurry
        dictionary // 读取字典: 新品牛排几成熟客人可以在字典文件中指定
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>  //建立一个临时的体积标量场
Foam::viscosityModels::timeSlurry::calcNu() const // calcNu()是函数名，前面都是在告诉calcNu()属于Foam商店，viscosityModels餐别，timeSlurry牛排
{
    scalar timeIndex = U_.time().value(); // 这里使用U_.time()或者p_.time()

    dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);

   // tmp<volScalarField> sr(strainRate()); //src/transportModels/incompressible/viscosityModels/viscosityModel/viscosityModel.C：
    // return sqrt(2.0)*mag(symm(fvc::grad(U_)));
    tmp<volScalarField> sr = strainRate();
    volScalarField srLimited = max(sr(), dimensionedScalar("VSMALL", dimless/dimTime, VSMALL));

    return
    min(
        nuMax_,
        (
            tau0_
            + k_*exp(timeCoeff_*timeIndex)*rtone*pow(tone*sr(), n_)
        )
        / srLimited  // 使用处理后的应变率
    );
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
    nu_
    (
        IOobject
        (
            name,  // 强制指定字段名称（原name参数）
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


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
