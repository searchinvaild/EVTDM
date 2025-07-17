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
#include "fvCFD.H"
#include "timeVaryingGrout.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(timeVaryingGrout, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        timeVaryingGrout,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::timeVaryingGrout::calcNu() const
{
    Info<< "timeVaryingGrout::calcNu() called at time = " << U_.time().value() << endl;
    
    scalar timeIndex = U_.time().value();

    dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);

    tmp<volScalarField> sr = strainRate();
    volScalarField srLimited = max(sr(), dimensionedScalar("VSMALL", dimless/dimTime, VSMALL));


    // // 临时测试：直接计算应变率
    // volScalarField srDirect = sqrt(2.0)*mag(symm(fvc::grad(U_)));
    // Info<< "    Direct strain rate range: " << min(srDirect).value() << " to " << max(srDirect).value() << endl;
    // Info<< "    strainRate() range: " << min(sr()).value() << " to " << max(sr()).value() << endl;

    // // 比较两者的差异
    // volScalarField srDiff = mag(srDirect - sr());
    // Info<< "    Difference range: " << min(srDiff).value() << " to " << max(srDiff).value() << endl;

    tmp<volScalarField> nuCalc = min(
        nuMax_,
        (tau0_ + k_*exp(timeCoeff_*timeIndex)*rtone*pow(tone*sr(), n_)) / srLimited
    );

    // 更新调试场
    const volScalarField* alpha1Ptr = U_.mesh().findObject<volScalarField>("alpha.grout");
    if (alpha1Ptr)
    {
        const volScalarField& alpha1 = *alpha1Ptr;
       // 定义 VSMALL 阈值
        dimensionedScalar alphaSmall("alphaSmall", dimless, 0.9);
        
        // 创建掩码场：alpha.grout > VSMALL 时为1，否则为0
        volScalarField mask = pos(alpha1 - alphaSmall);
        
        // 应用掩码到内部场
        Debug1_.primitiveFieldRef() = mask.primitiveField() * srLimited.primitiveField();
        Debug2_.primitiveFieldRef() = mask.primitiveField() * nuCalc().primitiveField();
        Debug3_.primitiveFieldRef() = mask.primitiveField() * nuCalc().primitiveField() * srLimited.primitiveField();
        
        // 应用掩码到边界场
        forAll(Debug1_.boundaryFieldRef(), patchi)
        {
            Debug1_.boundaryFieldRef()[patchi] = 
                mask.boundaryField()[patchi] * srLimited.boundaryField()[patchi];
            Debug2_.boundaryFieldRef()[patchi] = 
                mask.boundaryField()[patchi] * nuCalc().boundaryField()[patchi];
            Debug3_.boundaryFieldRef()[patchi] = 
                mask.boundaryField()[patchi] * nuCalc().boundaryField()[patchi]* srLimited.boundaryField()[patchi];
        }
        
        // 确保边界条件正确更新
        Debug1_.correctBoundaryConditions();
        Debug2_.correctBoundaryConditions();
        Debug3_.correctBoundaryConditions();
        
        Info<< "    Debug1 range: " << min(Debug1_).value() << " to " << max(Debug1_).value() << endl;
        Info<< "    Debug2 range: " << min(Debug2_).value() << " to " << max(Debug2_).value() << endl;
        Info<< "    Debug3 range: " << min(Debug3_).value() << " to " << max(Debug3_).value() << endl;
        Info<< "    Number of cells with alpha.grout > VSMALL: " 
            << returnReduce(sum(mask.primitiveField()), sumOp<scalar>()) << endl;
    }
    else
    {
        Warning<< "    alpha.grout not found!" << endl;
        // 如果找不到 alpha.grout，将调试场设为0
        Debug1_ = dimensionedScalar("zero", Debug1_.dimensions(), 0.0);
        Debug2_ = dimensionedScalar("zero", Debug2_.dimensions(), 0.0);
        Debug3_ = dimensionedScalar("zero", Debug3_.dimensions(), 0.0);
    }

    return nuCalc;
}

void Foam::viscosityModels::timeVaryingGrout::correct()
{
    Info<< "timeVaryingGrout::correct() called at time = " << U_.time().value() << endl;
    
    // 更新 nu_
    nu_ = calcNu();
    
    // 强制写出调试场（用于测试）
    if (U_.time().outputTime())
    {
        Info<< "Writing timeVaryingGrout debug fields" << endl;
        Debug1_.write();
        Debug2_.write();
        Debug3_.write();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::timeVaryingGrout::timeVaryingGrout
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    timeVaryingGroutCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
    k_("k", dimViscosity, timeVaryingGroutCoeffs_),
    n_("n", dimless, timeVaryingGroutCoeffs_),
    tau0_("tau0", dimViscosity/dimTime, timeVaryingGroutCoeffs_),
    nuMax_("nuMax", dimViscosity, timeVaryingGroutCoeffs_),
    timeCoeff_("timeCoeff", dimless, timeVaryingGroutCoeffs_),
    Debug1_
    (
        IOobject
        (
            "StrainRate_Debug1",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("zero", dimless/dimTime, 0.0),
        "zeroGradient"  // 设置边界条件类型
    ),
    Debug2_
    (
        IOobject
        (
            "CaluNu_Debug2",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("zero", dimViscosity, 0.0),
        "zeroGradient"  // 设置边界条件类型
    ),
    Debug3_
    (
        IOobject
        (
            "Nu_Debug3",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("zero", dimViscosity, 0.0),
        "zeroGradient"  // 设置边界条件类型
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
{
    Info<< "timeVaryingGrout constructor: Created for phase " << name << endl;
    Info<< "    k = " << k_.value() << endl;
    Info<< "    n = " << n_.value() << endl;
    Info<< "    tau0 = " << tau0_.value() << endl;
    Info<< "    nuMax = " << nuMax_.value() << endl;
    Info<< "    timeCoeff = " << timeCoeff_.value() << endl;
    
    // 初始调用 correct
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::timeVaryingGrout::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    timeVaryingGroutCoeffs_ =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    timeVaryingGroutCoeffs_.readEntry("k", k_);
    timeVaryingGroutCoeffs_.readEntry("n", n_);
    timeVaryingGroutCoeffs_.readEntry("tau0", tau0_);
    timeVaryingGroutCoeffs_.readEntry("nuMax", nuMax_);
    timeVaryingGroutCoeffs_.readEntry("timeCoeff", timeCoeff_);

    return true;
}


// ************************************************************************* //
