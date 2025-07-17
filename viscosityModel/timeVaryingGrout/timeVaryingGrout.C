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
    tmp<volScalarField> tsr = strainRate();
    const volScalarField& sr = tsr();
    
    // 时间相关系数
    scalar timeIndex = sr.time().timeOutputValue();
    dimensionedScalar kEffective = min(nuMax_, k_*exp(timeCoeff_*timeIndex));
    
    // Papanastasiou 参数
    dimensionedScalar m("m", dimTime, 1000.0);  // 可调整

    scalar rhoValue = rho_.value();  // 获取密度值

    // 创建粘度场
    tmp<volScalarField> tnu
    (
        new volScalarField
        (
            IOobject
            (
                "nu",
                sr.time().timeName(),
                sr.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sr.mesh(),
            dimensionedScalar("nu", dimViscosity, 0)
        )
    );
    
    volScalarField& nu = tnu.ref();
    
    // 计算每个单元的粘度
    forAll(nu, celli)
    {
        scalar sri = sr[celli];
        scalar srMag = mag(sri);
        
        if (srMag > VSMALL)
        {
            // Papanastasiou 正则化
            scalar expTerm = exp(-m.value()*srMag);
            scalar reg = (1.0 - expTerm)/srMag;
            
            // 屈服应力项（正则化）
            scalar yieldTerm = tau0_.value() * reg;
            
            // 粘性项
            scalar viscTerm = kEffective.value() * pow(srMag, n_.value() - 1.0);
            
            // 总粘度
            nu[celli] = (yieldTerm + viscTerm)/2*rhoValue;
        }
        else
        {
            // 极限情况：剪切率为零
            nu[celli] = tau0_.value() * m.value()/2*rhoValue;
        }
        
        // // 应用上限
        // nu[celli] = min(nu[celli], nuMax_.value());
    }
    
    // 更新边界条件
    nu.correctBoundaryConditions();

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
        Debug1_.primitiveFieldRef() = mask.primitiveField() * sr.primitiveField();
        Debug2_.primitiveFieldRef() = mask.primitiveField() * kEffective.value();
        Debug3_.primitiveFieldRef() = mask.primitiveField() * nu.primitiveField();
        
        // 应用掩码到边界场
         forAll(Debug1_.boundaryFieldRef(), patchi)
    {
        Debug1_.boundaryFieldRef()[patchi] = 
            mask.boundaryField()[patchi] * sr.boundaryField()[patchi];
        
        Debug2_.boundaryFieldRef()[patchi] = 
            mask.boundaryField()[patchi] * kEffective.value();
        
        Debug3_.boundaryFieldRef()[patchi] = 
            mask.boundaryField()[patchi] * nu.boundaryField()[patchi];
    }
    
    // 输出调试信息
    Info<< "timeVaryingGrout debug information:" << endl;
    Info<< "    Current time: " << timeIndex << " s" << endl;
    Info<< "    kEffective: " << kEffective << " m2/s" << endl;
    Info<< "    Debug1 (strain rate) range: " << min(Debug1_).value() 
        << " to " << max(Debug1_).value() << " 1/s" << endl;
    Info<< "    Debug2 (kEffective) range: " << min(Debug2_).value() 
        << " to " << max(Debug2_).value() << " m2/s" << endl;
    Info<< "    Debug3 (nu) range: " << min(Debug3_).value() 
        << " to " << max(Debug3_).value() << " m2/s" << endl;
    Info<< "    Number of cells with alpha.grout > " << alphaSmall.value() 
        << ": " << returnReduce(sum(mask.primitiveField()), sumOp<scalar>()) << endl;
    
    // 额外的调试信息：检查正则化效果
    scalar nu0 = tau0_.value() * 1000.0;  // m = 1000
    Info<< "    Theoretical nu at sr=0: " << nu0 << " m2/s" << endl;
    Info<< "    nuMax limit: " << nuMax_.value() << " m2/s" << endl;
    
    // 清理临时对象

}
else
{
    Warning<< "    alpha.grout not found!" << endl;
    // 如果找不到 alpha.grout，将调试场设为0
    Debug1_ = dimensionedScalar("zero", Debug1_.dimensions(), 0.0);
    Debug2_ = dimensionedScalar("zero", Debug2_.dimensions(), 0.0);
    Debug3_ = dimensionedScalar("zero", Debug3_.dimensions(), 0.0);
}
    
    return tnu;
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
    rho_("rho", dimDensity, viscosityProperties),
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
            "physicalNu_Debug3", // 修改名称以反映实际内容
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
