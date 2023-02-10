/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "DynSmag.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class BasicMomentumTransportModel>
void DynSmag<BasicMomentumTransportModel>::correctNut
(
    const tmp<volTensorField>& gradU
)
{
    this->nut_ = max(Cs_*sqr(this->delta())*mag(dev(symm(gradU))),-1.0*this->nu());
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}

template<class BasicTurbulenceModel>
void DynSmag<BasicTurbulenceModel>::correctNut()
{
    correctNut(fvc::grad(this->U_));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
DynSmag<BasicMomentumTransportModel>::DynSmag
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    LESeddyViscosity<BasicMomentumTransportModel>(
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),  
    k_(
        IOobject
        (
            IOobject::groupName("k", this->U_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    Cs_(
        IOobject
        (
            IOobject::groupName("Cs", this->U_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar ("Cs", dimless,SMALL)
    ),
    simpleFilter_(U.mesh()),
    filterPtr_(LESfilter::New(U.mesh(), this->coeffDict())),
    filter_(filterPtr_()),
    y_(wallDist::New(this->mesh_).y())
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool DynSmag<BasicMomentumTransportModel>::read()
{
    if (LESeddyViscosity<BasicMomentumTransportModel>::read())
    {
        filter_.read(this->coeffDict());        
        return true;
    }
    else 
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void DynSmag<BasicMomentumTransportModel>::correct()
{
    const volVectorField& U = this->U_;

    LESeddyViscosity<BasicMomentumTransportModel>::correct();
    
    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    volSymmTensorField S(dev(symm(gradU)));
    volScalarField magS(mag(S));

    volVectorField Uf(filter_(U));

    volSymmTensorField Sf(filter_(S));  
    
    volScalarField magSf(mag(Sf));
          
    
    volSymmTensorField LL = (dev(filter_(sqr(U)) - (sqr(filter_(U)))));

    volSymmTensorField MM(sqr(this->delta())*(filter_(magS*S) - 4.0*magSf*Sf));
    
    volScalarField MMMM = fvc::average(magSqr(MM));
    MMMM.max(VSMALL);

    Cs_= 0.5* fvc::average(LL && MM)/MMMM;

    volScalarField KK =
    0.5*(filter_(magSqr(U)) - magSqr(filter_(U)));
    
    volScalarField mm(sqr(this->delta())*(4.0*sqr(mag(Sf)) - filter_(sqr(magS))));

    volScalarField mmmm = fvc::average(magSqr(mm));
    mmmm.max(VSMALL);

    k_ = fvc::average(KK*mm)/mmmm * sqr(this->delta())*magSqr(S);

    correctNut(gradU);

    volScalarField filter_width(this->delta());

    if(this->runTime_.outputTime())
    {
        const char* Cs_name = "Cs";
        Cs_.rename(Cs_name);
        Cs_.write();

        const char* Uf_name = "Uf";
        Uf.rename(Uf_name);
        Uf.write();

        const char* S_name = "S_ij";
        S.rename(S_name);
        S.write();

        const char* delta_name = "del";
        filter_width.rename(delta_name);
        filter_width.write();

        const char* y_name = "yw";
        y_.rename(y_name);
        y_.write();

        const char* nut_name = "nut";
        this->nut_.rename(nut_name);
        this->nut_.write();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //