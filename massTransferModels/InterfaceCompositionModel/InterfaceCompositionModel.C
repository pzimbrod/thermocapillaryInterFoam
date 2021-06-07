/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "InterfaceCompositionModel.H"
#include "phaseModel.H"
#include "phasePair.H"
#include "pureMixture.H"
#include "multiComponentMixture.H"
#include "rhoThermo.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
template<class ThermoType>
const typename Foam::multiComponentMixture<ThermoType>::thermoType&
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::getLocalThermo
(
    const word& speciesName,
    const multiComponentMixture<ThermoType>& globalThermo
) const
{
    return
        globalThermo.getLocalThermo
        (
            globalThermo.species()
            [
                speciesName
            ]
        );
}


template<class Thermo, class OtherThermo>
template<class ThermoType>
const typename Foam::pureMixture<ThermoType>::thermoType&
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::getLocalThermo
(
    const word& speciesName,
    const pureMixture<ThermoType>& globalThermo
) const
{
    return globalThermo.cellMixture(0);
}


template<class Thermo, class OtherThermo>
template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::getSpecieMassFraction
(
    const word& speciesName,
    const multiComponentMixture<ThermoType>& mixture
) const
{
    const fvMesh& mesh = fromThermo_.p().mesh();

    auto tY = tmp<volScalarField>::New
    (
        IOobject
        (
            "tY",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    );

    auto& Ys = tY.ref();

    Ys = mixture.Y(speciesName);

    return tY;
}


template<class Thermo, class OtherThermo>
template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::getSpecieMassFraction
(
    const word& speciesName,
    const pureMixture<ThermoType>& mixture
) const
{
    const fvMesh& mesh = fromThermo_.p().mesh();

    return tmp<volScalarField>::New
    (
        IOobject
        (
            "tY",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, scalar(1)),
        zeroGradientFvPatchScalarField::typeName
    );
}


template<class Thermo, class OtherThermo>
template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::MwMixture
(
    const pureMixture<ThermoType>& mixture
) const
{
    const fvMesh& mesh = fromThermo_.p().mesh();

    return tmp<volScalarField>::New
    (
        IOobject
        (
            "tM",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "Mw",
            dimMass/dimMoles,
            1e-3*mixture.cellMixture(0).W()
        ),
        zeroGradientFvPatchScalarField::typeName
    );
}


template<class Thermo, class OtherThermo>
template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::MwMixture
(
    const multiComponentMixture<ThermoType>& mixture
) const
{
    return refCast<const basicSpecieMixture>(mixture).W();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::InterfaceCompositionModel
(
    const dictionary& dict,
    const phasePair& pair
)
:
    interfaceCompositionModel(dict, pair),
    fromThermo_
    (
        pair.from().mesh().lookupObject<Thermo>
        (
            IOobject::groupName
            (
                basicThermo::dictName,
                pair.from().name()
            )
        )
    ),
    toThermo_
    (
        pair.to().mesh().lookupObject<OtherThermo>
        (
            IOobject::groupName
            (
                basicThermo::dictName,
                pair.to().name()
            )
        )
    ),
    Le_("Le", dimless, 1.0, dict)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::D
(
    const word& speciesName
) const
{
    const typename Thermo::thermoType& fromThermo =
        getLocalThermo
        (
            speciesName,
            fromThermo_
        );

    const volScalarField& p(fromThermo_.p());

    const volScalarField& T(fromThermo_.T());

    auto tmpD = tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("D", pair_.name()),
            p.time().timeName(),
            p.mesh()
        ),
        p.mesh(),
        dimensionedScalar(dimArea/dimTime, Zero)
    );

    auto& D = tmpD.ref();

    forAll(p, cellI)
    {
        D[cellI] =
            fromThermo.alphah(p[cellI], T[cellI])
           /fromThermo.rho(p[cellI], T[cellI]);
    }

    D /= Le_;
    D.correctBoundaryConditions();

    return tmpD;
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::L
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    const typename Thermo::thermoType& fromThermo =
        getLocalThermo(speciesName, fromThermo_);
    const typename OtherThermo::thermoType& toThermo =
        getLocalThermo(speciesName, toThermo_);

    const volScalarField& p(fromThermo_.p());

    auto tmpL = tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("L", pair_.name()),
            p.time().timeName(),
            p.mesh()
        ),
        p.mesh(),
        dimensionedScalar(dimEnergy/dimMass, Zero),
        zeroGradientFvPatchScalarField::typeName
    );

    auto& L = tmpL.ref();

    // from Thermo (from) to Thermo (to)
    forAll(p, cellI)
    {
        L[cellI] = fromThermo.Hc() - toThermo.Hc();
    }

    L.correctBoundaryConditions();

    return tmpL;
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::dY
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    NotImplemented;
    return nullptr;
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::Yf
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    NotImplemented;
    return nullptr;
}


// ************************************************************************* //
