/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "phaseSystem.H"
#include "surfaceTensionModel.H"
#include "porousModel.H"

#include "HashPtrTable.H"

#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "fvMatrix.H"

#include "fvcReconstruct.H"
#include "Identity.H"

#include "zeroGradientFvPatchFields.H"
#include "fixedEnergyFvPatchScalarField.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseSystem, 0);
}

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

const Foam::word Foam::phaseSystem::phasePropertiesName("phaseProperties");

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::phaseSystem::calcMu()
{
    mu_ = mu()();
}


Foam::phaseSystem::phaseModelTable
Foam::phaseSystem::generatePhaseModels(const wordList& phaseNames) const
{
    phaseModelTable phaseModels;

    for (const word& phaseName : phaseNames)
    {
        phaseModels.insert
        (
            phaseName,
            phaseModel::New
            (
                *this,
                phaseName
            )
        );
    }

    return phaseModels;
}


Foam::tmp<Foam::surfaceScalarField> Foam::phaseSystem::generatePhi
(
    const phaseModelTable& phaseModels
) const
{
    auto iter = phaseModels.cbegin();

    auto tmpPhi = tmp<surfaceScalarField>::New
    (
        "phi",
        fvc::interpolate(iter()()) * iter()->phi()
    );

    for (++iter; iter != phaseModels.cend(); ++iter)
    {
        tmpPhi.ref() += fvc::interpolate(iter()()) * iter()->phi();
    }

    return tmpPhi;
}


void Foam::phaseSystem::generatePairs(const dictTable& modelDicts)
{
    forAllConstIters(modelDicts, iter)
    {
        const phasePairKey& key = iter.key();

        // pair already exists
        if (phasePairs_.found(key))
        {
            // do nothing ...
        }
        else if (key.ordered())
        {
            // New ordered pair
            phasePairs_.insert
            (
                key,
                autoPtr<phasePair>
                (
                    new orderedPhasePair
                    (
                        phaseModels_[key.first()],
                        phaseModels_[key.second()]
                    )
                )
            );
        }
        else
        {
            // New unordered pair
            phasePairs_.insert
            (
                key,
                autoPtr<phasePair>
                (
                    new phasePair
                    (
                        phaseModels_[key.first()],
                        phaseModels_[key.second()]
                    )
                )
            );
        }
    }
}


void Foam::phaseSystem::generatePairsTable()
{
    forAllConstIters(phaseModels_, phaseIter1)
    {
        forAllConstIters(phaseModels_, phaseIter2)
        {
            if (phaseIter1()->name() != phaseIter2()->name())
            {
                phasePairKey key
                (
                    phaseIter1()->name(),
                    phaseIter2()->name(),
                    true
                );

                phasePairKey keyInverse
                (
                    phaseIter2()->name(),
                    phaseIter1()->name(),
                    true
                );

                if
                (
                    !totalPhasePairs_.found(key)
                 && !totalPhasePairs_.found(keyInverse)
                )
                {
                    totalPhasePairs_.set
                    (
                        key,
                        autoPtr<phasePair>
                        (
                            new phasePair
                            (
                                phaseModels_[key.first()],
                                phaseModels_[key.second()]
                            )
                        )
                    );
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseSystem::phaseSystem
(
    const fvMesh& mesh
)
:
    basicThermo(mesh, word::null, phasePropertiesName),
    mesh_(mesh),
    mu_
    (
        IOobject
        (
            "mu",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimViscosity*dimDensity, Zero),
        calculatedFvPatchScalarField::typeName
    ),
    phaseNames_(get<wordList>("phases")),
    phi_
    (
        IOobject
        (
            "phi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimVolume/dimTime, Zero)
    ),
    rhoPhi_
    (
        IOobject
        (
            "rhoPhi",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimMass/dimTime, Zero)
    ),
    phaseModels_(generatePhaseModels(phaseNames_)),
    phasePairs_(),
    totalPhasePairs_(),
    Prt_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Prt", *this, 1.0
        )
    )
{
    rhoPhi_.setOriented();
    phi_.setOriented();

    // sub models
    if (found("surfaceTension"))
    {
        generatePairsAndSubModels
        (
            "surfaceTension",
            mesh_,
            surfaceTensionModels_
        );
    }
    if (found("interfacePorous"))
    {
        generatePairsAndSubModels
        (
            "interfacePorous",
            mesh_,
            interfacePorousModelTable_
        );
    }

    // Total phase pair
    generatePairsTable();

    // Update mu_
    calcMu();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseSystem::~phaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::phaseSystem::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::hc() const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tAlphaHc
    (
        iter()() * iter()->hc()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tAlphaHc.ref() += iter()() * iter()->hc();
    }

    return tAlphaHc;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::THE
(
    const scalarField& e,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::THE
(
    const scalarField& e,
    const scalarField& p,
    const scalarField& T0,
    const label patchI
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::rho() const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tmpRho
    (
        iter()() * iter()->rho()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpRho.ref() += iter()() * iter()->rho();
    }

    return tmpRho;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::rho(const label patchI) const
{
    auto iter = phaseModels_.cbegin();

    tmp<scalarField> tmpRho
    (
        iter()().boundaryField()[patchI]
      * iter()->rho()().boundaryField()[patchI]
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpRho.ref() +=
        (
            iter()().boundaryField()[patchI]
          * iter()->rho()().boundaryField()[patchI]
        );
    }

    return tmpRho;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::Cp() const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tmpCp
    (
        iter()() * iter()->Cp()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpCp.ref() += iter()() * iter()->Cp();
    }

    return tmpCp;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    auto iter = phaseModels_.cbegin();

    tmp<scalarField> tmpCp
    (
        iter()() * iter()->Cp(p, T, patchI)
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpCp.ref() += iter()() * iter()->Cp(p, T, patchI);
    }

    return tmpCp;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::Cv() const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tmpCv
    (
        iter()() * iter()->Cv()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpCv.ref() += iter()() * iter()->Cv();
    }

    return tmpCv;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    auto iter = phaseModels_.cbegin();

    tmp<scalarField> tmpCv
    (
        iter()() * iter()->Cv(p, T, patchI)
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpCv.ref() += iter()() * iter()->Cv(p, T, patchI);
    }

    return tmpCv;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::gamma() const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tmpCp
    (
        iter()() * iter()->Cp()
    );

    tmp<volScalarField> tmpCv
    (
        iter()() * iter()->Cv()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpCp.ref() += iter()() * iter()->Cp();
        tmpCv.ref() += iter()() * iter()->Cv();
    }

    return (tmpCp/tmpCv);
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    return
    (
        gamma()().boundaryField()[patchI]
    );
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::Cpv() const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tmpCpv
    (
        iter()() * iter()->Cpv()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpCpv.ref() += iter()() * iter()->Cpv();
    }

    return tmpCpv;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    auto iter = phaseModels_.cbegin();

    tmp<scalarField> tmpCpv
    (
        iter()() * iter()->Cpv(p, T, patchI)
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpCpv.ref() += iter()() * iter()->Cpv(p, T, patchI);
    }

    return tmpCpv;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::CpByCpv() const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tmpCpByCpv
    (
        iter()() * iter()->CpByCpv()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpCpByCpv.ref() += iter()() * iter()->CpByCpv();
    }

    return tmpCpByCpv;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    auto iter = phaseModels_.cbegin();

    tmp<scalarField> tmpCpv
    (
        iter()().boundaryField()[patchI]
      * iter()->CpByCpv(p, T, patchI)
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpCpv.ref() +=
        (
            iter()().boundaryField()[patchI]
          * iter()->CpByCpv(p, T, patchI)
        );
    }

    return tmpCpv;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::W() const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::kappa() const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tmpkappa
    (
        iter()() * iter()->kappa()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpkappa.ref() += iter()() * iter()->kappa();
    }

    return tmpkappa;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::kappa(const label patchI) const
{
    auto iter = phaseModels_.cbegin();

    tmp<scalarField> tmpKappa
    (
        iter()().boundaryField()[patchI]
      * iter()->kappa(patchI)
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpKappa.ref() +=
        (
            iter()().boundaryField()[patchI]
          * iter()->kappa(patchI)
        );
    }

    return tmpKappa;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::alphahe() const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<volScalarField> talphaEff
    (
        phaseModelIter()()*phaseModelIter()->alphahe()
    );

    for (; phaseModelIter != phaseModels_.end(); ++phaseModelIter)
    {
        talphaEff.ref() += phaseModelIter()()*phaseModelIter()->alphahe();
    }

    return talphaEff;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::alphahe
(
    const label patchi
) const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<scalarField> talphaEff
    (
        phaseModelIter()().boundaryField()[patchi]
       *phaseModelIter()->alphahe(patchi)
    );

    for (; phaseModelIter != phaseModels_.end(); ++phaseModelIter)
    {
        talphaEff.ref() +=
            phaseModelIter()().boundaryField()[patchi]
           *phaseModelIter()->alphahe(patchi);
    }

    return talphaEff;
}


Foam::tmp<Foam::volScalarField>Foam::phaseSystem::kappaEff
(
    const volScalarField& kappat
) const
{
    tmp<volScalarField> kappaEff(kappa() + kappat);
    kappaEff.ref().rename("kappaEff");
    return kappaEff;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::kappaEff
(
    const scalarField& kappat,
    const label patchI
) const
{
    return kappa(patchI) + kappat;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::alphaEff
(
    const volScalarField& alphat
) const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tmpAlpha
    (
        iter()() * iter()->alpha()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpAlpha.ref() += iter()() * iter()->alpha();
    }

    tmpAlpha.ref() += alphat;

    return tmpAlpha;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::alphaEff
(
    const scalarField& alphat,
    const label patchI
) const
{
    auto iter = phaseModels_.cbegin();

    tmp<scalarField> tmpAlpha
    (
        iter()().boundaryField()[patchI]
      * iter()->alpha(patchI)
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpAlpha.ref() +=
        (
            iter()().boundaryField()[patchI]
          * iter()->alpha(patchI)
        );
    }

    tmpAlpha.ref() += alphat;

    return tmpAlpha;
}


const Foam::dimensionedScalar& Foam::phaseSystem::Prt() const
{
    return Prt_;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::mu() const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tmpMu
    (
        iter()() * iter()->mu()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpMu.ref() += iter()() * iter()->mu();
    }

    return tmpMu;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::mu(const label patchI) const
{
    auto iter = phaseModels_.cbegin();

    tmp<scalarField> tmpMu
    (
        iter()().boundaryField()[patchI]
      * iter()->mu(patchI)
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpMu.ref() +=
        (
            iter()().boundaryField()[patchI]
          * iter()->mu(patchI)
        );
    }

    return tmpMu;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::nu() const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tmpNu
    (
        iter()() * iter()->nu()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpNu.ref() += iter()() * iter()->nu();
    }

    return tmpNu;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::nu(const label patchI) const
{
    auto iter = phaseModels_.cbegin();

    tmp<scalarField> tmpNu
    (
        iter()().boundaryField()[patchI]
      * iter()->nu(patchI)
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpNu.ref() +=
        (
            iter()().boundaryField()[patchI]
          * iter()->nu(patchI)
        );
    }

    return tmpNu;
}


const Foam::surfaceScalarField& Foam::phaseSystem::phi() const
{
    return phi_;
}


Foam::surfaceScalarField& Foam::phaseSystem::phi()
{
    return phi_;
}


const Foam::surfaceScalarField& Foam::phaseSystem::rhoPhi() const
{
    return rhoPhi_;
}


Foam::surfaceScalarField& Foam::phaseSystem::rhoPhi()
{
    return rhoPhi_;
}


void Foam::phaseSystem::correct()
{
    forAllIters(phaseModels_, iter)
    {
        iter()->correct();
    }

    calcMu();
}


void Foam::phaseSystem::correctTurbulence()
{
    forAllIters(phaseModels_, iter)
    {
        iter()->correctTurbulence();
    }
}


const Foam::phaseSystem::phaseModelTable& Foam::phaseSystem::phases() const
{
    return phaseModels_;
}


Foam::phaseSystem::phaseModelTable& Foam::phaseSystem::phases()
{
    return phaseModels_;
}


const Foam::phaseSystem::phasePairTable&
Foam::phaseSystem::totalPhasePairs() const
{
    return totalPhasePairs_;
}


Foam::phaseSystem::phasePairTable& Foam::phaseSystem::totalPhasePairs()
{
    return totalPhasePairs_;
}


bool Foam::phaseSystem::incompressible() const
{
    forAllConstIters(phaseModels_, iter)
    {
        if (!iter()->thermo().incompressible())
        {
            return false;
        }
    }

    return true;
}


bool Foam::phaseSystem::incompressible(const word phaseName) const
{
    return phaseModels_[phaseName]->thermo().incompressible();
}


bool Foam::phaseSystem::isochoric() const
{
    forAllConstIters(phaseModels_, iter)
    {
        if (!iter()->thermo().isochoric())
        {
            return false;
        }
    }

    return true;
}


const Foam::fvMesh& Foam::phaseSystem::mesh() const
{
    return mesh_;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::phaseSystem::surfaceTensionForce() const
{
    auto tstf = tmp<surfaceScalarField>::New
    (
        IOobject
        (
            "surfaceTensionForce",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar({1, -2, -2, 0, 0, 0}, Zero)
    );

    auto& stf = tstf.ref();
    stf.setOriented();

    if (surfaceTensionModels_.size())
    {
        forAllConstIters(phaseModels_, iter1)
        {
            const volScalarField& alpha1 = iter1()();

            auto iter2 = iter1;

            for (++iter2; iter2 != phaseModels_.cend(); ++iter2)
            {
                const volScalarField& alpha2 = iter2()();

                stf +=
                    fvc::interpolate
                    (
                        surfaceTensionCoeff
                        (
                            phasePairKey(iter1()->name(), iter2()->name())
                        )
                    )
                * fvc::interpolate(K(alpha1, alpha2))*
                    (
                        fvc::interpolate(alpha2)*fvc::snGrad(alpha1)
                      - fvc::interpolate(alpha1)*fvc::snGrad(alpha2)
                    );
            }
        }
    }

    return tstf;
}

Foam::tmp<Foam::volVectorField>
Foam::phaseSystem::divCapillaryStress() const
{
    auto tcst = tmp<volVectorField>::New
    (
        IOobject
        (
            "capillaryStressTensor",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector({1, -2, -2, 0, 0, 0}, Zero)
    );

    auto& cst = tcst.ref();

    /*
        Implementation of the capillary stress tensor
        with the Continuous Surface Stress Method
        iaw Brackbill et al., 1992, Lafaurie et al. 1994
        T = - sigma * (I - n x n) * mag(grad(alpha))
    */

    if (surfaceTensionModels_.size())
    {
        forAllConstIters(phaseModels_, iter1)
        {
            const volScalarField& alpha1 = iter1()();

            auto iter2 = iter1;

            for (++iter2; iter2 != phaseModels_.cend(); ++iter2)
            {
                const volScalarField& alpha2 = iter2()();

                volVectorField gradAlpha
                (
                    alpha2*fvc::grad(alpha1)
                    - alpha1*fvc::grad(alpha2)
                );

                cst +=
                    fvc::div
                    (
                        mag(gradAlpha)
                        *
                        (
                            surfaceTensionCoeff
                            (
                                phasePairKey(iter1()->name(), iter2()->name())
                            )
                            *
                            (
                                tensor::I
                                -
                                nHat(alpha1,alpha2) * nHat(alpha1,alpha2)
                            )
                        )
                    )
                    ;
            }
        }
    }

    return tcst;
    
}

Foam::tmp<Foam::volVectorField> Foam::phaseSystem::U() const
{
    auto tstf = tmp<volVectorField>::New
    (
        IOobject
        (
            "U",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimVelocity, Zero)
    );

    auto& stf = tstf.ref();

    forAllConstIters(phaseModels_, iter)
    {
        stf += iter()() * iter()->U();
    }

    return tstf;
}


Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::surfaceTensionCoeff(const phasePairKey& key) const
{
    return surfaceTensionModels_[key]->sigma();
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::coeffs
(
    const word& key
) const
{
    return 1.0/(phaseModels_[key]->thermo().rho());
}


void Foam::phaseSystem::addInterfacePorosity(fvVectorMatrix& UEqn)
{
    const scalarField& Vc = mesh_.V();
    scalarField& Udiag = UEqn.diag();

    forAllConstIters(phaseModels_, iteri)
    {
        const phaseModel& phasei = iteri()();

        auto iterk = iteri;

        for (++iterk; iterk != phaseModels_.cend(); ++iterk)
        {
            if (iteri()().name() != iterk()().name())
            {
                const phaseModel& phasek = iterk()();

                // Phase i and k
                const phasePairKey keyik
                (
                    phasei.name(),
                    phasek.name(),
                    false
                );

                if (interfacePorousModelTable_.found(keyik))
                {
                    autoPtr<porousModel>& interfacePtr =
                        interfacePorousModelTable_[keyik];

                    Udiag += Vc*interfacePtr->S();
                }
            }
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::K
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    tmp<surfaceVectorField> tnHatfv = nHatfv(alpha1, alpha2);

    // Simple expression for curvature
    return -fvc::div(tnHatfv.ref() & mesh_.Sf());
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::nearInterface
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    return
    (
        pos(alpha1 - 0.1)*pos(0.9 - alpha1)
       *pos(alpha2 - 0.1)*pos(0.9 - alpha2)
    );
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::nearInterface() const
{
    auto tnearInt = tmp<volScalarField>::New
    (
        IOobject
        (
            "nearInterface",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );

    auto& nearInt = tnearInt.ref();

    forAllConstIters(phaseModels_, iter1)
    {
        const volScalarField& alpha1 = iter1()();

        auto iter2 = iter1;

        for (++iter2; iter2 != phaseModels_.cend(); ++iter2)
        {
            const volScalarField& alpha2 = iter2()();

            nearInt +=
            (
                pos(alpha1 - 0.1)*pos(0.9 - alpha1)
               *pos(alpha2 - 0.1)*pos(0.9 - alpha2)
            );
        }
    }

    return tnearInt;
}


Foam::tmp<Foam::surfaceVectorField> Foam::phaseSystem::nHatfv
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{

    surfaceVectorField gradAlphaf
    (
        fvc::interpolate(alpha2)*fvc::interpolate(fvc::grad(alpha1))
      - fvc::interpolate(alpha1)*fvc::interpolate(fvc::grad(alpha2))
    );

    const dimensionedScalar deltaN
    (
        "deltaN",
        1e-8/cbrt(average(mesh_.V()))
    );

    // Face unit interface normal
    return gradAlphaf/(mag(gradAlphaf) + deltaN);
}

Foam::tmp<Foam::volVectorField> Foam::phaseSystem::nHat
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{

    volVectorField gradAlphaf
    (
        alpha2*fvc::grad(alpha1)
      - alpha1*fvc::grad(alpha2)
    );

    const dimensionedScalar deltaN
    (
        "deltaN",
        1e-8/cbrt(average(mesh_.V()))
    );

    // Face unit interface normal
    return gradAlphaf/(mag(gradAlphaf) + deltaN);
}

Foam::tmp<Foam::surfaceScalarField> Foam::phaseSystem::nHatf
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    // Face unit interface normal flux
    return nHatfv(alpha1, alpha2) & mesh_.Sf();
}


bool Foam::phaseSystem::read()
{
    if (regIOobject::read())
    {
        return true;
    }

    return false;
}


// ************************************************************************* //
