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

#include "localDensityAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(localDensityAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            localDensityAbsorptionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::volScalarField&
Foam::radiation::localDensityAbsorptionEmission::alpha(word alphaName) const
{
    if (!mesh_.foundObject<volScalarField>(alphaName))
    {
        FatalErrorInFunction
            << "Unable to retrieve density field " << alphaName << " from "
            << "database.  Available objects:" << mesh_.sortedNames()
            << exit(FatalError);
    }

    return mesh_.lookupObject<volScalarField>(alphaName);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::localDensityAbsorptionEmission::localDensityAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    alphaNames_(coeffsDict_.lookup("alphaNames")),
    aCoeff_(coeffsDict_.lookup("aCoeff")),
    eCoeff_(coeffsDict_.lookup("eCoeff")),
    ECoeff_(coeffsDict_.lookup("ECoeff"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::localDensityAbsorptionEmission::aCont(const label bandI) const
{
    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "a",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(inv(dimLength), Zero)
        )
    );

    volScalarField& a = ta.ref();

    forAll(alphaNames_, i)
    {
        dimensionedScalar aPhase("a", dimless/dimLength, aCoeff_[i]);
        a += max(alpha(alphaNames_[i]), scalar(0))*aPhase;
    }

    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::localDensityAbsorptionEmission::eCont(const label bandI) const
{
    tmp<volScalarField> te
    (
        new volScalarField
        (
            IOobject
            (
                "e",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(inv(dimLength), Zero)
        )
    );

    volScalarField& e = te.ref();

    forAll(alphaNames_, i)
    {
        dimensionedScalar ePhase("e", dimless/dimLength, eCoeff_[i]);
        e += max(alpha(alphaNames_[i]), scalar(0))*ePhase;
    }

    return te;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::localDensityAbsorptionEmission::ECont(const label bandI) const
{
    tmp<volScalarField> tE
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
        )
    );

    scalarField& E = tE.ref().primitiveFieldRef();

    forAll(alphaNames_, i)
    {
        dimensionedScalar EPhase
        (
            "E",
            dimMass/dimLength/pow3(dimTime),
            ECoeff_[i]
        );

        E += max(alpha(alphaNames_[i]), scalar(0))*EPhase;
    }

    return tE;
}


// ************************************************************************* //
