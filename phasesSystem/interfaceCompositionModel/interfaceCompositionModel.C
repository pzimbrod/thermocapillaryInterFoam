/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "interfaceCompositionModel.H"
#include "phaseModel.H"
#include "phasePair.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interfaceCompositionModel, 0);
    defineRunTimeSelectionTable(interfaceCompositionModel, dictionary);
}


const Foam::Enum<Foam::interfaceCompositionModel::modelVariable>
Foam::interfaceCompositionModel::modelVariableNames_
{
    { modelVariable::T, "temperature" },
    { modelVariable::P, "pressure" },
    { modelVariable::Y, "massFraction" },
    { modelVariable::alpha, "alphaVolumeFraction" },
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceCompositionModel::interfaceCompositionModel
(
    const dictionary& dict,
    const phasePair& pair
)
:
    modelVariable_
    (
        modelVariableNames_.getOrDefault
        (
            "variable",
            dict,
            modelVariable::T
        )
    ),
    includeVolChange_(dict.getOrDefault("includeVolChange", true)),
    pair_(pair),
    speciesName_(dict.getOrDefault<word>("species", "none")),
    mesh_(pair_.from().mesh())
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::interfaceCompositionModel>
Foam::interfaceCompositionModel::New
(
    const dictionary& dict,
    const phasePair& pair
)
{
    const word modelType
    (
        dict.get<word>("type")
      + "<"
      + pair.phase1().thermo().type()
      + ","
      + pair.phase2().thermo().type()
      + ">"
    );

    Info<< "Selecting interfaceCompositionModel for "
        << pair << ": " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "interfaceCompositionModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(dict, pair);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word Foam::interfaceCompositionModel::transferSpecie() const
{
    return speciesName_;
}


const Foam::phasePair& Foam::interfaceCompositionModel::pair() const
{
    return pair_;
}


const Foam::word& Foam::interfaceCompositionModel::variable() const
{
    return modelVariableNames_[modelVariable_];
}


bool Foam::interfaceCompositionModel::includeDivU() const noexcept
{
    return true;
}


bool Foam::interfaceCompositionModel::includeVolChange()
{
    return includeVolChange_;
}


// ************************************************************************* //
