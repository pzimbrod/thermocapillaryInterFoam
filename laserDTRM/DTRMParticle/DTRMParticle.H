/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 OpenCFD Ltd
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

Class
    Foam::DTRMParticle

Description
    Discrete Transfer Radiation Model (DTRM) particle

SourceFiles
    DTRMParticle.H

\*---------------------------------------------------------------------------*/

#ifndef DTRMParticle_H
#define DTRMParticle_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "interpolationCell.H"
#include "volFieldsFwd.H"
#include "reflectionModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class DTRMParticle;
using namespace Foam::radiation;

// Forward declaration of friend functions
Ostream& operator<<(Ostream&, const DTRMParticle&);

/*---------------------------------------------------------------------------*\
                        Class DTRMParticle Declaration
\*---------------------------------------------------------------------------*/

class DTRMParticle
:
    public particle
{
    // Private data

        //- Initial position
        point p0_;

        //- Target position
        point p1_;

        //- Initial radiation intensity [W/m2]
        scalar I0_;

        //- Radiation intensity [W/m2]
        scalar I_;

        //- Area of radiation
        scalar dA_;

        //- Trasnmissive index
        label transmissiveId_;


public:

    friend class Cloud<DTRMParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        // Interpolators for continuous phase fields

            const interpolationCell<scalar>& aInterp_;
            const interpolationCell<scalar>& eInterp_;
            const interpolationCell<scalar>& EInterp_;
            const interpolationCell<scalar>& TInterp_;
            const interpolationCellPoint<vector>& nHatInterp_;

            //- Reflected cells
            const labelField& relfectedCells_;

            //- Ptr to reflectiom model
            UPtrList<reflectionModel> reflection_;

            //- Heat source term
            volScalarField& Q_;


    public:

        // Constructors

            inline trackingData
            (
                Cloud<DTRMParticle>& spc,
                const interpolationCell<scalar>& aInterp,
                const interpolationCell<scalar>& eInterp,
                const interpolationCell<scalar>& EInterp,
                const interpolationCell<scalar>& TInterp,
                const interpolationCellPoint<vector>& nHatInterp,
                const labelField&,
                const UPtrList<reflectionModel>&,
                volScalarField& Q
            );


        // Member functions

            inline const interpolationCell<scalar>& aInterp() const;
            inline const interpolationCell<scalar>& eInterp() const;
            inline const interpolationCell<scalar>& EInterp() const;
            inline const interpolationCell<scalar>& TInterp() const;
            inline const interpolationCellPoint<vector>& nHatInterp() const;
            inline const labelField& relfectedCells() const;
            inline const UPtrList<reflectionModel>& reflection() const;

            inline scalar& Q(label celli);
    };


    // Static Data Members

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;


        //- String representation of properties
        AddToPropertyList
        (
            particle,
            " p0"
          + " p1"
          + " I0"
          + " I"
          + " dA"
          + " transmissiveId";
        );


    // Constructors

        //- Construct from components, with searching for tetFace and
        //  tetPt unless disabled by doCellFacePt = false.
        DTRMParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const vector& targetPosition,
            const scalar I,
            const label cellI,
            const scalar dA,
            const label transmissiveId
        );

        //- Construct from components
        DTRMParticle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const vector& position,
            const vector& targetPosition,
            const scalar I,
            const scalar dA,
            const label transmissiveId
        );

        //- Construct from Istream
        DTRMParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as copy
        DTRMParticle(const DTRMParticle& p);


        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<DTRMParticle> operator()(Istream& is) const
            {
                return autoPtr<DTRMParticle>
                (
                    new DTRMParticle(mesh_, is, true)
                );
            }
        };


    // Access

        //- Return const access to the initial position
        inline const point& p0() const;

        //- Return const access to the target position
        inline const point& p1() const;

        //- Return const access to the initial intensity
        inline scalar I0() const;

        //- Return const access to the current intensity
        inline scalar I() const;

        //- Return const access dA
        inline scalar dA() const;


    // Edit

        //- Return access to the target position
        inline point& p1();

        //- Return access to the initial intensity
        inline scalar& I0();

        //- Return access to the current intensity
        inline scalar& I();

        //- Return access to dA
        inline scalar& dA();

        //- Return access to reflectedId
        inline label& reflectedId();


   // Tracking

        //- Move
        bool move(Cloud<DTRMParticle>& , trackingData&, const scalar);


    // Member Operators

        //- Overridable function to handle the particle hitting a processorPatch
        void hitProcessorPatch
        (
            Cloud<DTRMParticle>&,
            trackingData& td
        );

        //- Overridable function to handle the particle hitting a wallPatch
        void hitWallPatch
        (
            Cloud<DTRMParticle>&,
            trackingData& td
        );

        bool hitPatch
        (
            Cloud<DTRMParticle>&,
            trackingData& td
        );


        // I-O

            //- Write individual parcel properties to stream
            void writeProperties
            (
                Ostream& os,
                const wordRes& filters,
                const word& delim,
                const bool namesOnly = false
            ) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const DTRMParticle& p);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DTRMParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
