/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::kOmegaSSTPDABase

Description
    
    This is the openFoam kOmegaSST model corrected with separation factor and secondary flow prediction,
    Progressive data augmented (PDA) 2023, A. Amarloo, M. Rincon

    More details and test cases in following publications:
    A. Amarloo, M. J. Rincon, M. Reclari, and M. Abkar, 
    "Progressive augmentation of RANS models for separated flow prediction
    by CFD-driven surrogate multi-objective optimisation." (2023).
    and
    M. J. Rincon, A. Amarloo,  M. Reclari, X.I.A. Yang, and M. Abkar, 
    "Progressive augmentation of Reynolds stress tensor models for secondary flow prediction 
    by computational fluid dynamics driven surrogate optimisation" (2023).
    
    The model coefficients  are
    \verbatim
        kOmegaSSTPDABaseCoeffs
        {
            //original coefficients of KOSST
            alphaK1         0.85;
            alphaK2         1.0;
            alphaOmega1     0.5;
            alphaOmega2     0.856;
            beta1           0.075;
            beta2           0.0828;
            betaStar        0.09;
            gamma1          5/9;
            gamma2          0.44;
            a1              0.31;
            b1              1.0;
            c1              10.0;
            F3              no;


            //Separation coefficients
            separationMode      4;              //optional - default:4 - off:0 | ModelI:1 | ModelII:2 | ModelIII:3 | ModelIV:4
            separationLambda1   1;              //optional - default taken from separationMode 4
            separationLambda2   1;              //optional - default taken from separationMode 4
            C0                  -1;             //optional - default taken from separationMode 4
            C1                  0;              //optional - default taken from separationMode 4
            C2                  0;              //optional - default taken from separationMode 4


            //Secodnary Flow coefficients
            secondaryMode       2;              //optional - default:2 - off:0 | ModelI:1 | ModelII:2
            A0                  -1;             //optional - default taken from secondaryMode 2
            A1                  0;              //optional - default taken from secondaryMode 2
            A2                  0;              //optional - default taken from secondaryMode 2


            // Optional decay control
            decayControl    yes;
            kInf            \<far-field k value\>;
            omegaInf        \<far-field omega value\>;
        }
    \endverbatim

SourceFiles
    kOmegaSSTPDABase.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTPDABase_H
#define kOmegaSSTPDABase_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class kOmegaSSTPDABase Declaration
\*---------------------------------------------------------------------------*/

template<class BasicEddyViscosityModel>
class kOmegaSSTPDABase
:
    public BasicEddyViscosityModel
{
    // Private Member Functions

        //- No copy construct
        kOmegaSSTPDABase(const kOmegaSSTPDABase&) = delete;

        //- No copy assignment
        void operator=(const kOmegaSSTPDABase&) = delete;


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar alphaK1_;
            dimensionedScalar alphaK2_;

            dimensionedScalar alphaOmega1_;
            dimensionedScalar alphaOmega2_;

            dimensionedScalar gamma1_;
            dimensionedScalar gamma2_;

            dimensionedScalar beta1_;
            dimensionedScalar beta2_;

            dimensionedScalar betaStar_;

            dimensionedScalar a1_;
            dimensionedScalar b1_;
            dimensionedScalar c1_;


            dimensionedScalar separationMode_;
            dimensionedScalar separationLambda1_;
            dimensionedScalar separationLambda2_;
            dimensionedScalar C0_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;


            dimensionedScalar secondaryMode_;
            dimensionedScalar secondary_relaxation_;
            dimensionedScalar A0_;
            dimensionedScalar A1_;
            dimensionedScalar A2_;


            //- Flag to include the F3 term
            Switch F3_;


        // Fields

            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField& y_;

            volScalarField k_;
            volScalarField omega_;
            volScalarField separationFactor_;
            volSymmTensorField bijDelta_;
            volSymmTensorField Rij_;
        // Decay control

            //- Flag to include the decay control
            Switch decayControl_;
            dimensionedScalar kInf_;
            dimensionedScalar omegaInf_;


    // Protected Member Functions

        void setDecayControl(const dictionary& dict);

        virtual tmp<volScalarField> F1(const volScalarField& CDkOmega) const;
        virtual tmp<volScalarField> F2() const;
        virtual tmp<volScalarField> F3() const;
        virtual tmp<volScalarField> F23() const;

        tmp<volScalarField> blend
        (
            const volScalarField& F1,
            const dimensionedScalar& psi1,
            const dimensionedScalar& psi2
        ) const
        {
            return F1*(psi1 - psi2) + psi2;
        }

        tmp<volScalarField::Internal> blend
        (
            const volScalarField::Internal& F1,
            const dimensionedScalar& psi1,
            const dimensionedScalar& psi2
        ) const
        {
            return F1*(psi1 - psi2) + psi2;
        }

        tmp<volScalarField> alphaK(const volScalarField& F1) const
        {
            return blend(F1, alphaK1_, alphaK2_);
        }

        tmp<volScalarField> alphaOmega(const volScalarField& F1) const
        {
            return blend(F1, alphaOmega1_, alphaOmega2_);
        }

        tmp<volScalarField::Internal> beta
        (
            const volScalarField::Internal& F1
        ) const
        {
            return tmp<volScalarField::Internal>::New
            (
                this->type() + ":beta",
                blend(F1, beta1_, beta2_)
            );
        }

        tmp<volScalarField::Internal> gamma
        (
            const volScalarField::Internal& F1
        ) const
        {
            return tmp<volScalarField::Internal>::New
            (
                this->type() + ":gamma",
                blend(F1, gamma1_, gamma2_)
            );
        }

        virtual void correctNut(const volScalarField& S2);

        virtual void correctNut();

        //- Return k production rate
        virtual tmp<volScalarField::Internal> Pk
        (
            const volScalarField::Internal& G
        ) const;

        //- Return epsilon/k which for standard RAS is betaStar*omega
        virtual tmp<volScalarField::Internal> epsilonByk
        (
            const volScalarField& F1,
            const volTensorField& gradU
        ) const;

        //- Return G/nu
        virtual tmp<volScalarField::Internal> GbyNu
        (
            const volScalarField::Internal& GbyNu0,
            const volScalarField::Internal& F2,
            const volScalarField::Internal& S2
        ) const;

        virtual tmp<fvScalarMatrix> kSource() const;

        virtual tmp<fvScalarMatrix> omegaSource() const;

        virtual tmp<fvScalarMatrix> Qsas
        (
            const volScalarField::Internal& S2,
            const volScalarField::Internal& gamma,
            const volScalarField::Internal& beta
        ) const;


public:

    typedef typename BasicEddyViscosityModel::alphaField alphaField;
    typedef typename BasicEddyViscosityModel::rhoField rhoField;
    typedef typename BasicEddyViscosityModel::transportModel transportModel;


    // Constructors

        //- Construct from components
        kOmegaSSTPDABase
        (
            const word& type,
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName
        );


    //- Destructor
    virtual ~kOmegaSSTPDABase() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", alphaK(F1)*this->nut_ + this->nu())
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DomegaEff",
                    alphaOmega(F1)*this->nut_ + this->nu()
                )
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kOmegaSSTPDABase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
