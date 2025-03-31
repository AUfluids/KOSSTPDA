/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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



\*---------------------------------------------------------------------------*/

#include "kOmegaSSTPDABase.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTPDABase<BasicEddyViscosityModel>::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTPDABase<BasicEddyViscosityModel>::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTPDABase<BasicEddyViscosityModel>::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omega_*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTPDABase<BasicEddyViscosityModel>::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23.ref() *= F3();
    }

    return f23;
}


template<class BasicEddyViscosityModel>
void kOmegaSSTPDABase<BasicEddyViscosityModel>::correctNut
(
    const volScalarField& S2
)
{
    // Correct the turbulence viscosity
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F23()*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
}


template<class BasicEddyViscosityModel>
void kOmegaSSTPDABase<BasicEddyViscosityModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}

template<class BasicEddyViscosityModel>
tmp<fvVectorMatrix> kOmegaSSTPDABase<BasicEddyViscosityModel>::divDevReff
(
    volVectorField& U
) const
{
    //Info << "In: nonlinearRST::divDevReff()" << endl;
    return
    (
      - fvc::div((this->alpha_*this->rho_*this->nuEff())*dev2(T(fvc::grad(U))))
      - fvm::laplacian(this->alpha_*this->rho_*this->nuEff(), U)    // linear part
      + fvc::div(dev(2.*this->k_*this->bijDelta_))  // non-linear part
    );
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal> kOmegaSSTPDABase<BasicEddyViscosityModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    return min(G, (c1_*betaStar_)*this->k_()*this->omega_());
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal>
kOmegaSSTPDABase<BasicEddyViscosityModel>::epsilonByk
(
    const volScalarField& F1,
    const volTensorField& gradU
) const
{
    return betaStar_*omega_();
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal> kOmegaSSTPDABase<BasicEddyViscosityModel>::GbyNu
(
    const volScalarField::Internal& GbyNu0,
    const volScalarField::Internal& F2,
    const volScalarField::Internal& S2
) const
{
    return min
    (
        GbyNu0,
        (c1_/a1_)*betaStar_*omega_()
       *max(a1_*omega_(), b1_*F2*sqrt(S2))
    );
}


template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> kOmegaSSTPDABase<BasicEddyViscosityModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
        )
    );
}


template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> kOmegaSSTPDABase<BasicEddyViscosityModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> kOmegaSSTPDABase<BasicEddyViscosityModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
kOmegaSSTPDABase<BasicEddyViscosityModel>::kOmegaSSTPDABase
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicEddyViscosityModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    // Added by Mario

    I1_
    (
        IOobject
        (
            IOobject::groupName("I1", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,                     // Use 'this->mesh_' to initialize the field
        dimensionedScalar("zero", dimless, 0)
    ),
    I2_
    (
        IOobject
        (
            IOobject::groupName("I2", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimless, 0)
    ),
    I3_
    (
        IOobject
        (
            IOobject::groupName("I3", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimless, 0)
    ),
    I4_
    (
        IOobject
        (
            IOobject::groupName("I4", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimless, 0)
    ),
    I5_
    (
        IOobject
        (
            IOobject::groupName("I5", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimless, 0)
    ),
    Sij_
    (
        IOobject
        (
            IOobject::groupName("Sij", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    Omegaij_
    (
        IOobject
        (
            IOobject::groupName("Omegaij", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedTensor("zero", dimless, tensor::zero)
    ),
    Tij1_
    (
        IOobject
        (
            IOobject::groupName("Tij1", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    Tij2_
    (
        IOobject
        (
            IOobject::groupName("Tij2", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    Tij3_
    (
        IOobject
        (
            IOobject::groupName("Tij3", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    Tij4_
    (
        IOobject
        (
            IOobject::groupName("Tij4", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    Tij5_
    (
        IOobject
        (
            IOobject::groupName("Tij5", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),

    // Normal model

    alphaK1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    separationMode_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "separationMode",
            this->coeffDict_,
            3.0
        )
    ),
    separationLambda1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "separationLambda1",
            this->coeffDict_,
            0.0
        )
    ),
    separationLambda2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "separationLambda2",
            this->coeffDict_,
            0.0
        )
    ),
    C0_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C0",
            this->coeffDict_,
            0.0
        )
    ),
    C1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C1",
            this->coeffDict_,
            0.0
        )
    ),
    C2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C2",
            this->coeffDict_,
            0.0
        )
    ),
    secondaryMode_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "secondaryMode",
            this->coeffDict_,
            2.0
        )
    ),
    bijDelta_
    (
        IOobject
        (
            IOobject::groupName("bijDelta", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        0*dev(symm(fvc::grad(U)))/(omega_ + this->omegaMin_)
    ),
    Rij_
    (
        IOobject
        (
            IOobject::groupName("Rij", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        0.0*(((2.0/3.0)*I)*k_ - this->nut_*twoSymm(fvc::grad(this->U_)))
    ),
    secondary_relaxation_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "secondary_relaxation",
            this->coeffDict_,
            0.5
        )
    ),
    A0_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A0",
            this->coeffDict_,
            0.0
        )
    ),
    A1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A1",
            this->coeffDict_,
            0.0
        )
    ),
    A2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A2",
            this->coeffDict_,
            0.0
        )
    ),
    F3_
    (
        Switch::getOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),

    y_(wallDist::New(this->mesh_).y()),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    decayControl_
    (
        Switch::getOrAddToDict
        (
            "decayControl",
            this->coeffDict_,
            false
        )
    ),
    kInf_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kInf",
            this->coeffDict_,
            k_.dimensions(),
            0
        )
    ),
    omegaInf_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "omegaInf",
            this->coeffDict_,
            omega_.dimensions(),
            0
        )
    ),
    separationFactor_
    (
        IOobject
        (
            IOobject::groupName("separationFactor", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        0*k_/k_
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    setDecayControl(this->coeffDict_);

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
void kOmegaSSTPDABase<BasicEddyViscosityModel>::setDecayControl
(
    const dictionary& dict
)
{
    decayControl_.readIfPresent("decayControl", dict);

    if (decayControl_)
    {
        kInf_.read(dict);
        omegaInf_.read(dict);

        Info<< "    Employing decay control with kInf:" << kInf_
            << " and omegaInf:" << omegaInf_ << endl;
    }
    else
    {
        kInf_.value() = 0;
        omegaInf_.value() = 0;
    }
}


template<class BasicEddyViscosityModel>
bool kOmegaSSTPDABase<BasicEddyViscosityModel>::read()
{
    if (BasicEddyViscosityModel::read())
    {
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());
        separationMode_.readIfPresent(this->coeffDict());
        separationLambda1_.readIfPresent(this->coeffDict());
        separationLambda2_.readIfPresent(this->coeffDict());
        C0_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        secondaryMode_.readIfPresent(this->coeffDict());
        A0_.readIfPresent(this->coeffDict());
        A1_.readIfPresent(this->coeffDict());
        A2_.readIfPresent(this->coeffDict());
        setDecayControl(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicEddyViscosityModel>
void kOmegaSSTPDABase<BasicEddyViscosityModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    BasicEddyViscosityModel::correct();

    volScalarField::Internal divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));


    // New calculations by Mario

    // Write the fields to the results directory
    // I1_.write();
    // Tij1_.write();
    // Sij_.write();
    // Omegaij_.write();

    // New calculations

    volTensorField gradU(fvc::grad(U));
    volSymmTensorField Sij(symm(gradU));
    volTensorField Omegaij(-0.5*(gradU - gradU.T()));
    volScalarField S(sqrt(2*magSqr(symm(fvc::grad(U)))));

    volScalarField tScale(1./max( S/a1_ + this->omegaMin_,omega_ + this->omegaMin_));
    volScalarField tScale2(tScale*tScale);
    volScalarField tScale3(tScale*tScale2);
    volScalarField tScale4(tScale*tScale3);
    volScalarField tScale5(tScale*tScale4);

    volScalarField I1_(tScale2 * tr(Sij & Sij));
    volScalarField I2_(tScale2 * tr(Omegaij & Omegaij));
    volScalarField I3_(tScale3 * tr((Sij & Sij) & Sij));
    volScalarField I4_(tScale3 * tr((Omegaij & Omegaij) & Sij));
    volScalarField I5_(tScale4 * tr((Omegaij & Omegaij) & (Sij & Sij)));

    volSymmTensorField Tij1_(tScale * Sij);
    volSymmTensorField Tij2_(tScale2 * symm((Sij & Omegaij) - (Omegaij & Sij)));
    volSymmTensorField Tij3_(tScale3 * symm((Sij & Sij) & Omegaij - Omegaij & (Sij & Sij)));
    volSymmTensorField Tij4_(tScale4 * symm((Omegaij & Omegaij) & Sij - Sij & (Omegaij & Omegaij)));
    volSymmTensorField Tij5_(tScale5 * symm((Omegaij & Omegaij) & (Sij & Sij) - (Sij & Sij) & (Omegaij & Omegaij)));

    forAll(Sij, CellI)
    {
        // Compute invariants
        I1_[CellI] = tScale2[CellI] * tr(Sij[CellI] & Sij[CellI]);
        I2_[CellI] = tScale2[CellI] * tr(Omegaij[CellI] & Omegaij[CellI]);
        I3_[CellI] = tScale3[CellI] * tr((Sij[CellI] & Sij[CellI]) & Sij[CellI]);
        I4_[CellI] = tScale3[CellI] * tr((Omegaij[CellI] & Omegaij[CellI]) & Sij[CellI]);
        I5_[CellI] = tScale4[CellI] * tr((Omegaij[CellI] & Omegaij[CellI]) & (Sij[CellI] & Sij[CellI]));

        // Compute tensors
        Tij1_[CellI] = tScale[CellI] * Sij[CellI];
        Tij2_[CellI] = tScale2[CellI] * symm((Sij[CellI] & Omegaij[CellI]) - (Omegaij[CellI] & Sij[CellI]));
        Tij3_[CellI] = tScale3[CellI] * symm((Sij[CellI] & Sij[CellI]) & Omegaij[CellI] - Omegaij[CellI] & (Sij[CellI] & Sij[CellI]));
        Tij4_[CellI] = tScale4[CellI] * symm((Omegaij[CellI] & Omegaij[CellI]) & Sij[CellI] - Sij[CellI] & (Omegaij[CellI] & Omegaij[CellI]));
        Tij5_[CellI] = tScale5[CellI] * symm((Omegaij[CellI] & Omegaij[CellI]) & (Sij[CellI] & Sij[CellI]) - (Sij[CellI] & Sij[CellI]) & (Omegaij[CellI] & Omegaij[CellI]));
    }


    dimensionedScalar nutMin("nutMin", dimensionSet(0, 2, -1, 0, 0, 0 ,0), 1e-9);

    volScalarField alpha1_separation(I1_ * 0.0);

    // Separation Factor

    // if (separationMode_.value() == 1.0)
    // {
    //     if (C0_.value() == 0.0 && C1_.value() == 0.0 && C2_.value() == 0.0)
    //     {
    //         C0_ = -1.08213;
    //         C1_ = -0.377233;
    //         C2_ = 0.062756;
    //     }
    //     separationLambda1_ = 1.0;
    //     separationLambda2_ = 1.0;
    // }

    // if (separationMode_.value() == 2.0)
    // {
    //     if (C0_.value() == 0.0 && C1_.value() == 0.0 && C2_.value() == 0.0)
    //     {
    //         C0_ = -1.85884;
    //         C1_ = 0.0626469;
    //         C2_ = 0.0157284;
    //     }
    //     separationLambda1_ = 1.0;
    //     separationLambda2_ = 1.0;
    // }

    // if (separationMode_.value() == 3.0)
    // {
    //     if (C0_.value() == 0.0 && C1_.value() == 0.0 && C2_.value() == 0.0 && separationLambda1_.value() == 0.0 && separationLambda2_.value() == 0.0)
    //     {
    //         C0_ = -0.252012;
    //         C1_ = -0.441849;
    //         C2_ = -0.0254661;
    //         separationLambda1_ = 16.4685;
    //         separationLambda2_ = 3.99033;
    //     }
    // }

    // if (separationMode_.value() == 4.0)
    // {
    //     if (C0_.value() == 0.0 && C1_.value() == 0.0 && C2_.value() == 0.0 && separationLambda1_.value() == 0.0 && separationLambda2_.value() == 0.0)
    //     {
    //         C0_ = -0.872209;
    //         C1_ = 0.0131861;
    //         C2_ = -0.0766894;
    //         separationLambda1_ = 20;
    //         separationLambda2_ = 7.2513;
    //     }
    // }



    // if (separationMode_.value() == 1.0 || separationMode_.value() == 3.0) // average between PH2800 and CBFS13700
    // {
        alpha1_separation = C0_
                 + C1_*(I1_- (2.86797085e-02 + 3.22109598e-02) / 2.0) / ((1.96630250e-02 + 1.69618216e-02) / 2.0)
                 + C2_*(I2_ - (-1.21140076e-02 -2.59699426e-02) / 2.0) / ((1.83587958e-02 + 1.65637539e-02) / 2.0);
    // }


    // if (separationMode_.value() == 2.0 || separationMode_.value() == 4.0)  // I1I2 for PH2800 + CBFS13700
    // {

    //     std::vector<float> Mean_Funcs_SEP ={0.0 , 3.04453342e-02, -1.90419751e-02,  1.26720478e-03, -7.63523411e-04,
    //                                             7.16295213e-04,  5.66460407e-05, -3.35790350e-05,  3.07466925e-05,
    //                                            -3.12439589e-05,  2.60333900e-06, -1.52808501e-06,  1.38771303e-06,
    //                                            -1.38617717e-06,  1.49221337e-06,  1.21315469e-07, -7.07001591e-08,
    //                                             6.39716068e-08, -6.34625036e-08,  6.73430685e-08, -7.66798242e-08,
    //                                             5.69918377e-09, -3.30261711e-09,  2.98224267e-09, -2.94830635e-09,
    //                                             3.10670530e-09,  3.10670530e-09, -3.48853434e-09};

    //     // scalarList Std_Funcs(21);
    //     std::vector<float> Std_Funcs_SEP ={1.0 , 1.84468535e-02, 1.88068710e-02, 9.98764762e-04, 8.97075822e-04,
    //                                             9.89512270e-04, 4.99040063e-05, 4.30661244e-05, 4.64902807e-05,
    //                                             5.69466498e-05, 2.44483400e-06, 2.06574320e-06, 2.20975847e-06,
    //                                             2.64454589e-06, 3.66260083e-06, 1.18875010e-07, 9.91138283e-08,
    //                                             1.05445583e-07, 1.24745863e-07, 1.68475795e-07, 2.61569539e-07,
    //                                             5.75863939e-09, 4.75739289e-09, 5.04256885e-09, 5.92425621e-09,
    //                                             7.89044848e-09, 7.89044848e-09, 1.19776568e-08};

    //     // scalarList PC1_Coef(21);
    //     std::vector<float> PC1_Coef_SEP =  {0.0 , 1.25216111e-01, -1.74322502e-01,  1.33197282e-01,
    //                                             -1.96654511e-01,  2.02496020e-01,  1.34887383e-01,
    //                                             -2.01363658e-01,  2.12796015e-01, -2.10940176e-01,
    //                                              1.34807484e-01, -2.02309964e-01,  2.15034456e-01,
    //                                             -2.17782965e-01,  2.04606963e-01,  1.34144487e-01,
    //                                             -2.02056548e-01,  2.15135671e-01, -2.19542292e-01,
    //                                              2.10492582e-01, -1.88316897e-01,  1.33293656e-01,
    //                                             -2.01363397e-01,  2.14521855e-01, -2.19644094e-01,
    //                                              2.12536909e-01,  2.12536909e-01, -1.93293536e-01};


    //     // scalarList PC2_Coef(21);
    //     std::vector<float> PC2_Coef_SEP =  {0.0 , -3.10667434e-01,  1.64116055e-03, -3.32279120e-01,
    //                                              7.07945518e-02,  7.06067963e-02, -3.38522248e-01,
    //                                              9.38708636e-02,  3.73264185e-02, -1.36784245e-01,
    //                                             -3.40273942e-01,  1.04402282e-01,  2.28540667e-02,
    //                                             -1.19414761e-01,  1.87050452e-01, -3.40405798e-01,
    //                                              1.10250461e-01,  1.48061924e-02, -1.09713269e-01,
    //                                              1.78667339e-01, -2.16418816e-01, -3.39892893e-01,
    //                                              1.13907811e-01,  9.56852155e-03, -1.03226515e-01,
    //                                              1.72454529e-01,  1.72454529e-01, -2.12910228e-01};

    //     std::vector<float> Theta_SEP(Mean_Funcs_SEP.size(), 0.0);
    //     for (int i = 0; i < Theta_SEP.size(); i++)
    //     {
    //         Theta_SEP[0] = Theta_SEP[0] - (C1_.value()*PC1_Coef_SEP[i] + C2_.value()*PC2_Coef_SEP[i])* Mean_Funcs_SEP[i]/Std_Funcs_SEP[i];
    //         Theta_SEP[i] = (C1_.value()*PC1_Coef_SEP[i] + C2_.value()*PC2_Coef_SEP[i])/Std_Funcs_SEP[i];
    //     }
    //     Theta_SEP[0] = Theta_SEP[0] + C0_.value();
    //     // Info << "Theta_SEP: " << Theta_SEP << endl;


    //     alpha1_separation =   Theta_SEP[0]
    //                         + Theta_SEP[1]*i1_ + Theta_SEP[2]*i2_
    //                         + Theta_SEP[3]*pow(i1_,2.0) + Theta_SEP[4]*i1_*i2_ + Theta_SEP[5]*pow(i2_,2.0)
    //                         + Theta_SEP[6]*pow(i1_,3.0) + Theta_SEP[7]*pow(i1_,2.0)*i2_ + Theta_SEP[8]*i1_*pow(i2_,2.0) + Theta_SEP[9]*pow(i2_,3.0)
    //                         + Theta_SEP[10]*pow(i1_,4.0) + Theta_SEP[11]*pow(i1_,3.0)*i2_ + Theta_SEP[12]*pow(i1_,2.0)*pow(i2_,2.0) + Theta_SEP[13]*i1_*pow(i2_,3.0) + Theta_SEP[14]*pow(i2_,4.0)
    //                         + Theta_SEP[15]*pow(i1_,5.0) + Theta_SEP[16]*pow(i1_,4.0)*i2_ + Theta_SEP[17]*pow(i1_,3.0)*pow(i2_,2.0) + Theta_SEP[18]*pow(i1_,2.0)*pow(i2_,3.0) + Theta_SEP[19]*i1_*pow(i2_,4.0) + Theta_SEP[20]*pow(i2_,5.0)
    //                         + Theta_SEP[21]*pow(i1_,6.0) + Theta_SEP[22]*pow(i1_,5.0)*i2_ + Theta_SEP[23]*pow(i1_,4.0)*pow(i2_,2.0) + Theta_SEP[24]*pow(i1_,3.0)*pow(i2_,3.0) + Theta_SEP[25]*pow(i1_,2.0)*pow(i2_,4.0) + Theta_SEP[26]*i1_*pow(i2_,5.0) + Theta_SEP[27]*pow(i2_,6.0);

    //     // separationFactor_ = pow(scalar(1)-nut*omega_/k_,SeparationPower_.value())*(alpha1_separation);
    // }

    separationFactor_ = pow(
        max
        (
            min
            (
                (scalar(1)-pow(nut*omega_/k_, separationLambda1_.value())),
                scalar(1)
            ),
            scalar(0)
        ),separationLambda2_.value())*(alpha1_separation);







//  Secondary Flow Calculations

    volScalarField alpha2_secondary(I1_ * 0.0);

    // if (secondaryMode_.value() == 1.0)
    // {
    //     if (A0_.value() == 0.0 && A1_.value() == 0.0 && A2_.value() == 0.0)
    //     {
    //         A0_ = -1.653279626929183;
    //         A1_ = 0.6251856319378204;
    //         A2_ = 0.9998609367402819;
    //     }
        alpha2_secondary = A0_
                 + A1_ * (I1_- (4.13641572e-02)) / 9.70441569e-03
                 + A2_ * (I2_ - (-4.13023579e-02)) / 9.75952414e-03;

    // }

    // if (secondaryMode_.value() == 2.0)
    // {
    //     if (A0_.value() == 0.0 && A1_.value() == 0.0 && A2_.value() == 0.0)
    //     {
    //         A0_ = -1.61278;
    //         A1_ = 0.0744613;
    //         A2_ = 0.0152856;
    //     }

    //     std::vector<float> Mean_Funcs ={0.0 , 4.07361637e-02, -4.06622382e-02, -2.26840875e-04,  7.56121923e-05,
    //                                        -8.86716851e-04,  1.78739817e-03,  1.78171730e-03,  3.21617000e-07,
    //                                         3.57357756e-08,  8.96210799e-07, -1.78451382e-03, -1.10427626e-05,
    //                                         3.68087574e-06, -3.97537088e-05,  1.11014338e-05, -3.70043344e-06,
    //                                         3.96981235e-05, -1.07206496e-07,  2.65303338e-07, -8.84337399e-08};

    //     // scalarList Std_Funcs(21);
    //     std::vector<float> Std_Funcs ={1.0 , 1.13120792e-02, 1.13269452e-02, 5.19769389e-04, 1.73258685e-04,
    //                                         3.31578083e-04, 6.64376406e-04, 6.65801773e-04, 3.76274007e-07,
    //                                         4.18091615e-08, 4.55532854e-07, 6.64990400e-04, 2.25476858e-05,
    //                                         7.51597840e-06, 1.77336316e-05, 2.24959571e-05, 7.49873406e-06,
    //                                         1.77460355e-05, 1.25426076e-07, 4.96889902e-07, 1.65631535e-07};

    //     // scalarList PC1_Coef(21);
    //     std::vector<float> PC1_Coef =  {0.0 , -2.02243623e-01,  2.03704260e-01,  2.38820097e-01,
    //                                         -2.38822527e-01,  2.30314951e-01, -2.30856619e-01,
    //                                         -2.31949199e-01, -1.30774177e-01, -1.30769020e-01,
    //                                         -2.43145321e-01,  2.31445084e-01,  2.45458420e-01,
    //                                         -2.45459807e-01,  2.40217847e-01, -2.44779753e-01,
    //                                          2.44781156e-01, -2.40374857e-01,  1.30771600e-01,
    //                                         -2.48576040e-01,  2.48577012e-01};


    //     // scalarList PC2_Coef(21);
    //     std::vector<float> PC2_Coef =  {0.0 , -2.65008897e-01,  2.61415953e-01, -2.33956132e-01,
    //                                          2.33951526e-01,  2.51420539e-01, -2.52546533e-01,
    //                                         -2.47383736e-01,  1.24902307e-01,  1.24895664e-01,
    //                                         -2.04347009e-01,  2.49984607e-01, -2.30083943e-01,
    //                                          2.30081804e-01,  2.27204329e-01,  2.29736516e-01,
    //                                         -2.29734398e-01, -2.25321733e-01, -1.24898987e-01,
    //                                          2.22372343e-01, -2.22371319e-01};




    //     std::vector<float> Theta_(21, 0.0);
    //     for (int i = 0; i < Theta_.size(); i++)
    //     {
    //         Theta_[0] = Theta_[0] - (A1_.value()*PC1_Coef[i] + A2_.value()*PC2_Coef[i])* Mean_Funcs[i]/Std_Funcs[i];
    //         Theta_[i] = (A1_.value()*PC1_Coef[i] + A2_.value()*PC2_Coef[i])/Std_Funcs[i];
    //     }
    //     Theta_[0] = Theta_[0] + A0_.value();
    //     //Info << "Theta_: " << Theta_ << endl;


    //     alpha2_secondary =    Theta_[0]
    //         + Theta_[1]*i1_ + Theta_[2]*i2_ + Theta_[3]*i3_ + Theta_[4]*i4_ + Theta_[5]*i5_
    //         + Theta_[6]*i1_*i1_ + Theta_[7]*i2_*i2_ + Theta_[8]*i3_*i3_ + Theta_[9]*i4_*i4_ + Theta_[10]*i5_*i5_
    //         + Theta_[11]*i1_*i2_ + Theta_[12]*i1_*i3_ + Theta_[13]*i1_*i4_ + Theta_[14]*i1_*i5_
    //         + Theta_[15]*i2_*i3_ + Theta_[16]*i2_*i4_ + Theta_[17]*i2_*i5_
    //         + Theta_[18]*i3_*i4_ + Theta_[19]*i3_*i5_
    //         + Theta_[20]*i4_*i5_;
    // }

    bijDelta_ = bijDelta_ + ((nut*omega_/(k_ + this->kMin_))*(alpha2_secondary*Tij2_) - bijDelta_)*secondary_relaxation_;

    Rij_ = ((2.0/3.0)*I)*k_ - 2.0*nut*Sij + 2*k_*bijDelta_;
    volSymmTensorField dAij(2*k_*bijDelta_);
    volSymmTensorField P(-twoSymm(dAij & gradU));
    volScalarField Pk_bijDelta_(0.5*tr(P));



//  Continue with kOSST


    // volScalarField::Internal GbyNu0
    // (
    //     this->type() + ":GbyNu",
    //     // ((tgradU() && dev(twoSymm(tgradU()))) )
    //     ( (tgradU() && dev(twoSymm(tgradU()))) )
    // );

    volScalarField::Internal GbyNu0
    (
        this->type() + ":GbyNu",
        // ((tgradU() && dev(twoSymm(tgradU()))) )
        ( (tgradU() && dev(twoSymm(tgradU()))) + Pk_bijDelta_/(nut + nutMin) )
    );


    volScalarField::Internal G(this->GName(), nut*GbyNu0);
    // const fvMesh& mesh1 = this->mesh_;
    // volSymmTensorField& bijDelta = mesh1.lookupObjectRef<volSymmTensorField>("bijDelta");

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();




    volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());

    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        GbyNu0 = GbyNu(GbyNu0, F23(), S2());
        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
         ==
            // alpha()*rho()*gamma*GbyNu0
            alpha()*rho()*gamma*(GbyNu0+GbyNu0*separationFactor_)  // Cw1 F1 omega/k Pk
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
          - fvm::Sp(alpha()*rho()*beta*omega_(), omega_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
                omega_
            )
          + alpha()*rho()*beta*sqr(omegaInf_)
          + Qsas(S2(), gamma, beta)
          + omegaSource()
          + fvOptions(alpha, rho, omega_)
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);
    }


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(F1), k_)
     ==
        alpha()*rho()*Pk(G)
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilonByk(F1, tgradU()), k_)
      + alpha()*rho()*betaStar_*omegaInf_*kInf_
      + kSource()
      + fvOptions(alpha, rho, k_)
    );


    tgradU.clear();

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut(S2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
