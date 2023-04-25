/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
#include "epsilonGWFFvPatchScalarField.H"
#include "nutkGWFFvPatchScalarField.H"
#include "momentumTransportModel.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "fvCFD.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::epsilonGWFFvPatchScalarField::tolerance_ = 1e-5;

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::epsilonGWFFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<epsilonGWFFvPatchScalarField>(bf[patchi]))
        {
            epsilonGWFFvPatchScalarField& epf = epsilonPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            epf.master() = master;
        }
    }
}


void Foam::epsilonGWFFvPatchScalarField::createAveragingWeights()
{
    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    const fvMesh& mesh = epsilon.mesh();

    if (initialised_ && !mesh.changing())
    {
        return;
    }

    volScalarField weights
    (
        IOobject
        (
            "weights",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false // do not register
        ),
        mesh,
        dimensionedScalar(dimless, 0)
    );

    DynamicList<label> epsilonPatches(bf.size());
    forAll(bf, patchi)
    {
        if (isA<epsilonGWFFvPatchScalarField>(bf[patchi]))
        {
            epsilonPatches.append(patchi);

            const labelUList& faceCells = bf[patchi].patch().faceCells();
            forAll(faceCells, i)
            {
                weights[faceCells[i]]++;
            }
        }
    }

    cornerWeights_.setSize(bf.size());
    forAll(epsilonPatches, i)
    {
        label patchi = epsilonPatches[i];
        const fvPatchScalarField& wf = weights.boundaryField()[patchi];
        cornerWeights_[patchi] = 1.0/wf.patchInternalField();
    }

    G_.setSize(internalField().size(), 0.0);
    epsilon_.setSize(internalField().size(), 0.0);

    initialised_ = true;
}


Foam::epsilonGWFFvPatchScalarField&
Foam::epsilonGWFFvPatchScalarField::epsilonPatch(const label patchi)
{
    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    const epsilonGWFFvPatchScalarField& epf =
        refCast<const epsilonGWFFvPatchScalarField>(bf[patchi]);

    return const_cast<epsilonGWFFvPatchScalarField&>(epf);
}


void Foam::epsilonGWFFvPatchScalarField::calculateTurbulenceFields
(
    const momentumTransportModel& turbulence,
    scalarField& G0,
    scalarField& epsilon0
)
{
    // Accumulate all of the G and epsilon contributions
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            epsilonGWFFvPatchScalarField& epf = epsilonPatch(patchi);

            const List<scalar>& w = cornerWeights_[patchi];

            epf.calculate(turbulence, w, epf.patch(), G0, epsilon0);
        }
    }

    // Apply zero-gradient condition for epsilon
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            epsilonGWFFvPatchScalarField& epf = epsilonPatch(patchi);

            epf == scalarField(epsilon0, epf.patch().faceCells());
        }
    }
}


void Foam::epsilonGWFFvPatchScalarField::calculate
(
    const momentumTransportModel& turbModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G0,
    scalarField& epsilon0
)
{
    const label patchi = patch.index();// e.g.looking at a patch at index i lower wall patch

    const nutWallFunctionFvPatchScalarField& nutw = nutWallFunctionFvPatchScalarField::nutw(turbModel, patchi);
    const scalarField& y = turbModel.y()[patchi];// distance to the wall from the patch

    const tmp<scalarField> tnuw = turbModel.nu(patchi);//temporary new nu
    const scalarField& nuw = tnuw();

    const tmp<volScalarField> tk = turbModel.k();//temporary scalar-efficient way of handling and defining tk for k
    const volScalarField& k = tk();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];// velocity at the patch//reference is an alias

    const scalarField magGradUw(mag(Uw.snGrad()));

    const scalar Cmu25 = pow025(nutw.Cmu());
    const scalar Cmu75 = pow(nutw.Cmu(), 0.75);
    const scalar Cmu_  = 0.09;
    const scalar kappa_= 0.41;
    // Set epsilon and G
    /*Generalised Wall Function*/
    //const labelList cells = patch().faceCells();// get adjacent cell list for a face
    //https://www.cfd-online.com/Forums/openfoam-meshing/82881-getting-cells-next-patch.html
    const labelList cells = patch.faceCells();//the cells adjacent to patch
    tmp<vectorField> tnv = patch.Sf()/patch.magSf(); // normal face vector
    
    const volVectorField vU = turbModel.U(); // take Velocity Field
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));
    vectorField& nv = tnv.ref();//reference-access variable and value to volume vector & to tnv - alias
    vectorField tv = nv;       //init tangential vector
    scalarField alpha = nutw;  // init alpha, angle between coordinate systems
    scalarField beta = nutw;
    scalarField gamma = nutw;
    scalarField Ut = nutw;
    scalarField Un = nutw;
    scalarField txn = nutw;

    ///////////////////////////////////////////////////////

    const volVectorField& U = this->db().lookupObject<volVectorField>("U");
    const volScalarField& P = this->db().lookupObject<volScalarField>("p"); // take Pressure Field
    volTensorField dUd = fvc::grad(U);
    volVectorField GradP=fvc::grad(P);
        

    vectorField dUtd = nv;
    scalarField dUtdt = nutw;
    scalarField dUtdn = nutw;
    scalarField dPdt = nutw;
    scalarField Cu = nutw;
    
    //access P at the boundary 
    vectorField gradP = GradP.boundaryField()[patchi].patchInternalField();
    scalarField nutw2 = nutw;

    forAll(nv,celli)//for(int celli=0; i<dUtdt.size(); celli++)

    {
        //Un=Upx*unx+Upy*uny+Upz*unz
        Un[celli] = vU[cells[celli]].component(0) * int(nv[celli].component(0))+vU[cells[celli]].component(1) * int(nv[celli].component(1))+ vU[cells[celli]].component(2) * int(nv[celli].component(2)); //Normal Wall Velocity ComponentNormal Wall Velocity Component
       //tx = Upx-Upx*unx
       tv[celli].component(0) = vU[celli].component(0)-Un[celli]*nv[celli].component(0);
       tv[celli].component(1) = vU[celli].component(1)-Un[celli]*nv[celli].component(1);
       tv[celli].component(2) = vU[celli].component(2)-Un[celli]*nv[celli].component(2);

       //mag
       txn[celli] = sqrt(pow(tv[celli].component(0),2)+pow(tv[celli].component(1),2)+pow(tv[celli].component(2),2));

       tv[celli].component(0) = int(tv[celli].component(0)/txn[celli]);
       tv[celli].component(1) = int(tv[celli].component(1)/txn[celli]);
       tv[celli].component(2) = int(tv[celli].component(2)/txn[celli]);

       Ut[celli] = vU[celli].component(0) * tv[celli].component(0) + vU[celli].component(1) * tv[celli].component(1)+vU[celli].component(2) * tv[celli].component(2); 
       
      //dUtdt[celli] = Ut[celli] * (
      //           tv[celli].component(0) * (tv[celli].component(0) * dUd[cells[celli]].component(0) + tv[celli].component(1) *dUd[cells[celli]].component(3) + tv[celli].component(2)* dUd[cells[celli]].component(6)) +
      //           tv[celli].component(1) * (tv[celli].component(0) * dUd[cells[celli]].component(1) + tv[celli].component(1) * dUd[cells[celli]].component(4) + tv[celli].component(2) * dUd[cells[celli]].component(7)) +
      //           tv[celli].component(2) * (tv[celli].component(0)* dUd[cells[celli]].component(2) + tv[celli].component(1)* dUd[cells[celli]].component(5) + tv[celli].component(2)* dUd[cells[celli]].component(8))
      //                 );

       //dPdt[celli] = gradP[celli].component(0) * tv[celli].component(0) + gradP[celli].component(1) *tv[celli].component(1) +gradP[celli].component(2) *tv[celli].component(2);// dP/ dw
       dPdt[celli] =gradP[celli] & tv[celli];

       dUtd[celli].component(0) =tv[celli].component(0) * (tv[celli].component(0) * dUd[celli].component(0) + tv[celli].component(1) *dUd[celli].component(3) + tv[celli].component(2)* dUd[celli].component(6));
       dUtd[celli].component(1) = tv[celli].component(1) * (tv[celli].component(0) * dUd[celli].component(1) + tv[celli].component(1) * dUd[celli].component(4) + tv[celli].component(2) * dUd[celli].component(7));
       dUtd[celli].component(2) = tv[celli].component(2) * (tv[celli].component(0)* dUd[celli].component(2) + tv[celli].component(1)* dUd[celli].component(5) + tv[celli].component(2)* dUd[celli].component(8));
       Cu[celli] = dPdt[celli]+ Ut[celli]*dUtdt[celli];


    }


    ///////////////////////////////////////////////////////
    forAll(nutw, facei)// loop over patch faces
    {
        const label celli = patch.faceCells()[facei];

        const scalar yPlus = Cmu25*y[facei]*sqrt(k[celli])/nuw[facei];

        const scalar w = cornerWeights[facei];

        //if (yPlus > nutw.yPlusLam())
        //{
        //    epsilon0[celli] +=
        //        w*Cmu75*pow(k[celli], 1.5)/(nutw.kappa()*y[facei]);

        //    G0[celli] +=
        //        w
        //       *(nutw[facei] + nuw[facei])
        //       *magGradUw[facei]
        //       *Cmu25*sqrt(k[celli])
        //      /(nutw.kappa()*y[facei]);
        //}
        //else
        //{
        //    epsilon0[celli] += w*2.0*k[celli]*nuw[facei]/sqr(y[facei]);
        //}

        scalar kstar = kappa_*Cmu25;

        scalar uk = Cmu25*sqrt(k[celli]);
        //scalar utau = sqrt(nuw[facei]*mag(dUtdn[facei]))+1e-12;
        scalar Ustar = mag(Ut[facei]) /uk;
        scalar Custar = Cu[facei] * nuw[facei] / (sqr(uk)*sqrt(k[celli]));
        scalar psi = 1 - Custar * yPlus / (kstar * Ustar );
        if (psi <= 0.1 )
        {
            psi = 0.1;
        }
        
        // Compound Wall Treatment
        scalar Gamma = 0.01 * pow(yPlus,4) / (1 + 5 * yPlus);
        scalar Gamma2 = 0.001 * pow(yPlus,4) / (1 + yPlus);
        epsilon0[celli] += w*(2.0*k[celli]*nuw[facei]/sqr(y[facei]) * Foam::exp(-Gamma2)+ Cmu75*pow(k[celli], 1.5)/(psi*kappa_*y[facei]) * Foam::exp(-1/Gamma2));
        G0[celli] += w*(Cmu_*sqr(k[celli])/epsilon0[celli]*sqr(magGradUw[facei])*Foam::exp(-Gamma)+(nutw[facei] + nuw[facei])*magGradUw[facei]*Cmu25*sqrt(k[celli])/(nutw.kappa()*y[facei]) *
Foam::exp(-1/Gamma));
        
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::epsilonGWFFvPatchScalarField::
epsilonGWFFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    G_(),
    epsilon_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{}


Foam::epsilonGWFFvPatchScalarField::
epsilonGWFFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    G_(),
    epsilon_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    // Apply zero-gradient condition on start-up
    this->operator==(patchInternalField());
}


Foam::epsilonGWFFvPatchScalarField::
epsilonGWFFvPatchScalarField
(
    const epsilonGWFFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    G_(),
    epsilon_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{}


Foam::epsilonGWFFvPatchScalarField::
epsilonGWFFvPatchScalarField
(
    const epsilonGWFFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ewfpsf, iF),
    G_(),
    epsilon_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField& Foam::epsilonGWFFvPatchScalarField::G(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            G_ = 0.0;
        }

        return G_;
    }

    return epsilonPatch(master_).G();
}


Foam::scalarField& Foam::epsilonGWFFvPatchScalarField::epsilon
(
    bool init
)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            epsilon_ = 0.0;
        }

        return epsilon_;
    }

    return epsilonPatch(master_).epsilon(init);
}


void Foam::epsilonGWFFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const momentumTransportModel& turbModel =
        db().lookupObject<momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                internalField().group()
            )
        );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), epsilon(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& epsilon0 = this->epsilon();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbModel.GName())
        );

    FieldType& epsilon = const_cast<FieldType&>(internalField());

    forAll(*this, facei)
    {
        label celli = patch().faceCells()[facei];

        G[celli] = G0[celli];
        epsilon[celli] = epsilon0[celli];
    }

    fvPatchField<scalar>::updateCoeffs();
}


void Foam::epsilonGWFFvPatchScalarField::updateWeightedCoeffs
(
    const scalarField& weights
)
{
    if (updated())
    {
        return;
    }

    const momentumTransportModel& turbModel =
        db().lookupObject<momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                internalField().group()
            )
        );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), epsilon(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& epsilon0 = this->epsilon();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbModel.GName())
        );

    FieldType& epsilon = const_cast<FieldType&>(internalField());

    scalarField& epsilonf = *this;

    // Only set the values if the weights are > tolerance
    forAll(weights, facei)
    {
        scalar w = weights[facei];

        if (w > tolerance_)
        {
            label celli = patch().faceCells()[facei];

            G[celli] = (1.0 - w)*G[celli] + w*G0[celli];
            epsilon[celli] = (1.0 - w)*epsilon[celli] + w*epsilon0[celli];
            epsilonf[facei] = epsilon[celli];
        }
    }

    fvPatchField<scalar>::updateCoeffs();
}


void Foam::epsilonGWFFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    matrix.setValues(patch().faceCells(), patchInternalField());

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void Foam::epsilonGWFFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const Field<scalar>& weights
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    DynamicList<label> constraintCells(weights.size());
    DynamicList<scalar> constraintEpsilon(weights.size());
    const labelUList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& epsilon
        = internalField();

    label nConstrainedCells = 0;


    forAll(weights, facei)
    {
        // Only set the values if the weights are > tolerance
        if (weights[facei] > tolerance_)
        {
            nConstrainedCells++;

            label celli = faceCells[facei];

            constraintCells.append(celli);
            constraintEpsilon.append(epsilon[celli]);
        }
    }

    if (debug)
    {
        Pout<< "Patch: " << patch().name()
            << ": number of constrained cells = " << nConstrainedCells
            << " out of " << patch().size()
            << endl;
    }

    matrix.setValues
    (
        constraintCells,
        scalarField(constraintEpsilon)
    );

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        epsilonGWFFvPatchScalarField
    );
}


// ************************************************************************* //
