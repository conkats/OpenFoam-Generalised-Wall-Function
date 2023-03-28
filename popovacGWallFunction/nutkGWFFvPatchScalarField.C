/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "nutkGWFFvPatchScalarField.H"
#include "momentumTransportModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> nutkGWFFvPatchScalarField::nut() const
                 
{
    const label patchi = patch().index(); // e.g.looking at a patch at index i lower wall patch 

    const momentumTransportModel& turbModel =
        db().lookupObject<momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                internalField().group()
            )
        );
        
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];// velocity at the patch//reference is an alias
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));
    const scalarField magGradUw(mag(Uw.snGrad()));
    

    const scalarField& y = turbModel.y()[patchi];// distance to the wall from the patch
    const tmp<volScalarField> tk = turbModel.k();//temporary scalar-efficient way of handling and defining tk for k 
    const volScalarField& k = tk();
    const tmp<scalarField> tnuw = turbModel.nu(patchi);//temporary new nu
    const scalarField& nuw = tnuw();

    const scalar Cmu25 = pow025(Cmu_);

    tmp<scalarField> tnutw(new scalarField(patch().size(), 0.0));
    scalarField& nutw = tnutw.ref();
    
    
    /*Generalised Wall Function*/
    
    volVectorField vU = turbModel.U(); // take Velocity Field
    
    const labelList cells = patch().faceCells();//the cells adjacent to patch
    //const labelList& cells;
    tmp<vectorField> tnv = patch().Sf()/patch().magSf(); // normal face vector

    vectorField& nv = tnv.ref();//reference-access variable and value to volume vector & to tnv - alias
    vectorField tv = nv;       //init tangential vector
    scalarField alpha = nutw;  // init alpha, angle between coordinate systems
    scalarField beta = nutw;
    scalarField gamma = nutw;
    scalarField Ut = nutw;
    scalarField Un = nutw;
    scalarField txn = nutw;

    ///////////////////////////////////////////////////////
   /*forAll(nv,celli)// loop over the cells of patch, total is cells
    
    {
        if (nv[celli].component(2) < 0)// if the z component is less than zero
        {
            tv[celli].component(0) = -nv[celli].component(2);//the z component of the normal component
            tv[celli].component(2) = nv[celli].component(0);
            
        
        }//u.n=0
        if (nv[celli].component(2) >= 0)
        {
            tv[celli].component(0) = nv[celli].component(2);
            tv[celli].component(2) = -nv[celli].component(0);
            
        
        }//tangential vector is 1
        //DEBUGGING vectors
         //tv[celli].component(0) = 1;
        //tv[celli].component(2) = 1;//-nv[celli].component(0);
        //Info<< "ComponentX = " << tv[celli].component(0) << endl;//output in log file
        
    //Ut=U-U.n.n
    
    alpha[celli] = atan2(tv[celli].component(2),tv[celli].component(0));
    beta[celli] = atan2(vU[cells[celli]].component(2),(vU[cells[celli]].component(0)));
    // angle between velocity and x-axis
    gamma[celli] = beta[celli] - alpha[celli]; // angle between velocity and wall
    //Ut[celli] = magUp[celli].component[0]=magUp[celli]

    Ut[celli] = magUp[celli] * cos(gamma[celli]); // Wall Velocity Component
    Un[celli] = magUp[celli] * sin(gamma[celli]); // Normal Wall Velocity Component

    //Ut[celli] = magUp[celli] ;
    //Un[celli] = magUp[celli] * nv[celli].component(2);
    
    
    }
     //Needs to be checked since cannot find gradU from lookupObject 
    //using lookupObjectFunction
    //const volVectorField& GradU = this->db().lookupObject<volVectorField>("gradU");
    //const volVectorField& GradV = this->db().lookupObject<volVectorField>("gradV");
    //const volVectorField& GradP = this->db().lookupObject<volVectorField>("gradP"); // take Pressure Field
    //volVectorField dUd = GradU;
    //volVectorField dVd = GradV;
    
    //tested but still complains!!!
    const volVectorField& U = this->db().lookupObject<volVectorField>("U"); //efficient way
    const volScalarField& P = this->db().lookupObject<volScalarField>("p"); // take Pressure Field
    //calculate the gradien
    volTensorField dUd = fvc::grad(U);
    //volTensorField dVd = fvc::grad(V);
    volVectorField GradP = fvc::grad(P);
    vectorField    dUtd = nv;
    
    scalarField dUtdt = nutw;
    scalarField dUtdn = nutw;
    scalarField dPdt = nutw;
    scalarField Cu = nutw;
    //access P at the boundary
    vectorField gradP = GradP.boundaryField()[patchi].patchInternalField();

    scalarField nutw2 = nutw;
    //for [celli] loop of each celli in the domain and compute the components
    // the cells[celli] loops over the boundary cells
    forAll(dUtdt,celli)
    {
    
//       dUtd[celli].component(0) = cos(alpha[celli]) *         dUd[cells[celli]].component(0) + sin(alpha[celli]) * dVd[cells[celli]].component(0); // dUw / dx

//       dUtd[celli].component(2) = cos(alpha[celli]) * dUd[cells[celli]].component(2) + sin(alpha[celli]) * dVd[cells[celli]].component(2); // dUw / dy

       // ref:https://www.cfd-online.com/Forums/openfoam-programming-development/230642-arrangement-components-gradu.html
       //  component of gradU  
       //  gradU.compoent(0)=dudx
       //  gradU.compoent(1)=dvdx
       //  gradU.compoent(2)=dwdx
       //  gradU.compoent(3)=dudy
       // gradU.compoent(4)=dvdy
       //  gradU.compoent(5)=dwdy
       //  gradU.compoent(6)=dudz
       //  gradU.compoent(7)=dvdz
       //  gradU.compoent(8)=dwdz
       //  can be also access the component of grad(U) using component(tensor:XX) i.e
       //  dUd[cells[celli]].component(tensor::XX)
       
       dUtd[celli].component(0) = cos(alpha[celli]) *         dUd[cells[celli]].component(0) + sin(alpha[celli]) * dUd[cells[celli]].component(1); // dUw / dx

       dUtd[celli].component(2) = cos(alpha[celli]) * dUd[cells[celli]].component(2) + sin(alpha[celli]) * dUd[cells[celli]].component(4); // dUw / dy

       dUtdt[celli] = dUtd[celli] & tv[celli]; // dUw /dw 
       dUtdn[celli] = dUtd[celli] & (-nv[celli]); // dUw / dz
    
       dPdt[celli] =gradP[celli] & tv[celli];// dP/ dw
       Cu[celli] = dUtdt[celli] * Ut[celli] + dUtdn[celli] * Un[celli] + dPdt[celli];

    }
*/
    ///////////////////////////////////////////////////////
    //Code_Saturne version

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
       //Info<< "Banana3 = " <<  dUd[cells[celli]]<< endl;//output in log file
       //Info<< "Banana4 = " <<  gradP[celli]<< endl;//output in log file
       Cu[celli] = dPdt[celli]+ Ut[celli]*dUtdt[celli];

    }

    ///////////////////////////////////////////////////////

    forAll(nutw, facei)//patch
    {
        label celli = patch().faceCells()[facei];

        scalar yPlus = Cmu25*y[facei]*sqrt(k[celli])/nuw[facei];
        
        // Non dimensionless Psi Calculation
        scalar kstar = kappa_*Cmu25;
        scalar uk = Cmu25*sqrt(k[celli]);
        //scalar utau = sqrt(nuw[facei]*mag(dUtdn[facei]))+1e-12;
        scalar Ustar = mag(Ut[facei]) /uk;

        scalar Custar = Cu[facei] * nuw[facei] / (sqr(uk)*sqrt(k[celli]));
        scalar psi = 1 - Custar * yPlus / (kstar * Ustar );
        if (psi <= 0.1)
        {
           psi = 0.1;
        }
        
        // Compound Wall Treatment
        scalar Gamma = 0.01 * pow(yPlus,4) / (1 + 5 * yPlus);
        nutw[facei] = nuw[facei]*(yPlus*kappa_*psi/log(E_*yPlus) - 1.0)*Foam::exp(-1/Gamma);
        //if (yPlus > yPlusLam_)
        //{
        //    nutw[facei] = nuw[facei]*(yPlus*kappa_/log(E_*yPlus) - 1.0);
        //}
    }

    return tnutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutkGWFFvPatchScalarField::nutkGWFFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF)
{}


nutkGWFFvPatchScalarField::nutkGWFFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict)
{}


nutkGWFFvPatchScalarField::nutkGWFFvPatchScalarField
(
    const nutkGWFFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


nutkGWFFvPatchScalarField::nutkGWFFvPatchScalarField
(
    const nutkGWFFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(wfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> nutkGWFFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();

    const momentumTransportModel& turbModel =
        db().lookupObject<momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                internalField().group()
            )
        );

    const scalarField& y = turbModel.y()[patchi];

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    tmp<scalarField> kwc = k.boundaryField()[patchi].patchInternalField();
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();
    return pow025(Cmu_)*y*sqrt(kwc)/nuw;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutkGWFFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
