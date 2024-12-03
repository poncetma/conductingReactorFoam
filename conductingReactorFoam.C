/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Application
    laplacianFoam

Group
    grpBasicSolvers

Description
    Laplace equation solver for a scalar quantity.

    \heading Solver details
    The solver is applicable to, e.g. for thermal diffusion in a solid.  The
    equation is given by:

    \f[
        \ddt{T}  = \div \left( D_T \grad T \right)
    \f]

    Where:
    \vartable
        T     | Scalar field which is solved for, e.g. temperature
        D_T   | Diffusion coefficient
    \endvartable

    \heading Required fields
    \plaintable
        T     | Scalar field which is solved for, e.g. temperature
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	/*
    argList::addNote
    (
        "Laplace equation solver for a scalar quantity."
    );
	*/

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;
	
	//Take some inspiration from GeN-Foam here (eigenvalue vs transient neutronics selection, see the neutronics class)
	bool isTransient;
	
	IOdictionary simulationType
    (
        IOobject
        (
            "simulationType",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

	isTransient = simulationType.lookupOrDefault("isTransient", false);
	if (isTransient) {
		Info << "We have a transient simulation." << nl << endl;
	}
	else {
		Info << "Running a steady-state simulation." << nl << endl;
	}
	
	//read the desired nominal power from the same input dict 
	double nominal_power;	
	nominal_power = simulationType.lookupOrDefault("nominalPower", 0.0);
	Info << "Read the nominal power: " << nominal_power << " W" << nl << endl;
	
	
	//Compute total mesh volume
	scalar total_volume(0.0);
	forAll(mesh.cells(),cellI)
	{
		total_volume += mesh.V()[cellI];
	}
	Info << "The volume of the mesh (fuel) is:" << total_volume << " m3" << nl << endl;


	double integral_power;
	//This loop performs multiple solves of the heat equation at each time-step, since the input of fvOptions becomes an explicit source
    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
		
        while (simple.correctNonOrthogonal()) //Give the option of non-orthogonal correction iterations. Defaults to just 1 iteration.
        {
			//solve a dummy "heat source equation" which just takes an input from fvOptions 
			fvScalarMatrix qsEqn
			(
				fvm::ddt(q_s) //OF knows to make this equal to 0				
			);
			fvOptions.correct(q_s); 
			
			/*			
			q_s should be rescaled based on the nominal power (do this directly in the source term).
			FF will pass a power density field which is normalised to the initial power (as computed in a proportional but inexact way by FF),
			so it can be scaled by the true nominal power to recover the true power density. 	
			*/
			
			//check that the integral power is properly recovered (should match the nominal power at t=0)
			integral_power = 0.0; 
			forAll(mesh.cells(),cellI)
			{
				integral_power += mesh.V()[cellI]*(q_s[cellI]*nominal_power/total_volume);
			}			
			Info << "The integral power is:" << integral_power << " W" << nl << endl;
			
			//solve the heat equation
			fvScalarMatrix TEqn
            (
                -fvm::laplacian(DT, T) == (q_s*nominal_power/total_volume)*DT/k //Here's the chance to amplify the field given by the neutronics solver				
            );
			if (isTransient) {
				TEqn += fvm::ddt(T);
			}			
            
            SolverPerformance<double> solverPerf = TEqn.solve();             
			
			//Try manually fetching and printing the residual
			Info<< "Checking solver performance..." << nl << endl;
			Info<< "T residual = " << solverPerf.finalResidual() << nl << endl;
        }

        #include "write.H"

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
