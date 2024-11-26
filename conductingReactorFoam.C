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
	

	//This loop performs multiple solves of the heat equation at each time-step, since the input of fvOptions becomes an explicit source
    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
		
        while (simple.correctNonOrthogonal()) //Give the option of non-orthogonal correction iterations. Currently defaulting to just 1 iteration
        {
			//solve a dummy "heat source equation" which just takes an input from fvOptions 
			fvScalarMatrix qsEqn
			(
				fvm::ddt(q_s) //OF knows to make this equal to 0
				//fvm::ddt(S)
			);
			fvOptions.correct(q_s); //Should assume that the powershape passed in is normalised. 
			//fvOptions.correct(S); 
			
			//solve the heat equation
			//Try to formulate this in a more general way
			fvScalarMatrix TEqn
            (
                -fvm::laplacian(DT, T) == q_s*DT/k //Here's the chance to amplify the field given by the neutronics solver
				//+ fvOptions(T)  
				//S*DT/k //+ fvOptions(T)  //added the temperature change due to the heat source
            );
			if (isTransient) {
				TEqn += fvm::ddt(T);
			}
			
            //fvOptions.constrain(TEqn);
            SolverPerformance<double> solverPerf = TEqn.solve(); 
            //fvOptions.correct(T);
			
			//Try manually fetching and printing the residual
			Info<< "ATTN" << nl << endl;
			Info<< "T residual = " << solverPerf.finalResidual() << nl << endl;
        }

        #include "write.H"

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
