/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
	    fvScalarMatrix TEqn
            (
	        fvm::ddt(T) - fvm::laplacian(DT, T)
	    );
	    
	    TEqn.solve();
	    
	    List<List<scalar> > A;
	    List<scalar> b;
                
	    A.resize(T.size());
	    b.resize(T.size());
	    forAll(A, i)
	    {
            	    A[i].resize(T.size());
            	    forAll(A[i],j)
            	    {
            		A[i][j] = 0.0;
            	    }//for j
            	    b[i] = 0.0;
	    }//for i
                
	    forAll(A,i)
	    {
            	    A[i][i] = TEqn.diag()[i];
            	    b[i]    = TEqn.source()[i];
	    }
                
	    const lduAddressing& addr = TEqn.lduAddr();
	    const labelList& lowerAddr = addr.lowerAddr();
	    const labelList& upperAddr = addr.upperAddr();                
                
	    forAll(lowerAddr, i)
	    {
            	    A[lowerAddr[i]][upperAddr[i]] = TEqn.upper()[i];
            	    A[upperAddr[i]][lowerAddr[i]] = TEqn.lower()[i];           	                	    
	    }
                
	    forAll(T.boundaryField(),I)
	    {
            	    const fvPatch &ptch=T.boundaryField()[I].patch();
            	    forAll(ptch,J)
            	    {
            		int w=ptch.faceCells()[J];
            		A[w][w]+=TEqn.internalCoeffs()[I][J];
            		b[w]   +=TEqn.boundaryCoeffs()[I][J];            		
            	    }
                
	    }
                
	    //Info << "A=" <<nl<<A<<nl<<endl;
	    //Info << "b=" <<nl<<b<<nl<<endl;
                
	    OFstream fileA("mA.txt");
	    fileA.precision(18);
	    forAll(A,i)
	    {
		    if (i!=0)
			fileA << endl;
		    forAll(A[i],j)
		    {
		        if (j!=0)
		    	    fileA<<" ";
		        fileA<<A[i][j];
		    }

	    }              
	    //fileA.close();

	    OFstream fileb("mb.txt");
	    fileb.precision(18);
	    forAll(b,i)
	    {
		    if (i!=0)
			fileb << " ";
		    fileb << b[i]; 
	    }              

        }

        #include "write.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
