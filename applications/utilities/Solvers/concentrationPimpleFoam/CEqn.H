
{

    Deff = Dmol + turbulence->nut()/Sct;
    
    fvScalarMatrix CEqn
    (
 		fvm::ddt(C)
              + fvm::div(phi, C)
              - fvm::laplacian(Deff, C)
             ==
                fvOptions(C)
    );

    CEqn.relax();

    fvOptions.constrain(CEqn);

    CEqn.solve();

    fvOptions.correct(C);

    if(runTime.outputTime())
    {
    	Deff.write();	
    }


// Calculates flux on cathode
        //Info<< "Calculating Flux" << endl;

	label patchi = mesh.boundaryMesh().findPatchID("wall_cathode");

        Flux.boundaryFieldRef()[patchi] = -1*Deff.boundaryFieldRef()[patchi]*C.boundaryFieldRef()[patchi].snGrad();

			    
// Calculates average Flux

	scalar area = gSum(mesh.magSf().boundaryField()[patchi]);
        scalar sumField = 0;

        if (area > 0)
        {
        	sumField = gSum
                (
                	mesh.magSf().boundaryField()[patchi]
                      * Flux.boundaryField()[patchi]
                )/ area;
        }

        Info<< "Average flux over cathode = " << sumField << endl;
}
