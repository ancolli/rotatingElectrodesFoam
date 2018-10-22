Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//From here

Info<< "Reading field C\n" << endl;
    volScalarField C
    (
        IOobject
        (
            "C",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    volScalarField Flux
    (
        IOobject
        (
            "Flux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE//AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("0", C.dimensions()/dimTime*dimLength, 0.0)
    );

Info<< "Reading field Deff\n" << endl;

    volScalarField Deff
    (
        IOobject
        (
            "Deff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,//MUST_READ,/
            IOobject::AUTO_WRITE//NO_WRITE
        ),
        mesh,
        dimensionedScalar("Deff",dimensionSet(0,2,-1,0,0,0,0),0)
    );

Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );


 Info<< "Reading diffusivity D\n" << endl;

    dimensionedScalar Dmol
    (
        transportProperties.lookup("Dmol")
    );

    Info<< "Reading turbulent Schmidt number\n" << endl;

    dimensionedScalar Sct
    (
        transportProperties.lookup("Sct")
    );

//Until here

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"

