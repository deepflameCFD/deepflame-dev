#include "readGravitationalAcceleration.H"

Info<< "Reading thermophysical properties\n" << endl;
fluidThermo* pThermo = new hePsiThermo<psiThermo, CanteraMixture>(mesh, word::null);
fluidThermo& thermo = *pThermo;
thermo.validate(args.executable(), "ha");

SLGThermo slgThermo(mesh, thermo);

const volScalarField& psi = thermo.psi();
volScalarField& p = thermo.p();
volScalarField& T = thermo.T();
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

dfChemistryModel<basicThermo> chemistry(thermo);
PtrList<volScalarField>& Y = chemistry.Y();
const word inertSpecie(chemistry.lookup("inertSpecie"));
const label inertIndex(chemistry.species()[inertSpecie]);

chemistry.correctThermo();
Info<< "At initial time, min/max(T) = " << min(T).value() << ", " << max(T).value() << endl;

Info<< "\nReading field U\n" << endl;
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

#include "compressibleCreatePhi.H"

mesh.setFluxRequired(p.name());

dimensionedScalar rhoMax
(
    dimensionedScalar::lookupOrDefault
    (
        "rhoMax",
        pimple.dict(),
        dimDensity,
        great
    )
);

dimensionedScalar rhoMin
(
    dimensionedScalar::lookupOrDefault
    (
        "rhoMin",
        pimple.dict(),
        dimDensity,
        0
    )
);

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));
multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

#include "createMRF.H"
#include "createClouds.H"
//#include "createRadiationModel.H"
//#include "createFvOptions.H"
