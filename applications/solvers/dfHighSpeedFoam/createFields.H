#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;

psiThermo* pThermo = new hePsiThermo<psiThermo, CanteraMixture>(mesh, word::null);
psiThermo& thermo = *pThermo;

//move from creatFieldRefs.H to createFields.H
//p needed to be created before e
volScalarField& p = thermo.p();
const volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();
const volScalarField& mu = thermo.mu();

bool inviscid(true);
if (max(mu.primitiveField()) > 0.0)
{
    inviscid = false;
}

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

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volScalarField e = thermo.he() - p/rho;
volScalarField& ha = thermo.he();

volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);

volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(e + 0.5*magSqr(U))
);

surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dimless, -1.0)
);

surfaceScalarField phi("phi", fvc::flux(rhoU));

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

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

dfChemistryModel<basicThermo> chemistry(thermo);
PtrList<volScalarField>& Y = chemistry.Y();
const word inertSpecie(chemistry.lookup("inertSpecie"));
const label inertIndex(chemistry.species()[inertSpecie]);

chemistry.correctThermo();
Info<< "At initial time, min/max(T) = " << min(T).value() << ", " << max(T).value() << endl;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());