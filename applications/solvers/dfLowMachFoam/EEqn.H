{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      - dpdt
      - fvm::laplacian(turbulence->alphaEff(), he)
      + diffAlphaD
     ==
        fvc::div(hDiffCorrFlux)
    );

    EEqn.relax();

    EEqn.solve();
}
