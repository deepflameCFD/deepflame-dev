{
    volScalarField& he = thermo.he();
    if (constProp == "volume") he[0] = u0 + p[0]/rho[0];
    chemistry.correctThermo();
}
