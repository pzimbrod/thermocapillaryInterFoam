{
    radiation->correct();
    rhoCp =  rho*fluid.Cp();
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- Debug Information
    Info<< "min/max(Cp) = "
	<< min(fluid.Cp()).value() << ", " << max(fluid.Cp()).value() << endl;

    Info<< "min/max(rhoCp) = "
	<< min(rhoCp).value() << ", " << max(rhoCp).value() << endl;
	
    Info<< "min/max(Q) = "
	<< min(radiation->Ru()).value() << ", " << max(radiation->Ru()).value() << endl;
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const surfaceScalarField rhoCpPhi(fvc::interpolate(fluid.Cp())*rhoPhi);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- Debug Information
    Info<< "min/max(rhoCpPhi) = "
	<< min(rhoCpPhi).value() << ", " << max(rhoCpPhi).value() << endl;
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const volScalarField kappaEff
    (
        "kappaEff",
        fluid.kappa() + fluid.Cp()*turbulence->mut()/fluid.Prt()
    );
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- Debug Information
    Info<< "min/max(kappa) = "
	<< min(fluid.kappa()).value() << ", " << max(fluid.kappa()).value() << endl;
	
    Info<< "min/max(mut) = "
	<< min(turbulence->mut()).value() << ", " << max(turbulence->mut()).value() << endl;
	
    Info<< "min/max(kappaEff) = "
	<< min(kappaEff).value() << ", " << max(kappaEff).value() << endl;
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoCp, T)
      + fvm::div(rhoCpPhi, T, "div(phi,T)")
      - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoCpPhi), T)
      - fvm::laplacian(kappaEff, T,  "laplacian(kappa,T)")
      ==
        fluid.heatTransfer(T)
      + radiation->ST(T)
      + fvOptions(rhoCp, T)
    );

    TEqn.relax();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //- Debug Information    
    Info<< "min/max(T) after relax() = "
        << min(T).value() << ", " << max(T).value() << endl;

    fvOptions.constrain(TEqn);
    
    Info<< "min/max(T) after constrain() = "
        << min(T).value() << ", " << max(T).value() << endl;

    TEqn.solve();
    
    Info<< "min/max(T) after solve() = "
        << min(T).value() << ", " << max(T).value() << endl;

    fvOptions.correct(T);
    
    Info<< "min/max(T) after fvOptions.correct() = "
        << min(T).value() << ", " << max(T).value() << endl;

    fluid.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
}
