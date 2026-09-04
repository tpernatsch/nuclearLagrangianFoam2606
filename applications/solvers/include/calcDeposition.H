// Steady-state fission products deposition

dimensionedScalar fpiTot("fpiTot", fvc::domainIntegrate(fpi));  // Total number of fission products atoms
//dimensionedScalar fpsiTot("fpsiTot", fvc::domainIntegrate(fpsi));  // Total number of dissolved fission products atoms
//dimensionedScalar fppiTot("fppiTot", fvc::domainIntegrate(fppi));  // Total number of precipited fission products atoms
dimensionedScalar fpiSrcTot("fpiSourceTot", fvc::domainIntegrate(FPYield[i]*Rf));  // Global fission products production rate
dimensionedScalar fpiDepTot("fpiDepTot", fpiTot.dimensions()/dimTime, 0.0);  // Global fission products deposition rate
dimensionedScalar fpiDecTot("fpiDecTot", fpiTot*FPLambda[i]);  // Global fission products decay rate

dimensionedScalar fpiDerTot("fpiDerTot", fvc::domainIntegrate(fvc::ddt(fpi)));  // Global fission products concentration derivative
dimensionedScalar fpiAvg("fpiAvg", fpiTot/gSum(mesh.V()));  // Average fission products concentration

// New variables name for each element of the arrays
const volScalarField Lambda = (Dfpi/Dfpi)*FPLambda[i];  // FPLambda from dimensionedScalar to volScalarField
const volScalarField PD = (Dfpi/Dfpi)*FPPd[i];  // FPPd from dimensionedScalar to volScalarField
const volScalarField FpiAvg = (Dfpi/Dfpi)*fpiAvg;
volScalarField& di = d[i];
volScalarField& qdecay_walli = qdecay_wall[i];
volScalarField& gammai = gamma[i];

// Assign names to boundary fields, constant properties or assume constant fields from previous steps
const volScalarField::Boundary& LambdaBf = Lambda.boundaryField();
const volScalarField::Boundary& PDBf = PD.boundaryField();
const volScalarField::Boundary& fpiBf = fpi.boundaryField();
const volScalarField::Boundary& DfpBf = Dfpi.boundaryField();
const volScalarField::Boundary& FpiAvgBf = FpiAvg.boundaryField();

// Assign names to boundary fields, to be computed
volScalarField::Boundary& fpdiBf = fpdi.boundaryFieldRef();
volScalarField::Boundary& diBf = di.boundaryFieldRef();
volScalarField::Boundary& qdecay_walliBf = qdecay_walli.boundaryFieldRef();
volScalarField::Boundary& gammaiBf = gammai.boundaryFieldRef();

// For all points in the prescribed patches
forAll(fpiBf, patchi)
{
    if (fpiBf[patchi].patch().type() == "wall")
    {
        fpdiBf[patchi] = - DfpBf[patchi]*fpiBf[patchi].snGrad();  // Evaluate fission products deposition rate
	    fpiDepTot.value() += gSum(fpdiBf[patchi]*fpiBf[patchi].patch().magSf());  // Evaluate global fission products deposition rate
        gammaiBf[patchi] = fpdiBf[patchi]/FpiAvgBf[patchi];
        diBf[patchi] = fpdiBf[patchi]/LambdaBf[patchi];  // Evaluation of radionuclides deposited
        qdecay_walliBf[patchi] = diBf[patchi] * PDBf[patchi] * 1.6022e-16;  // Evaluate decay heat flux
    }
}

scalar balanceError(1.0 - ((fpiDepTot + fpiDecTot + fpiDerTot)/fpiSrcTot).value());  // Evaluate balance error for Van't Hoff equilibrium equation
//scalar equilibriumError(max(mag(1.0 - Ki*fpi/(fpsi + dimensionedScalar("VSMALL", fpsi.dimensions(), VSMALL)))).value());

Info<< "Specie balance: " //<< fpsiTot.value() << " at. (sol.), "
                          //<< fppiTot.value() << " at. (prec.), "
                          << fpiSrcTot.value() << " at./s (source), "
                          << fpiDepTot.value() << " at./s (dep.), "
                          << fpiDecTot.value() << " at./s (decay), "
                          << fpiDerTot.value() << " at./s (acc.), "
                          //<< -1 << " equilibrium rel. error, "
                          << balanceError << " balance rel. error" << endl;
