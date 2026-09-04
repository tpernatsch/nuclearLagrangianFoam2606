volScalarField Uy=U2.component(1);
volScalarField phiPump=(pump.component(1))*Uy;
//la geometria EVOL 2D ha y come coordinata verticale, il 3D ha z
//per passare da uno all'altro cambiare component(1) in component(2) e Uy in Uz

scalar flowRate = (fvc::domainIntegrate(phiPump).value()) / pumpHeigth;
Info << nl << "Volumetric flow rate = " << flowRate << endl;
