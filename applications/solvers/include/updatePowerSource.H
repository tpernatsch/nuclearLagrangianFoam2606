{
    // Reset fields
    qPrompt *= 0.0;
    qDecay *= 0.0;

    // Evaluate prompt heat source from current fluxes and cross-sections
    forAll(flux, i)
    {
        qPrompt += (1 - decBetaTot)*Ef[i]*Sf[i]*flux[i];
    }

    // Evaluate delayed heat source from current precursors
    forAll(dec, i)
    {
        qDecay += decLambda[i]*dec[i];
    }

    q = qPrompt + qDecay;
    
    Qprompt = unitField*fvc::domainIntegrate(qPrompt).value()*dimensionedScalar("", dimPower, 1);

    Qdecay = unitField*fvc::domainIntegrate(qDecay).value()*dimensionedScalar("", dimPower, 1);

    Q = Qprompt + Qdecay;
}
