{
    forAll(flux, i)
    {
        // Temp zero fields
        volScalarField tSNS = zeroField*dimensionedScalar("", dimless/dimVolume/dimTime, 0);
        volScalarField tFNS = zeroField*dimensionedScalar("", dimless/dimVolume/dimTime, 0);
        volScalarField tDNS = zeroField*dimensionedScalar("", dimless/dimVolume/dimTime, 0);

        // Evaluate neutron sources based on current fluxes and cross sections
        forAll(flux, j)
        {
            if (j!=i)
            {
                tSNS += Ss[nEnergyGroups*j + i]*flux[j];
                tFNS += chiPrompt[i]*(1 - precBetaTot)*Nu[j]*Sf[j]*flux[j];
            }
        }

        forAll(prec, j)
        {
            tDNS += chiDelayed[i]*precLambda[j]*prec[j];
        }

        // Repopulate pointer lists
        neutronSource[i] = tSNS + tFNS + tDNS;
    }
}
