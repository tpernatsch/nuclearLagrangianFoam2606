{
    // Evaluate density and temperature fields
    rhok = 1.0 - beta*(T - TRef);
    rhoNorm = rhok/(1.0 - beta*(TRefXS - TRef));
    logT = log(max(0.5, T/TRefXS));

    // Evaluate temperature- and density-dependent nuclear properties
    forAll(flux, i)
    {
        D[i] = (D_ref[i] + DopplerControl*alphaD_ref[i]*logT)/rhoNorm;
        Sa[i] = (Sa_ref[i] + DopplerControl*alphaSa_ref[i]*logT)*rhoNorm;
        Sf[i] = (1/Keff)*(Sf_ref[i] + DopplerControl*alphaSf_ref[i]*logT)*rhoNorm;
    }

    forAll(Ss, i)
    {
        Ss[i] = (Ss_ref[i] + DopplerControl*alphaSs_ref[i]*logT)*rhoNorm;
    }

    // Evaluate net neutron production operator
    forAll(flux, i)
    {
        // Total removal cross-section of i-th group due to scattering to other groups
        volScalarField Ssi = 0*Ss[0];

        forAll(flux, j)
        {
            if (j!=i)
            {
                Ssi += Ss[nEnergyGroups*i + j];
            }
        }

        // Neutrons produced from fissions in i-th group minus neutrons removed
        M[i] = chiPrompt[i]*(1-precBetaTot)*Nu[i]*Sf[i] - (Sa[i] + Ssi);
    }
}
