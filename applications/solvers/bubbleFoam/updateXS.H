//Update XS

// TR[eI] (used as 1/(3*TR) in the flux diffusion coefficient) scales
// directly with rhok, so rhok must stay bounded away from zero -- an
// unfloored alpha2 can locally undercut phase2's own residualAlpha,
// driving TR toward zero and the diffusion coefficient toward infinity.
rhok=voidManual*max(alpha2, phase2.residualAlpha())*rho2/rhoReference;
rhod=(alpha1*rho1 + alpha2*rho2)/rhoReference*voidManual;
Tfuel=thermo2.T();

logT = log(max(0.5,(Tfuel/TrefXS)));

for (int eI=0; eI<energyGroups ; ++eI)
{
		SP[eI] = (SP0[eI] + alfaSP[eI] * logT)*rhok;
		NSF[eI] = (NSF0[eI] + alfaNSF[eI] * logT)*rhok;
		F[eI] = NSF[eI] * (1.0-precBetaTot) * chiPro[eI];
		M[eI] = (M0[eI] + alfaM[eI] * logT)*rhok;
		TR[eI] = (TR0[eI] + alfaTR[eI] * logT)*rhok;

    for (int eJ=0; eJ<energyGroups ; ++eJ)
		{
				int eIJ = energyGroups*eI + eJ;
				S[eI][eJ] = (S0[eIJ] + alfaS[eIJ] * logT)*rhok;
		}
}
