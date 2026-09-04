//Update Xenon inventory and cross sections

xenonSolute = rho2*phase2.Y("xenon")*alpha2;
xenonGas = rho1*phase1.Y("xenon")*alpha1;

xenonMass = xenonSolute + xenonGas;

for (int eI=0; eI<energyGroups; ++eI)
{
  XE[eI] = xenonMass*(N_Av/MW_xe) * alfaXE[eI];
}
