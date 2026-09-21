MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
  ==
  //parcels.SU(U, true)
    FP1.SU(U, true)
  //+ FP2.SU(U, true)
  + fvOptions(U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p_rgh));
    fvOptions.correct(U);
}
