// CSF surface-tension force (multiphaseInterSystem::surfaceTensionForce), with
// optional curvature from alpha smoothed nSmooth times by cell-face-cell
// averaging (Lafaurie et al. 1994); nSmooth 0 is the stock force.
// The force direction still uses the sharp snGrad(alpha).

Foam::tmp<Foam::surfaceScalarField> surfaceTensionForce
(
    const Foam::multiphaseInter::multiphaseSystem& fluid,
    const Foam::label nSmooth
)
{
    using namespace Foam;
    using namespace Foam::multiphaseInter;

    if (nSmooth <= 0)
    {
        return fluid.surfaceTensionForce();
    }

    auto smooth = [nSmooth](const volScalarField& alpha)
    {
        tmp<volScalarField> ta(alpha);
        for (label i = 0; i < nSmooth; ++i)
        {
            ta = fvc::average(fvc::interpolate(ta()));
        }
        return ta;
    };

    auto tstf = surfaceScalarField::New
    (
        "surfaceTensionForce",
        IOobject::NO_REGISTER,
        fluid.mesh(),
        dimensionedScalar({1, -2, -2, 0, 0, 0}, Zero)
    );
    auto& stf = tstf.ref();
    stf.setOriented();

    const auto& phases = fluid.phases();

    forAllConstIters(phases, iter1)
    {
        const phaseModel& phase1 = iter1.val();
        const volScalarField& alpha1 = phase1;

        auto iter2 = iter1;
        for (++iter2; iter2 != phases.cend(); ++iter2)
        {
            const phaseModel& phase2 = iter2.val();
            const volScalarField& alpha2 = phase2;

            stf +=
                fvc::interpolate
                (
                    fluid.surfaceTensionCoeff
                    (
                        phasePairKey(phase1.name(), phase2.name())
                    )
                )
              * fvc::interpolate(fluid.K(smooth(alpha1)(), smooth(alpha2)()))
              * (
                    fvc::interpolate(alpha2)*fvc::snGrad(alpha1)
                  - fvc::interpolate(alpha1)*fvc::snGrad(alpha2)
                );
        }
    }

    return tstf;
}
