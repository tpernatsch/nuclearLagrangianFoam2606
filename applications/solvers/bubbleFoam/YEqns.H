{
    autoPtr<phaseSystem::massTransferTable>
        massTransferPtr(fluid.massTransfer());

    phaseSystem::massTransferTable&
        massTransfer(massTransferPtr());

    if (!phase1.pure())
    {
        UPtrList<volScalarField>& Y1 = phase1.YActiveRef();

        forAll(Y1, i)
        {
          if (Y1[i].name() == "xenon.air")
          {
            fvScalarMatrix Y1iEqn
            (
                phase1.YiEqn(Y1[i])
             ==
               *massTransfer[Y1[i].name()]
              + fvOptions(alpha1, rho1, Y1[i])
              - alpha1*rho1*fvm::Sp(lambda_xe,Y1[i])
            );

            Y1iEqn.relax();
            Y1iEqn.solve(mesh.solver("Yi"));
          }
          else
          {
            fvScalarMatrix Y1iEqn
            (
                phase1.YiEqn(Y1[i])
             ==
               *massTransfer[Y1[i].name()]
              + fvOptions(alpha1, rho1, Y1[i])
            );

            Y1iEqn.relax();
            Y1iEqn.solve(mesh.solver("Yi"));
          }
        }
    }

    if (!phase2.pure())
    {
        UPtrList<volScalarField>& Y2 = phase2.YActiveRef();

        forAll(Y2, i)
        {
          if (Y2[i].name() == "xenon.water")
          {
            fvScalarMatrix Y2iEqn
            (
                phase2.YiEqn(Y2[i])
             ==
               *massTransfer[Y2[i].name()]
              + fvOptions(alpha2, rho2, Y2[i])
              + fissionRate*y_xe*MW_xe/N_Av
              - alpha2*rho2*fvm::Sp(lambda_xe,Y2[i])
              - alpha2*rho2*fvm::Sp(xeCapture,Y2[i])
            );

            Y2iEqn.relax();
            Y2iEqn.solve(mesh.solver("Yi"));
          }
          else
          {
            fvScalarMatrix Y2iEqn
            (
                phase2.YiEqn(Y2[i])
             ==
               *massTransfer[Y2[i].name()]
              + fvOptions(alpha2, rho2, Y2[i])
            );

            Y2iEqn.relax();
            Y2iEqn.solve(mesh.solver("Yi"));
          }
        }
    }
}
