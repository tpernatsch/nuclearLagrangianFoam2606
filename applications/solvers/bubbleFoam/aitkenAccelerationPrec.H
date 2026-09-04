if (aitkenIterCounter==0) // only happening at timeStepIteration==1
{
	Info << "Aitken step 0 " << endl;
	integralFlux1_0=(fvc::domainIntegrate(flux1)).value();
	integralFlux2_0=(fvc::domainIntegrate(flux2)).value();
	integralFlux3_0=(fvc::domainIntegrate(flux3)).value();
	integralFlux4_0=(fvc::domainIntegrate(flux4)).value();
	integralFlux5_0=(fvc::domainIntegrate(flux5)).value();
	integralFlux6_0=(fvc::domainIntegrate(flux6)).value();

	integralPrec1_0=(fvc::domainIntegrate(prec1)).value();
	integralPrec2_0=(fvc::domainIntegrate(prec2)).value();
	integralPrec3_0=(fvc::domainIntegrate(prec3)).value();
	integralPrec4_0=(fvc::domainIntegrate(prec4)).value();
	integralPrec5_0=(fvc::domainIntegrate(prec5)).value();
	integralPrec6_0=(fvc::domainIntegrate(prec6)).value();
	integralPrec7_0=(fvc::domainIntegrate(prec7)).value();
	integralPrec8_0=(fvc::domainIntegrate(prec8)).value();

	aitkenIterCounter++ ;
}
else
{
	if (aitkenIterCounter==1)
	{
		Info << "Aitken step 1 " << endl;
		integralFlux1_1=(fvc::domainIntegrate(flux1)).value();
		integralFlux2_1=(fvc::domainIntegrate(flux2)).value();
		integralFlux3_1=(fvc::domainIntegrate(flux3)).value();
		integralFlux4_1=(fvc::domainIntegrate(flux4)).value();
		integralFlux5_1=(fvc::domainIntegrate(flux5)).value();
		integralFlux6_1=(fvc::domainIntegrate(flux6)).value();

		integralPrec1_1=(fvc::domainIntegrate(prec1)).value();
		integralPrec2_1=(fvc::domainIntegrate(prec2)).value();
		integralPrec3_1=(fvc::domainIntegrate(prec3)).value();
		integralPrec4_1=(fvc::domainIntegrate(prec4)).value();
		integralPrec5_1=(fvc::domainIntegrate(prec5)).value();
		integralPrec6_1=(fvc::domainIntegrate(prec6)).value();
		integralPrec7_1=(fvc::domainIntegrate(prec7)).value();
		integralPrec8_1=(fvc::domainIntegrate(prec8)).value();

		aitkenIterCounter++;
	}
	else
	{
		if (aitkenIterCounter==2)
		{
			integralFlux1_2=(fvc::domainIntegrate(flux1)).value();
			integralFlux2_2=(fvc::domainIntegrate(flux2)).value();
			integralFlux3_2=(fvc::domainIntegrate(flux3)).value();
			integralFlux4_2=(fvc::domainIntegrate(flux4)).value();
			integralFlux5_2=(fvc::domainIntegrate(flux5)).value();
			integralFlux6_2=(fvc::domainIntegrate(flux6)).value();

			integralPrec1_2=(fvc::domainIntegrate(prec1)).value();
			integralPrec2_2=(fvc::domainIntegrate(prec2)).value();
			integralPrec3_2=(fvc::domainIntegrate(prec3)).value();
			integralPrec4_2=(fvc::domainIntegrate(prec4)).value();
			integralPrec5_2=(fvc::domainIntegrate(prec5)).value();
			integralPrec6_2=(fvc::domainIntegrate(prec6)).value();
			integralPrec7_2=(fvc::domainIntegrate(prec7)).value();
			integralPrec8_2=(fvc::domainIntegrate(prec8)).value();

			bool aitkenYesNo=true;

			if
			(
				mag(integralFlux1_2-integralFlux1_1)/(max(mag(integralFlux1_1-integralFlux1_0),1e-20)) >0.975 ||
				mag(integralFlux2_2-integralFlux2_1)/(max(mag(integralFlux2_1-integralFlux2_0),1e-20)) >0.975 ||
				mag(integralFlux3_2-integralFlux3_1)/(max(mag(integralFlux3_1-integralFlux3_0),1e-20)) >0.975 ||
				mag(integralFlux4_2-integralFlux4_1)/(max(mag(integralFlux4_1-integralFlux4_0),1e-20)) >0.975 ||
				mag(integralFlux5_2-integralFlux5_1)/(max(mag(integralFlux5_1-integralFlux5_0),1e-20)) >0.975 ||
				mag(integralFlux6_2-integralFlux6_1)/(max(mag(integralFlux6_1-integralFlux6_0),1e-20)) >0.975 ||
				mag(integralPrec1_2-integralPrec1_1)/(max(mag(integralPrec1_1-integralPrec1_0),1e-20)) >0.975 ||
				mag(integralPrec2_2-integralPrec2_1)/(max(mag(integralPrec2_1-integralPrec2_0),1e-20)) >0.975 ||
				mag(integralPrec3_2-integralPrec3_1)/(max(mag(integralPrec3_1-integralPrec3_0),1e-20)) >0.975 ||
				mag(integralPrec4_2-integralPrec4_1)/(max(mag(integralPrec4_1-integralPrec4_0),1e-20)) >0.975 ||
				mag(integralPrec5_2-integralPrec5_1)/(max(mag(integralPrec5_1-integralPrec5_0),1e-20)) >0.975 ||
				mag(integralPrec6_2-integralPrec6_1)/(max(mag(integralPrec6_1-integralPrec6_0),1e-20)) >0.975 ||
				mag(integralPrec7_2-integralPrec7_1)/(max(mag(integralPrec7_1-integralPrec7_0),1e-20)) >0.975 ||
				mag(integralPrec8_2-integralPrec8_1)/(max(mag(integralPrec8_1-integralPrec8_0),1e-20)) >0.975 ||
				( integralFlux1_2-2*integralFlux1_1+integralFlux1_0 == 0 ) ||
				( integralFlux2_2-2*integralFlux2_1+integralFlux2_0 == 0 ) ||
				( integralFlux3_2-2*integralFlux3_1+integralFlux3_0 == 0 ) ||
				( integralFlux4_2-2*integralFlux4_1+integralFlux4_0 == 0 ) ||
				( integralFlux5_2-2*integralFlux5_1+integralFlux5_0 == 0 ) ||
				( integralFlux6_2-2*integralFlux6_1+integralFlux6_0 == 0 ) ||
				( integralPrec1_2-2*integralPrec1_1+integralPrec1_0 == 0 ) ||
				( integralPrec2_2-2*integralPrec2_1+integralPrec2_0 == 0 ) ||
				( integralPrec3_2-2*integralPrec3_1+integralPrec3_0 == 0 ) ||
				( integralPrec4_2-2*integralPrec4_1+integralPrec4_0 == 0 ) ||
				( integralPrec5_2-2*integralPrec5_1+integralPrec5_0 == 0 ) ||
				( integralPrec6_2-2*integralPrec6_1+integralPrec6_0 == 0 ) ||
				( integralPrec7_2-2*integralPrec7_1+integralPrec7_0 == 0 ) ||
				( integralPrec8_2-2*integralPrec8_1+integralPrec8_0 == 0 )
			)
			{
				aitkenYesNo=false;
			}
			else
			{
				aitkenYesNo=true;
			}

			if (aitkenYesNo==true)
			{
				Info << "Aitken step 2 " << endl;

				scalar aitkenScaling1 =	((integralFlux1_0 - (integralFlux1_1-integralFlux1_0)*(integralFlux1_1-integralFlux1_0)/(integralFlux1_2-2*integralFlux1_1+integralFlux1_0))/integralFlux1_2);
				scalar aitkenScaling2 = ((integralFlux2_0 - (integralFlux2_1-integralFlux2_0)*(integralFlux2_1-integralFlux2_0)/(integralFlux2_2-2*integralFlux2_1+integralFlux2_0))/integralFlux2_2);
				scalar aitkenScaling3 =	((integralFlux3_0 - (integralFlux3_1-integralFlux3_0)*(integralFlux3_1-integralFlux3_0)/(integralFlux3_2-2*integralFlux3_1+integralFlux3_0))/integralFlux3_2);
				scalar aitkenScaling4 =	((integralFlux4_0 - (integralFlux4_1-integralFlux4_0)*(integralFlux4_1-integralFlux4_0)/(integralFlux4_2-2*integralFlux4_1+integralFlux4_0))/integralFlux4_2);
				scalar aitkenScaling5 =	((integralFlux5_0 - (integralFlux5_1-integralFlux5_0)*(integralFlux5_1-integralFlux5_0)/(integralFlux5_2-2*integralFlux5_1+integralFlux5_0))/integralFlux5_2);
				scalar aitkenScaling6 =	((integralFlux6_0 - (integralFlux6_1-integralFlux6_0)*(integralFlux6_1-integralFlux6_0)/(integralFlux6_2-2*integralFlux6_1+integralFlux6_0))/integralFlux6_2);

				scalar aitkenScalingPrec1 = ((integralPrec1_0 - (integralPrec1_1-integralPrec1_0)*(integralPrec1_1-integralPrec1_0)/(integralPrec1_2-2*integralPrec1_1+integralPrec1_0))/integralPrec1_2);
				scalar aitkenScalingPrec2 =	((integralPrec2_0 - (integralPrec2_1-integralPrec2_0)*(integralPrec2_1-integralPrec2_0)/(integralPrec2_2-2*integralPrec2_1+integralPrec2_0))/integralPrec2_2);
				scalar aitkenScalingPrec3 =	((integralPrec3_0 - (integralPrec3_1-integralPrec3_0)*(integralPrec3_1-integralPrec3_0)/(integralPrec3_2-2*integralPrec3_1+integralPrec3_0))/integralPrec3_2);
				scalar aitkenScalingPrec4 =	((integralPrec4_0 - (integralPrec4_1-integralPrec4_0)*(integralPrec4_1-integralPrec4_0)/(integralPrec4_2-2*integralPrec4_1+integralPrec4_0))/integralPrec4_2);
				scalar aitkenScalingPrec5 =	((integralPrec5_0 - (integralPrec5_1-integralPrec5_0)*(integralPrec5_1-integralPrec5_0)/(integralPrec5_2-2*integralPrec5_1+integralPrec5_0))/integralPrec5_2);
				scalar aitkenScalingPrec6 =	((integralPrec6_0 - (integralPrec6_1-integralPrec6_0)*(integralPrec6_1-integralPrec6_0)/(integralPrec6_2-2*integralPrec6_1+integralPrec6_0))/integralPrec6_2);
				scalar aitkenScalingPrec7 =	((integralPrec7_0 - (integralPrec7_1-integralPrec7_0)*(integralPrec7_1-integralPrec7_0)/(integralPrec7_2-2*integralPrec7_1+integralPrec7_0))/integralPrec7_2);
				scalar aitkenScalingPrec8 =	((integralPrec8_0 - (integralPrec8_1-integralPrec8_0)*(integralPrec8_1-integralPrec8_0)/(integralPrec8_2-2*integralPrec8_1+integralPrec8_0))/integralPrec8_2);

				flux1 = flux1*aitkenScaling1;
				flux2 = flux2*aitkenScaling2;
				flux3 = flux3*aitkenScaling3;
				flux4 = flux4*aitkenScaling4;
				flux5 = flux5*aitkenScaling5;
				flux6 = flux6*aitkenScaling6;

				prec1 = prec1*aitkenScalingPrec1;
				prec2 = prec2*aitkenScalingPrec2;
				prec3 = prec3*aitkenScalingPrec3;
				prec4 = prec4*aitkenScalingPrec4;
				prec5 = prec5*aitkenScalingPrec5;
				prec6 = prec6*aitkenScalingPrec6;
				prec7 = prec7*aitkenScalingPrec7;
				prec8 = prec8*aitkenScalingPrec8;

				/*
				integralFlux1_0 = (fvc::domainIntegrate(flux1)).value();
				integralFlux2_0 = (fvc::domainIntegrate(flux2)).value();
				integralFlux3_0 = (fvc::domainIntegrate(flux3)).value();
				integralFlux4_0 = (fvc::domainIntegrate(flux4)).value();
				integralFlux5_0 = (fvc::domainIntegrate(flux5)).value();
				integralFlux6_0 = (fvc::domainIntegrate(flux6)).value();
				*/

				promptHeatSource = (1 - decBetaTot)*(flux1*SP1 + flux2*SP2 + flux3*SP3 + flux4*SP4 + flux5*SP5 + flux6*SP6);  //dipende da come è calcolato Sigma Power
				NSFflu = NSF1*flux1 + NSF2*flux2 + NSF3*flux3 + NSF4*flux4 + NSF5*flux5 + NSF6*flux6; //nu*sigmaf*flux serve per precursori
				SPflu = SP1*flux1 + SP2*flux2 + SP3*flux3 + SP4*flux4 + SP5*flux5 + SP6*flux6;  //Sigma power * flux serve per il decay

                //neutroni da scattering
				scaNeu1=S21*flux2 + S31*flux3 + S41*flux4 + S51*flux5 + S61*flux6;
				scaNeu2=S12*flux1 + S32*flux3 + S42*flux4 + S52*flux5 + S62*flux6;
				scaNeu3=S13*flux1 + S23*flux2 + S43*flux4 + S53*flux5 + S63*flux6;
				scaNeu4=S14*flux1 + S24*flux2 + S34*flux3 + S54*flux5 + S64*flux6;
				scaNeu5=S15*flux1 + S25*flux2 + S35*flux3 + S45*flux4 + S65*flux6;
				scaNeu6=S16*flux1 + S26*flux2 + S36*flux3 + S46*flux4 + S56*flux5;

				//neutroni da fissioni avvenute in altri gruppi
				othGroFisNeu1=(1-betaTot)*(NSF2*flux2 + NSF3*flux3 + NSF4*flux4 + NSF5*flux5 + NSF6*flux6)*chiPro1;
				othGroFisNeu2=(1-betaTot)*(NSF1*flux1 + NSF3*flux3 + NSF4*flux4 + NSF5*flux5 + NSF6*flux6)*chiPro2;
				othGroFisNeu3=(1-betaTot)*(NSF1*flux1 + NSF2*flux2 + NSF4*flux4 + NSF5*flux5 + NSF6*flux6)*chiPro3;
				othGroFisNeu4=(1-betaTot)*(NSF1*flux1 + NSF2*flux2 + NSF3*flux3 + NSF5*flux5 + NSF6*flux6)*chiPro4;
				othGroFisNeu5=(1-betaTot)*(NSF1*flux1 + NSF2*flux2 + NSF3*flux3 + NSF4*flux4 + NSF6*flux6)*chiPro5;
				othGroFisNeu6=(1-betaTot)*(NSF1*flux1 + NSF2*flux2 + NSF3*flux3 + NSF4*flux4 + NSF5*flux5)*chiPro6;

				//flusso totale
				flux=flux1+flux2+flux3+flux4+flux5+flux6;

				delNeuSource1 = (prec1*lam1 + prec2*lam2 + prec3*lam3 + prec4*lam4 + prec5*lam5 + prec6*lam6 + prec7*lam7 + prec8*lam8)*chiDel1;
				delNeuSource2 = (prec1*lam1 + prec2*lam2 + prec3*lam3 + prec4*lam4 + prec5*lam5 + prec6*lam6 + prec7*lam7 + prec8*lam8)*chiDel2;
				delNeuSource3 = (prec1*lam1 + prec2*lam2 + prec3*lam3 + prec4*lam4 + prec5*lam5 + prec6*lam6 + prec7*lam7 + prec8*lam8)*chiDel3;
				delNeuSource4 = (prec1*lam1 + prec2*lam2 + prec3*lam3 + prec4*lam4 + prec5*lam5 + prec6*lam6 + prec7*lam7 + prec8*lam8)*chiDel4;
				delNeuSource5 = (prec1*lam1 + prec2*lam2 + prec3*lam3 + prec4*lam4 + prec5*lam5 + prec6*lam6 + prec7*lam7 + prec8*lam8)*chiDel5;
				delNeuSource6 = (prec1*lam1 + prec2*lam2 + prec3*lam3 + prec4*lam4 + prec5*lam5 + prec6*lam6 + prec7*lam7 + prec8*lam8)*chiDel6;

				//	aitkenIterNo++ ;
				aitkenIterCounter = 0;

				Info<< "Done "  << endl;
			}
			else
			{
				Info << "Aitken not performed because of diverging series" << endl;

				integralFlux1_0 = integralFlux1_1;
				integralFlux2_0 = integralFlux2_1;
				integralFlux3_0 = integralFlux3_1;
				integralFlux4_0 = integralFlux4_1;
				integralFlux5_0 = integralFlux5_1;
				integralFlux6_0 = integralFlux6_1;

				integralPrec1_0 = integralPrec1_1;
				integralPrec2_0 = integralPrec2_1;
				integralPrec3_0 = integralPrec3_1;
				integralPrec4_0 = integralPrec4_1;
				integralPrec5_0 = integralPrec5_1;
				integralPrec6_0 = integralPrec6_1;
				integralPrec7_0 = integralPrec7_1;
				integralPrec8_0 = integralPrec8_1;

				integralFlux1_1 = integralFlux1_2;
				integralFlux2_1 = integralFlux2_2;
				integralFlux3_1 = integralFlux3_2;
				integralFlux4_1 = integralFlux4_2;
				integralFlux5_1 = integralFlux5_2;
				integralFlux6_1 = integralFlux6_2;

				integralPrec1_1 = integralPrec1_2;
				integralPrec2_1 = integralPrec2_2;
				integralPrec3_1 = integralPrec3_2;
				integralPrec4_1 = integralPrec4_2;
				integralPrec5_1 = integralPrec5_2;
				integralPrec6_1 = integralPrec6_2;
				integralPrec7_1 = integralPrec7_2;
				integralPrec8_1 = integralPrec8_2;

				aitkenIterCounter--;
			}
		}
	}
}
