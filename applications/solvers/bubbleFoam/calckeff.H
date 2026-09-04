FS_integral_old = fvc::domainIntegrate(NSFflu_old).value();

FS_integral = fvc::domainIntegrate(NSFflu).value();

Info  << "source_old = " <<  FS_integral_old << endl;

Info  << "source = " <<  FS_integral << endl;

Keff_new=Keff*(FS_integral/FS_integral_old);

res=(Keff_new-Keff)/Keff;

reactivity = (Keff_new-1)/Keff_new*1e5;

//Keff_new=0.1*Keff_new+0.9*Keff;

Info  << "Keff = " << Keff << endl;

Info  << "Keff_new = " << Keff_new << endl;

Info  << "res = " << res << endl;

Info << "reactivity = " << reactivity << " pcm" << endl;

Keff=Keff_new;
