%let Rep = D:\RECHERCHE\RUPTURES\E-Budinska\FINAL RESULTS\analyses by EM\Chromosome 8\balanced\EM\Franck 0.5 outer;

data CHROMO8;
	infile "&Rep/rupt.CGH.chr.8.apres.txt" firstobs=2;
	input Group Patient Begin End;
proc Sort data=CHROMO8;
	by Group Patient Begin;
proc GChart data=CHROMO8;
	vbar Begin;
	by Group;
run;
