goptions reset;

data PVAL;
	n = 5000;
	EP1 = 0.05;
	a = 0.2;
	do i=1 to n;
		H = 0;
		P = ranuni(2);
		if i < a*n then do;
			H = 1;
			P = -log(P)*EP1;
		end;
		output;
	end;
run;
*ods epsc file="D:/RECHERCHE/EXPRESSION/EXPOSES/Figures/HistoPvalRef.eps";
proc GChart data=PVAL;
	pattern1 c=green;
	pattern2 c=red;
	vbar P / midpoints=(0.025 to 0.975 by 0.05) 
		subgroup=H nospace href=200 lref=2 woutline=3;
*ods epsc close;
run; *quit;


data PVAL;
	n = 5000;
	EP1 = 0.01;
	a = 0.2;
	do i=1 to n;
		H = 0;
		P = ranuni(2);
		if i < a*n then do;
			H = 1;
			P = -log(P)*EP1;
		end;
		F = probit(P);
		output;
	end;
run;
*ods epsc file="D:/RECHERCHE/EXPRESSION/EXPOSES/Figures/HistoPval.eps";
proc GChart data=PVAL;
	pattern1 c=green;
	pattern2 c=red;
	vbar P / midpoints=(0.025 to 0.975 by 0.05) 
		subgroup=H nospace woutline=3;
*ods epsc close;
run; *quit;

*ods epsc file="D:/RECHERCHE/EXPRESSION/EXPOSES/Figures/HistoProbit.eps";
proc GChart data=PVAL;
	pattern1 c=green;
	pattern2 c=red;
	vbar P / midpoints=(0.025 to 0.975 by 0.05) 
		subgroup=H nospace woutline=3;
*ods epsc close;
run; *quit;
*ods epsc file="D:/RECHERCHE/EXPRESSION/EXPOSES/Figures/HistoPvalRef.eps";
proc GChart data=PVAL;
	pattern1 c=green;
	pattern2 c=red;
	vbar F / /*midpoints=(0.025 to 0.975 by 0.05) */
		subgroup=H nospace woutline=3;
*ods epsc close;
run; *quit;
