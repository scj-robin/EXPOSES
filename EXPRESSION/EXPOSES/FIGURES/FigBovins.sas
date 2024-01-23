goptions reset;

data EXPRESS;
	infile 'D:\RECHERCHE\EXPRESSION\ETUDES\BOVINS\ModeleC-Hennequet\hybrid.txt' 
		firstobs=2 expandtabs;
	input ampli$ Tissue$ numampli nummemb SigMean gene$ 
		LogExp taille contig gc nature$ nonpris$ membrane$ 
		tube hybridation four typememb$;
	if NumMemb = 1;
	if Ampli = 'A+' then delete;
	keep Gene Ampli Tissue NumAmpli NumMemb LogExp;
proc Sort data=EXPRESS;
	by Gene Ampli Tissue NumAmpli NumMemb;
proc Freq data=EXPRESS noprint;
	tables Gene / out=FREQ;
run;

proc Anova data=EXPRESS noprint outstat=ANOVA;
	class Gene Ampli Tissue;
	model LogExp = Tissue Ampli Tissue*Ampli;
	by Gene;
run;
quit;
proc Sort data=ANOVA out=TMP;
	by _source_;
proc GChart data=TMP;
	vbar F;
	by _source_;
proc Transpose data=ANOVA out=SELECT;
	by Gene;
	var Prob;
	id _source_;
proc GPlot data=SELECT;
	plot Tissue * Ampli_Tissue;
data SELECT;
	set SELECT;
	if Ampli + Tissue + Ampli_Tissue < 0.5;
proc Sort data=SELECT;
	by Tissue;
run;

data EXEMPLE;
	set EXPRESS;
	if Gene = '14a04';
proc Anova data=EXEMPLE;
	class Ampli Tissue;
	model LogExp = Tissue Ampli Tissue*Ampli;
	means Tissue Ampli Ampli*Tissue;
run;

