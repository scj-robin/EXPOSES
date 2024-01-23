data;
	infile "D:\RECHERCHE\METAGENOME\EXPOSES\FIGURES\MICA-AnovaMixed.txt" firstobs=2 expandtabs;
	input x$ Trt$ Indiv$ Response;
proc Mixed;
     class Trt Indiv;
     model Response = Trt;
	 random Indiv(Trt);
     estimate 'Trt1-Trt2' Trt 1 -1 0;
run;
