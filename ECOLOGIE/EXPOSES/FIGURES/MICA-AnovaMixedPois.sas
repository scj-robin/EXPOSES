data;
	infile "D:\RECHERCHE\METAGENOME\EXPOSES\FIGURES\MICA-AnovaMixedPois.txt" firstobs=2 expandtabs;
	input x$ Trt$ Indiv$ Response;
proc glimmix;
	class Trt Indiv;
	model Response = Trt / dist=Poisson link=log;
	random Indiv(Trt);
	estimate 'Trt1-Trt2' Trt 1 -1 0;
run;
