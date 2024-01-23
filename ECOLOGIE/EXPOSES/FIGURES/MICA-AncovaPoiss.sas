data;
	infile "D:\RECHERCHE\METAGENOME\EXPOSES\FIGURES\MICA-AncovaPois.txt" firstobs=2 expandtabs;
	input x$ Strain$ Prep$ Time$ Response;
proc GenMod;
     class Strain Time Prep;
     model Response = Strain Time Strain*Time Prep / type1 dist=Poisson link=log;
     estimate 'TimeEnd' Strain 1 -1 Strain*Time 0 0 0 0 1 0 0 0 0 -1;
     estimate 'TakeOff' Strain*Time 1 -2 1 0 0 -1 2 -1 0 0;
run;
