proc Factex;
	factors A B C D E;
	size design = 8;
	model est=(A B C D E);
	examine aliasing confounding design;
run;
