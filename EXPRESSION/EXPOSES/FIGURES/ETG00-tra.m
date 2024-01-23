clear  all, clc 
LOADPATH = "../../occurrrences/PgmOccur:";
EDITOR = 'emacs';
FORMAT = 'COMPACT';

disp('ESTIMATION DU RAPPORT f0 / f DANS ETGC00');

m = 50
n1m = 1*m, mu1m = -2, sigma1m = 1
n0 = 8*m,  mu0 = 0,   sigma0 = 1
n1p = 1*m, mu1p = +2, sigma1p = 1
n = n1m + n1p + n0

				#SIMULATI0N DES SCOREs
T = zeros(1, n); Z = T;
for i=1:n1p, T(i) = -1; Z(i) = mu1p + sigma1p*randn(1,1); end
for i=n1p+1:n1p+n1m, T(i) = 1; Z(i) = mu1m + sigma1m*randn(1,1); end
for i=n1p+n1m+1:n, T(i) = 0; Z(i) = mu0 + sigma0*randn(1,1); end

#hist(Z(T==+1)); input('...');#hist(Z(T==0)); input('...');
#hist(Z(T==-1)); input('...');#hist(Z); input('...');

[FF1m ZZ1m] = hist(Z(T==-1)); 
[FF1p ZZ1p] = hist(Z(T==+1)); 
[FF0 ZZ0] = hist(Z(T==0)); 
plot(ZZ1m, FF1m, '1-', ZZ0, FF0, '2-', ZZ1p, FF1p, '3-'); 
    title('Répartitions séparées')
    input('...');	
[FF ZZ] = hist(Z); 
plot(ZZ, FF, '1-');
    title('Répartitions confondues')
    input('...');	
		
				# TRI PAR SCORE CROISSANT
[Z, I] = sort(Z);
T =T(I);
#[T' Z']
