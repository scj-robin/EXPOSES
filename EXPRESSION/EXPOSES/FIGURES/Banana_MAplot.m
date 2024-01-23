disp('Banana shape in M-A plots');

cd D:\RECHERCHE\EXPRESSION\EXPOSES\Figures

X = (3:0.005:8);
n = length(X)
s = .1
delta = 50

disp('Data simulation');
R = X + s*randn(1, n);
G = log(exp(X)+delta) + s*randn(1, n);
A = R-G; M = R+G;

subplot(211), plot(M, A, 'b+', [min(M) max(M)], [0 0], 'r-')
subplot(212), plot(M, exp(G) - exp(R), 'b+')

deltahat = abs(exp(G) - exp(R))*(exp(G) - exp(R))' / sum(abs(exp(G) - exp(R)))

disp('Max likelihood estimate of delta');

d = delta;
for i=1:n
   
end
