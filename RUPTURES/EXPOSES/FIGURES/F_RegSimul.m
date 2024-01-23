function [a, b, sigma2] = F_RegSimul(Y, X)

% Régression simultanées de toutes les colones de Y sur X
% Même pente pout toutes les colonnes mais constantes différentes

[n, p] = size(Y);
YY = reshape(Y, n*p, 1);
XX = [];
for j=1:p
   XX = [XX; ...
         zeros(n, j-1) ones(n, 1) zeros(n, p-j)];
end
XX = [XX repmat(X, p, 1)];
theta = linsolve(XX, YY);
a = theta(1:p);
b = theta(p+1);
sigma2 = norm(YY - XX*theta)/n/p

clear XX YY theta

