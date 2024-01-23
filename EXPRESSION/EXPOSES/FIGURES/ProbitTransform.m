% Figure illustrant la transformation probit

% paramètres
n = 1000, EF = 1e-3, p1 = 0.1
n1 = round(p1*n), n0 = n-n1

% simulation des P
P = [-EF*log(rand(1, n1)) rand(1, n0)];
P = sort(P);
[NN, PP] = hist(P, ceil(sqrt(n)));
PasHist = mean(diff(P));
[PP, NN] = stairs([min(PP)-PasHist PP max(PP)+PasHist], [0 NN 0]);
figure(1)
NNmax = max(NN)/n/PasHist;
plot(PP, NN/n/PasHist, 'k', 'LineWidth', 2);
saveas(1, 'D:\RECHERCHE\EXPRESSION\EXPOSES\Figures\ProbitTransform-P.eps', 'epsc')
   
% calcul des X
X = sqrt(2)*erfinv(2*P-1);
[NN, XX] = hist(X, ceil(sqrt(n)));
PasHist = mean(diff(XX));
[XX, NN] = stairs([min(XX)-PasHist XX max(XX)+PasHist], [0 NN 0]);
figure(2)
plot(XX, NN/n/PasHist, 'k', 'LineWidth', 2),
saveas(2, 'D:\RECHERCHE\EXPRESSION\EXPOSES\Figures\ProbitTransform-X.eps', 'epsc')

% Loi normale
figure(3)
Phi = 0.5*(erf(X/sqrt(2))+1);
plot(X, Phi, 'k', 'LineWidth', 2),
saveas(3, 'D:\RECHERCHE\EXPRESSION\EXPOSES\Figures\ProbitTransform-Phi.eps', 'epsc')



