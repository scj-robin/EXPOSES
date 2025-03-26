% FIGURES POUR L'EXPOSE
clear all, close all, clc, clf

% Loi exponentielle
lambda = 10000;
lambda2 = 25000
x = (1:1:10*lambda);
plot(x, exp(-x/lambda)/lambda, 'b-', [lambda lambda], [0, 1/lambda], 'b:', ...
   x, exp(-x/lambda2)/lambda2, 'r-', [lambda2 lambda2], [0, 1/lambda2], 'r:', ...
   'LineWidth', 2)
%saveas(1, 'DensExpo.eps', 'epsc')
print -depsc DensExpo.eps -f1
plot(log10(x), log10(x).*exp(-x/lambda)/lambda, 'b-', ...
   [log10(lambda) log10(lambda)], [0, log10(lambda)/lambda], 'b:', ...
   log10(x), log10(x).*exp(-x/lambda2)/lambda2, 'r-', ...
   [log10(lambda2) log10(lambda2)], [0, log10(lambda2)/lambda2], 'r:', ...
   'LineWidth', 2)
%saveas(1, 'DensExpo-EchLog.eps', 'epsc')
print -depsc DensExpo-EchLog.eps -f1

% Critères pour les mots
orient landscape
Souche = 'K12'
m = 5, P = 3, 
Seuil = 1
Rep = ['Q:\math_recherche\StatGenome\boucles\data_boucles\' Souche '\liste_boucles\'];
NomFicMots = [Rep 'Loop_' Souche '_triplemarkov_M' num2str(m) '_P' num2str(P) '_liste_mots.lst']
FicMots =fopen(NomFicMots, 'r');
Mot = [];
Data = [];
while feof(FicMots) == 0
%while length(Mot) < 500
   Ligne = fgetl(FicMots);
   Blanc = min(findstr(Ligne, ' '));
   Mot = [Mot; Ligne(1:Blanc-1)];
   Data = [Data; str2num(Ligne(Blanc:end))];
end
Mu = Data(:, 1);
MuP = Data(:, 2:1+P);
C2 = Data(:, 2+P:1+2*P);
C4 = Data(:, 2+2*P:1+3*P);
for p=1:P
%   subplot(ceil(P/3), 3, p),
   subplot(2, 3, p),
   Carac = ((abs(C2(:, p))>Seuil)+(abs(C4(:, p))>Seuil));
   plot(C2(:, p), C4(:, p), '+', [min(C2(:, p)) max(C2(:, p))], [0 0], '-.k', ...
      [0 0], [min(C4(:, p)) max(C4(:, p))], '-.k'), ...
      text(C2(Carac>=1, p), C4(Carac>=1, p), Mot(Carac>=1, :), 'fontsize', 5)
end
saveas(1, ['MotsCarac-' Souche '-M' num2str(m) '-P' num2str(P) '.eps'], 'epsc')
