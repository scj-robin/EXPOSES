% Figure exposé

clear all, clc, clf
figure(1)

% Données
DataDir = 'D:\RECHERCHE\RUPTURES\ChIP-Chip\Resultat\', 
FigDir = 'D:\RECHERCHE\RUPTURES\Exposes\Figures\', 
DataName = 'tfl2_1_ResReg'

% Paramètres
LineWidth = 3
FontSize = 20

% Lecture des données
Data = load([DataDir DataName '.dat']);
size(Data)
X = Data(:, 5);
Y = Data(:, 6);
Tau = Data(:, 7);
n = length(X);
clear Data

% Données brutes
plot(X, Y, '.k'), xlabel('Input', 'FontSize', FontSize), ylabel('IP', 'FontSize', FontSize)
saveas(1, [FigDir DataName '-RawIPInput.eps'], 'epsc');

% LogRatio
hist(Y-X, ceil(sqrt(n))), xlabel('IP/Input', 'FontSize', FontSize)
saveas(1, [FigDir DataName '-LogRatio.eps'], 'epsc');

% DistPosterior
T = sort(Tau);
plot((1:n), T, 'r', LineWidth', LineWidth), axis([1 n 0 1]), ...
   xlabel('Number of probes', 'FontSize', FontSize), ylabel('Posterior \tau', 'FontSize', FontSize)
saveas(1, [FigDir DataName '-PosteriorDist.eps'], 'epsc');

return

% Posterior
pcolor(X, Y, Tau), colorbar, xlabel('Input', 'FontSize', FontSize), ylabel('IP', 'FontSize', FontSize)
saveas(1, [FigDir DataName '-Posterior.eps'], 'epsc');

