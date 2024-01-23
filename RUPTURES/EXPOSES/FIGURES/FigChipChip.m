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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMUL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
pi = 0.2
sigma = 0.5
Data = load([DataDir DataName '.dat']);
size(Data)
X = 1.5+Data(:, 5)';
X = sort(X);
n = length(X);
%X = 7+9*rand(1, n);
Z = (rand(1, n) < pi);
E = sigma*randn(1, n);

% Bon cas
a0 = 2.8, a1 = a0+2, b0 = 0.7, b1 = b0
Y = (a0+b0*X).*(1-Z) + (a1+b1*X).*Z + E;
R = Y-X;
subplot(221), plot(X, Y, '.k'), hold on, axis([5 16 5 16]), ...
   plot(X, a0+b0*X, '-g', X, a1+b1*X, '-r', 'LineWidth', 1.5), hold off
%subplot(221), plot(X(Z==1), Y(Z==1), '.r', X(Z==0), Y(Z==0), '.k')
subplot(222), hist(R, ceil(sqrt(n))), axis([-4 4 0 600])

% Mauvais cas
a0 = 2.8, a1 = a0-0.1, b0 = 0.6, b1 = 0.8
Y = (a0+b0*X).*(1-Z) + (a1+b1*X).*Z + E;
R = Y-X;
subplot(223), plot(X, Y, '.k'), hold on, axis([5 16 5 16]), ...
   plot(X, a0+b0*X, '-g', X, a1+b1*X, '-r', 'LineWidth', 1.5), hold off
%subplot(223), plot(X(Z==1), Y(Z==1), '.r', X(Z==0), Y(Z==0), '.k')
subplot(224), hist(R, ceil(sqrt(n))), axis([-4 4 0 600])

saveas(1, [FigDir 'SimHistoChIPChip.eps'], 'epsc')

close all, return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

