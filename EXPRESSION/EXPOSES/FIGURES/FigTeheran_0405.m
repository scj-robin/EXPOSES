clear all, clc, clf
'Figure Exposes-Teheran'

orient landscape

% paramètres
n = 500
a = 0.3
t = 0.5

% simulation
P1 = rand(1, a*n)*10^(-4);
P0 = rand(1, (1-a)*n);
P = sort([P0 P1]);

% estimation sur histogramme
nbin = ceil(sqrt(n))
p0_est1 = (sum(P>t)/n)/(1-t)
[X N] = hist(P);
string = ['pi0 = ', num2str(p0_est1)]
subplot(1, 2, 1), hist(P, nbin), axis([0 1 0 max(X)]), hold on, ...
   plot([0 1], [p0_est1*n/nbin p0_est1*n/nbin], 'r-', 'LineWidth', 2), ...
   text(0.55, 0.6*max(X), string), ...
   plot([t t], [0 max(X)], 'b:', 'LineWidth', 2), hold off;


% estimation sur f° de répartition
G = (1:n)/n;
[PP GG] = stairs(P, G);
Theta = polyfit(P(P>t), G(P>t), 1)
p0_est2 = Theta(1)
string = ['pi0 = ', num2str(p0_est2)]
subplot(1, 2, 2), ...
   plot(PP, GG, '-k', 'LineWidth', 2), axis([0 1 0 1]), hold on, ...
   plot(P, polyval(Theta, P), '-r', 'LineWidth', 2), ...
   text(0.55, 0.4, string), ...
   plot([t t], [0 1], ':b', 'LineWidth', 2), hold off;


saveas(1, 'RegGenoWas.eps', 'epsc')

return

figure(2)
[N0, X0] = hist(P0, (0.05:0.05:0.95));
[XX0, NN0] = stairs([min(X0)-mean(diff(X0)) X0 max(X0)+mean(diff(X0))], [0 N0 0]);
[N, X] = hist(P, (0.05:0.05:0.95));
[XX, NN] = stairs([min(X)-mean(diff(X)) X max(X)+mean(diff(X))], [0 N 0]);
plot(XX0, NN0, 'g-', XX, NN, 'r-', 'LineWidth', 1)
