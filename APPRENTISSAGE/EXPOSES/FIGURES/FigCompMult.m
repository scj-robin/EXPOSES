% Figure exposés comparaison multiples

clear all, clc
cd('D:\RECHERCHE\APPRENTISSAGE\EXPOSES\Figures')

% Histogramme Storey
figure(1)
m = 10000
pi0 = .9
m0 = pi0*m;
E0 = 3e-2
m1 = m-m0;
P = [-E0*log(rand(1, m1)) rand(1, m0)];
P = sort(P);
[H, B] = hist(P, ceil(sqrt(m)));
delta = mean(diff(B))
[HH, BB] = F_Histo(H, B);
plot(BB, HH, 'k-', 'LineWidth', 2), hold on
plot([0 1], [m*delta*pi0 m*delta*pi0], 'r-', 'LineWidth', 2)
plot([.5 .5], [0 .9*max(H)], 'b-.', 'LineWidth', 2)
text(.48, .95*max(H), '\lambda^*', 'FontSize', 30)
saveas(1, 'FigStorey.fig')
saveas(1, 'FigStorey.eps', 'epsc')
saveas(1, 'FigStorey.png', 'png')
hold off

% Histogramme Celisse
figure(2)
lambda = 15*delta
lambda = max(BB(BB <= lambda))
N0 = sum(P>lambda)
h0 = delta*N0/(1-lambda);
plot([BB(BB<lambda); lambda], [HH(BB<lambda); h0], 'k-', ...
    [lambda 1], h0*[1 1], 'k-', ...
    'LineWidth', 2), hold on
plot([lambda lambda], [0 .9*max(H)], 'b-.', 'LineWidth', 2)
text(lambda-.02, .95*max(H), '\lambda', 'FontSize', 30)
text(.5*(1+lambda)-.02, .5*delta*m0, '(1-\lambda)\pi_o', 'FontSize', 30)
saveas(2, 'FigCelisse.fig')
saveas(2, 'FigCelisse.eps', 'epsc')
saveas(2, 'FigCelisse.png', 'png')
hold off

% Histogramme Celisse2
figure(3)
m00 = m0-m1;
P2 = [-E0*log(rand(1, m1)) rand(1, m00) 1+E0*log(rand(1, m1))];
[H, B] = hist(P2, ceil(sqrt(m)));
delta = mean(diff(B))
[HH, BB] = F_Histo(H, B);
lambda = 15*delta
lambda = max(BB(BB <= lambda))
mu = 1-12*delta
mu = min(BB(BB >= mu))
N00 = sum((lambda<=P2).*(P2<=mu))
h00 = delta*N00/(mu-lambda);
plot([BB(BB<lambda); lambda], [HH(BB<lambda); h00], 'k-', ...
    [lambda mu], h00*[1 1], 'k-', ...
    [mu; BB(BB>mu)], [h00; HH(BB>mu)], 'k-', ...
    'LineWidth', 2), hold on
plot([lambda lambda], [0 .9*max(H)], 'b-.', 'LineWidth', 2)
text(lambda-.02, .95*max(H), '\lambda', 'FontSize', 30)
plot([mu mu], [0 .9*max(H)], 'b-.', 'LineWidth', 2)
text(mu-.02, .95*max(H), '\mu', 'FontSize', 30)
text(.5*(mu+lambda)-.02, .5*delta*m0, '(\mu-\lambda)\pi_o', 'FontSize', 30)
saveas(3, 'FigCelisse2.fig')
saveas(3, 'FigCelisse2.eps', 'epsc')
saveas(3, 'FigCelisse2.png', 'png')
hold off

