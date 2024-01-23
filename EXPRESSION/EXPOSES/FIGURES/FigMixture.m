% Mélange de 2 gaussiennes
clear all, clc, clf

%cd D:\RECHERCHE\EXPRESSION\EXPOSES\Figures

orient landscape

% paramètres
m = 1.2, rho = 0.3, v1 = 1, v2 = 1.5
mu1 = [+m +m],  Sigma1 = v1*[1 +rho; +rho 1],  n1 = 20
mu2 = [-m -m],  Sigma2 = v2*[1 -rho; -rho 1],  n2 = 20
n = n1+n2
pi = [n1 n2]/n

% Simlation
X1 = randn(n1, 2)*chol(Sigma1) + repmat(mu1, n1, 1);
X2 = randn(n2, 2)*chol(Sigma2) + repmat(mu2, n2, 1);
%save('FigMixtureGauss2.mat', 'mu1', 'mu2', 'Sigma1', 'Sigma2', 'X1', 'X2')
Tmp = load('FigMixtureGauss2.mat', 'mu1', 'mu2', 'Sigma1', 'Sigma2', 'X1', 'X2')
X1 = Tmp.X1; X2 = Tmp.X2;
%mean(X1), cov(X1), mean(X2), cov(X2)

% Estimation mélange Gaussien
X = [X1; X2];
%[criterion, model, strategy] = mixmodInput();
strategy.tab.Algo.name = 'EM';
Gaus2 = mixmod(X, 2, [], {'Gaussian_pk_Lk_Ck'})%, strategy)
p = Gaus2.modelOutput.param.proportion
m1 = Gaus2.modelOutput.param.mean(1, :)
m2 = Gaus2.modelOutput.param.mean(2, :)
S1 = Gaus2.modelOutput.param.variance(:, :, 1)
S2 = Gaus2.modelOutput.param.variance(:, :, 2)

% K-means (en fait mélange gaussien V = I)
strategy = Gaus2.condExe.strategy;
strategy.tabAlgo.name = 'CEM';
strategy.initType.tabLabel = Gaus2.modelOutput.proba.partition;
Kmean2 = mixmod(X, 2, [], {'Gaussian_p_L_I'}, strategy);
S = Kmean2.modelOutput.param.variance(:, :, 1)
k1 = Kmean2.modelOutput.param.mean(1, :)
k2 = Kmean2.modelOutput.param.mean(2, :)

% Elipses
Coef = 1;
theta = (0:0.1:6.5)';
Cercle = Coef*[cos(theta) sin(theta)];
E1 = Cercle*chol(S1) + repmat(m1, size(Cercle, 1), 1);
E2 = Cercle*chol(S2) + repmat(m2, size(Cercle, 1), 1);
K1 = Cercle*chol(S) + repmat(k1, size(Cercle, 1), 1);
K2 = Cercle*chol(S) + repmat(k2, size(Cercle, 1), 1);

% proba a posteriori
iS1 = inv(S1); iS2 = inv(S2); dS1 = det(S1); dS2 = det(S2);
xmin = floor(min([X1(:, 1); X2(:, 1)])); xmax = ceil(max([X1(:, 1); X2(:, 1)])); 
ymin = floor(min([X1(:, 2); X2(:, 2)])); ymax = ceil(max([X1(:, 2); X2(:, 2)])); 
lim = max(abs([xmin xmax ymin ymax]));
x = (-lim : 0.1 : lim); 
y = (-lim : 0.1 : lim); 
for i=1:length(x)
   for j=1:length(y)
      X = [x(i) y(j)];
      f1 = p(1)*exp(-0.5*(X-m1)*iS1*X')/sqrt(dS1);
      f2 = p(2)*exp(-0.5*(X-m2)*iS2*X')/sqrt(dS2);
      f(i, j) = f1+f2;
      tau1(i, j) = f1 / f(i, j);
      tau2(i, j) = f2 / f(i, j);
      Kmean(i, j) = (k2-k1)*([x(i) y(j)]-(0.5*k1+0.5*k2))';
   end
end
if sum(m1)< sum(m2), Kmean = -Kmean; end


% Figures
LineWidth = 2
Axes = [xmin xmax ymin ymax];
subplot(221), surf(x, y, f'), axis([Axes 0 max(max(f))])
subplot(222), ...
   contourf(x, y, tau1',9), axis(Axes), colorbar, shading flat, hold on, ...
   plot(X1(:, 1), X1(:, 2), 'ok', X2(:, 1), X2(:, 2), 'ok', 'LineWidth', LineWidth), axis(Axes), ...
   plot(m1(1), m1(2), '+r', m2(1), m2(2), '+b', 'LineWidth', LineWidth), axis(Axes), ...
   hold off 
subplot(223), ...
   contourf(x, y, tau1', [0.5 0.5]), axis(Axes), hold on, ... 
   plot(X1(:, 1), X1(:, 2), 'ok', X2(:, 1), X2(:, 2), 'ok', 'LineWidth', LineWidth), axis(Axes), ...
   plot(m1(1), m1(2), '+r', m2(1), m2(2), '+b', 'LineWidth', LineWidth), axis(Axes), ...
   plot(E1(:, 1), E1(:, 2), '-r', E2(:, 1), E2(:, 2), '-b', 'LineWidth', LineWidth), axis(Axes), ...
   hold off
subplot(224), ...
   contourf(x, y, Kmean',[0 0]), axis(Axes), hold on, ...
   plot(X1(:, 1), X1(:, 2), 'ok', X2(:, 1), X2(:, 2), 'ok', 'LineWidth', LineWidth), axis(Axes), ...
   plot(k1(1), k1(2), '+r', k2(1), k2(2), '+b', 'LineWidth', LineWidth), axis(Axes), ...
   plot(K1(:, 1), K1(:, 2), ':r', K2(:, 1), K2(:, 2), ':b', 'LineWidth', LineWidth), axis(Axes), ...
   hold off
saveas(1, 'FigMixtureGauss2.eps', 'epsc')



