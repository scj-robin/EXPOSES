clear all, clc, close all

% Param
n = 100, K = 4
tau = [ceil(n*sort(rand(1, K-1))) n]
mu = .5*(1+(-1).^(1:K))
%tau = [10 50 70 100]
%mu = 2*[0 1 0 1]
sigma = .2

% Simulation
l = diff([0 tau])
K = size(tau, 2)
n = max(tau)
X = sigma*randn(1, n);
m = [];
for k=1:K
  m = [m mu(k)*ones(1, l(k))];
end
X = X+m;

% Plot
subplot(221), 
plot((1:n), X, '.', (1:n), m, '-k', 'LineWidth', 2), hold on
F_PlotCGH(X, tau)
hold off

% Cost matrix
C = zeros(n, n);
for i=1:n-1
  for j=i:n
    C(i, j) = sum((X(i:j) - mean(X(i:j))).^2);
  end
end

% Dynamic programming
Kmax = 10
S = zeros(Kmax, n);
T = S;
S(1, :) = C(1, :);
for k=2:Kmax
  for i=k:n
    Stmp = S(k-1, 1:i-1) + C(2:i, i)';
    [S(k, i) T(k, i)] = min(Stmp);
  end  
end

% Model selection
Smin = min(S(:, n));
Smax = max(S(:, n));
subplot(222)
BICpen = log(n)*(1:Kmax)';
plot((1:Kmax), S(:, n), '-b', (1:Kmax), S(:, n)+BICpen, '-r', 'LineWidth', 2)
hold on
[Crit KBIC] = min(S(:, n)+BICpen);
plot([KBIC KBIC], [Smin Smax], ':k', 'LineWidth', 2)
D = (gammaln(n) - gammaln(n-(1:Kmax)+1) - gammaln((1:Kmax)))';
plot((1:Kmax), S(:, n), '-b', (1:Kmax), S(:, n)+BICpen+D, '-g', 'LineWidth', 2)
[Crit Kopt] = min(S(:, n)+BICpen+D);
plot([Kopt Kopt], [Smin Smax], ':k', 'LineWidth', 2)
hold off

% Best segmentation
Kopt = K
t = [T(Kopt, n) n];
for k=Kopt-1:-1:2
  t = [T(k, t(1)) t];
end
t
subplot(223), F_PlotCGH(X, t)

% MDS
C = C+C';
save('CGHsim.txt', 'C', 'ascii')
%C = max(max(-C))-C;
N = sum(C);
NN = sum(N);
N = (N - NN/2/n)/n;
XpX = .5*(repmat(N, n, 1) + repmat(N', 1, n) - C);
[V delta] = eig(XpX);
delta = diag(delta);
%subplot(224), plot(delta);
Id = num2str((1:n)');
subplot(224), plot(V(:, 1), V(:, 2), '.k'), hold on
text(V([1 t], 1), V([1 t], 2), Id([1 t], :)), hold off 

