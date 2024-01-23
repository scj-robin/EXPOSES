%clf(1), clf(2), , clc
% Exemple pour model segmentation/classif

%cd D:\RECHERCHE\CGH\Exposes\Figures
Rep = 'D:\RECHERCHE\CGH\Exposes\Figures\'

%Param simul
L = [10 5 5 3 15] % lg des segments
Ltot = sum(L)
ZK = [1 2 1 2 1] % groupe des segments
K = length(ZK), P = max(ZK)
mu = [0 2.5]
sigma= [1 1]
X = (1:Ltot)

% param figure
Color = 'krbmg'
Marker = 'ov+sx'
LineWidth = 2
MarkerWidth = 2
MarkerSize = 5
t = 1

%simulation
Y = []; Z = []; S = []; 
for k=1:K
   Z = [Z ZK(k)*ones(1, L(k))];
   S = [S k*ones(1, L(k))];
   Ytmp = mu(ZK(k))+sigma(ZK(k))*randn(1, L(k));
   Y = [Y Ytmp];
end
Y

% calcul des moments par segments
for k=1:K
   m(k) = mean(Y(S==k));
   s(k) = std(Y(S==k));
end
[m; s]

% calcul des moments par population
for p=1:P
   mu(p) = mean(Y(Z==p));
   sigma(p) = std(Y(Z==p));
end
[mu; sigma]

figure(1)
hold on
for k=1:K
   plot(X(S==k), Y(S==k), [Color(k) Marker(k)], 'MarkerSize', MarkerSize, 'LineWidth', MarkerWidth)
   plot(X(S==k), m(k)*ones(1, L(k)), [Color(k) '-'], 'LineWidth', LineWidth)
   plot(X(S==k), (m(k)+t*s(k))*ones(1, L(k)), [Color(k) '-.'], 'LineWidth', LineWidth)
   plot(X(S==k), (m(k)-t*s(k))*ones(1, L(k)), [Color(k) '-.'], 'LineWidth', LineWidth)
end
hold off
saveas(1, [Rep 'FigSegClas-1.eps'], 'epsc')

figure(2)
hold on
for k=1:K
   p = Z(sum(L(1:k))-1);
   plot(X(S==k), Y(S==k), [Color(p) Marker(p)], 'MarkerSize', MarkerSize, 'LineWidth', MarkerWidth)
   plot(X(S==k), mu(p)*ones(1, L(k)), [Color(p) '-'], 'LineWidth', LineWidth)
   plot(X(S==k), (mu(p)+t*sigma(p))*ones(1, L(k)), [Color(p) '-.'], 'LineWidth', LineWidth)
   plot(X(S==k), (mu(p)-t*sigma(p))*ones(1, L(k)), [Color(p) '-.'], 'LineWidth', LineWidth)
end
hold off
saveas(2, [Rep 'FigSegClas-2.eps'], 'epsc')

