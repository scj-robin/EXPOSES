% Figure lois Beta

clear all, clf, clc

% Exemples de densités
x = (eps:0.01:1-eps);
r = [1 1 3 1/3]
s = [1 10 10 1/3]
for i=1:length(r)
   b(i, :) = (x.^(r(i)-1)).*(1-x).^(s(i)-1)/Beta(r(i), s(i));
end
plot(x, b(1, :), 'k-', x, b(2, :), 'r-', x, b(3, :), 'b-', ...
   x, b(4, :), 'g-', 'LineWidth', 3), axis([0 1 0 10])
saveas(1, 'FigBeta.eps', 'epsc')
%pause

% probas a posteriori
a = 0.1, r = 1, s = 10
b = (x.^(r-1)).*(1-x).^(s-1)/Beta(r, s);
B = Betainc(x, r, s);
f = a * b + (1-a);
F = a*B + (1-a)*x;
tau = a * b ./ f;
FDR = (1-a) * x ./ F;
FNR = a * (1-B) ./ (1-F);
plot(x, a*b, x, ones(size(x))*(1-a), x, f, 'LineWidth', 3)
plot(x, tau, x, FDR, x, FNR, 'LineWidth', 3), axis([0 1 0 1]);
saveas(1, 'FigBetaPosterior.eps', 'epsc')

