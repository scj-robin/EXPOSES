function F_PlotCGH(X, tau)
  
n = size(X, 2);
K = size(tau, 2);
l = diff([0 tau]);
m = mean(X(1:tau(1))) * ones(1, l(1));
Xmin = min(X);
Xmax = max(X);
for k=2:K
  m = [m mean(X(tau(k-1)+1:tau(k))) * ones(1, l(k))];
end

plot((1:n), X, '.', (1:n), m, '-r', 'LineWidth', 2), hold on
for k=1:K-1
  plot([tau(k)+.5 tau(k)+.5], [Xmin Xmax], 'b:', 'LineWidth', 2)
end
hold off