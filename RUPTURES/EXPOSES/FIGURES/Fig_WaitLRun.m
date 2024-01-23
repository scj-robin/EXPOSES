% Fig waiting time for an $\ell$_run

Pi = [99.72 0.28; 2.26 97.74]/100
P = zeros(l+1, l+1)
P = [Pi(1, 1) Pi(1, 2) zeros(1, l-1); ...
    Pi(2, 1)*ones(l-1, 1)  zeros(l-1, 1) Pi(2, 2)*eye(l-1); ...
    zeros(1, l) 1]
    

l = 10
n = 500
mu1 = [1 zeros(1, l)];
mu = zeros(n, l+1);
mu(1, 1) = 1;
for t=2:n
    mu(t, :) = mu(t-1, :) * P;
end
plot((1:n), mu(:, l+1), 'b-', 'LineWidth', 2)
axis([0 n 0 1])
saveas(1, 'Fig_WaitLRun.eps', 'epsc')