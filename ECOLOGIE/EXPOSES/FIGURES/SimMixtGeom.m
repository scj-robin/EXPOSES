% Param
K = 3;  
n = 1000
pi = (K:-1:1);   
pi = pi/sum(pi)
gamma = .6.^(1:K)
LineWidth = 3
Color = ['grb']

% Simul
Z = rand(n, 1);
Z = (repmat(Z, 1, K) < repmat(cumsum(pi), n, 1));
Z = K-Z*ones(K, 1)+1;
X = floor(-log(rand(1, n))./gamma(Z));

% Distrib
Xmax = max(X)
np = sum(X > 0);
f = zeros(K, 1+Xmax);
g = zeros(1, 1+Xmax);
for k=1:K
    f(k, :) = gamma(k)*(1-gamma(k)).^(0:Xmax);
    g = g + pi(k) * f(k, :);
end

% Plot
figure(1)
C = zeros(1, 1+Xmax);
for x=0:Xmax
    C(x+1) = sum(X==x);
end
plot((1:Xmax)', n*g(2:end), '-k', 'LineWidth', LineWidth), axis([0 Xmax 0 ceil(1.1*max(C))]), hold on
plot((0:1), n*g(1:2), '.--r', 'LineWidth', LineWidth), 
plot((1:Xmax), C(2:end), 'ok', 'LineWidth', LineWidth), 
plot((0), C(1), '*r', 'LineWidth', LineWidth), hold off
saveas(1, 'SimMixtGeom.eps', 'epsc')
saveas(1, 'SimMixtGeom.jpg', 'jpeg')
saveas(1, 'SimMixtGeom.png', 'png')

% Plot
figure(2)
plot((1:Xmax)', n*g(2:end), '-k', 'LineWidth', LineWidth), axis([0 Xmax 0 ceil(1.1*max(C(2:end)))]), hold on
plot((1:Xmax), C(2:end), 'ok', 'LineWidth', LineWidth), 
for k=1:K
	plot((1:Xmax)', n*f(k, 2:end), ['--' Color(k)], 'LineWidth', LineWidth), 
end
hold off
saveas(2, 'SimMixtGeomComp.eps', 'epsc')
saveas(2, 'SimMixtGeomComp.jpg', 'jpeg')
saveas(2, 'SimMixtGeomComp.png', 'png')

return 

% Plot log
figure(3)
C = zeros(1, 1+Xmax);
for x=0:Xmax
    C(x+1) = sum(X==x);
end
semilogy((1:Xmax), n*g(2:end), '--k', 'LineWidth', LineWidth), axis([0 Xmax 0 ceil(1.1*max(C))]), hold on
semilogy((0:1), n*g(1:2), '.--r', 'LineWidth', LineWidth), 
semilogy((1:Xmax), C(2:end), 'ok', 'LineWidth', LineWidth), 
semilogy((0), C(1), '*r', 'LineWidth', LineWidth), hold off
saveas(3, 'LogSimMixtGeom.eps', 'epsc')
saveas(3, 'LogSimMixtGeom.jpg', 'jpeg')
saveas(3, 'LogSimMixtGeom.png', 'png')

