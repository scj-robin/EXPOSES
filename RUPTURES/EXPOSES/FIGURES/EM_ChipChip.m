% Simulation de chip-chip (Buck & al.)

clear all, clc
cd D:\ENSEIGN\COURS\MELANGE\Exemples\ChipChip

% Simul
%mu = [-0.25 2.1]
%sigma = 0.6
%p = [0.8 0.2]
%n = 10000

%X = sigma*randn(n, 1);
%Z = (rand(n, 1) < p(1));
%X = Z*mu(1) + (1-Z)*mu(2) + X;
%X = sort(X);
%hist(X, ceil(sqrt(n)))

%save('EM_ChipChip.mat', 'X')

% Load
clear all
load('EM_ChipChip.mat', 'X');
n = length(X)
X2 = X.^2;
[HH XX] = hist(X, 2*ceil(sqrt(n)));
Max = max(HH)
Pas = mean(diff(XX))
subplot(211), hist(X, 2*ceil(sqrt(n)));
pause

% Init
eps = 1e-4
Z = (1:n)';
Xmin = min(X); Xmax = max(X); Xmed = median(X);
T = Xmed;
Z = (X < T);
%Z = 1+Z;
%[p, m, s2, tau, LI, LC, Iter] = F_EMGauss1(X, Z)
%return
Tau = [Z 1-Z];
subplot(212), plot(X, Tau(:, 1), 'g-', X, Tau(:, 2), 'r-', ...
    'LineWidth', 2)
pause

% E-M
Diff = 2*eps
while Diff > eps
    % M
    p_tmp = sum(Tau)/n;
    m_tmp = (X'*Tau) ./ (n*p_tmp);
    s2_tmp = (X2'*Tau) ./ (n*p_tmp) - m_tmp.^2
    s2_tmp = s2_tmp*p_tmp';
    
    % E
    P = ((repmat(X, 1, 2) - repmat(m_tmp, n, 1)).^2) / (2*s2_tmp);
    P = (exp(-P)/sqrt(2*pi*s2_tmp)) .* (repmat(p_tmp, n, 1));
    
    subplot(211), hist(X, 2*ceil(sqrt(n))), hold on;
    title(sprintf('Normal: p=%4.2f, m=%4.2f, s=%4.2f   /   Enriched: p=%4.2f, m=%4.2f, s=%4.2f', ...
        p_tmp(1), m_tmp(1), sqrt(s2_tmp), p_tmp(2), m_tmp(2), sqrt(s2_tmp)), ...
        'FontSize', 15), 
    %title(sprintf('prop=(%4.2f %4.2f),   mean=(%4.2f %4.2f),   std=%4.2f', ...
    %    p_tmp, m_tmp, sqrt(s2_tmp)), 'FontSize', 15), 
    pause
    
    plot(X, n*P(:, 1)*Pas, 'g-', X, n*P(:, 2)*Pas, 'r-', ...
        X, n*sum(P, 2)*Pas, 'y:', 'LineWidth', 2), 
    hold off;
    pause

    Tau_tmp = P ./ repmat(sum(P, 2), 1, 2);   
    subplot(212), plot(X, Tau_tmp(:, 1), 'g-', X, Tau_tmp(:, 2), 'r-', ...
        'LineWidth', 2)
    pause
        
    
    % Test
    Diff = max(max(abs(Tau - Tau_tmp)))
    
    % Update
    Tau = Tau_tmp;
    p = p_tmp;
    m = m_tmp;
    s2 = s2_tmp;
    disp([p m sqrt(s2)]);
    %pause(1)
end

    

