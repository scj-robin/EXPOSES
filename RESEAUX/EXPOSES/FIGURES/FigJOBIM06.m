%function F_Select(Rep, Q, Init)

Rep = '../EColi/Complet'
Q = 21
Init = 'Ward'

%close all
LineWidth = 2

if Rep(end) ~= '/',  Rep = [Rep '/']; end
G = load([Rep 'Graphe'])
Id = G.Id;
P = load([Rep 'Erdos2_' Init '_Q' num2str(Q)])
%P.alpha, P.pi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Effectifs des classes');
ZMAP = F_MAP(P.tau);
n = sum(ZMAP)
N = sum(n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'Mean degree'
Kmean = G.K * P.tau ./ sum(P.tau);
EK = (N-1) * (P.pi * P.alpha')'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Degrees');
P.pi, P.alpha
Kmax = max(full(G.K))
KK = (0:Kmax);
PoisMixt = repmat(EK', 1, Kmax+1).^repmat(KK, Q, 1);
PoisMixt = PoisMixt ./ repmat(gamma(KK+1), Q, 1);
PoisMixt = PoisMixt .* repmat(exp(-EK'), 1, Kmax+1);
PoisMixt = P.alpha * PoisMixt
% Histogram
[Nhist Khist]= F_HistoSR(G.K);
subplot(121), plot(Khist, Nhist, 'b-'), hold on
plot(KK, N*PoisMixt, 'k-', 'LineWidth', LineWidth), hold off
% PP plot
[Freq, KK] = hist(full(G.K), Kmax+1)
Freq = Freq / N
subplot(122), plot([0 1], [0 1], 'b-', 'LineWidth', 1), hold on
plot(cumsum(Freq), cumsum(PoisMixt), 'k-+', 'LineWidth', LineWidth);, hold off
 
