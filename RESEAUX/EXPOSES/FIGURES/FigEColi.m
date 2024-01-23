% Figure EColo-AvecClique -> Article

Rep = '..\EColi-AvecClique\'
Qmax = 13

G = load([Rep 'Graphe'])
n = length(G.X)

% LC, BIC
for Q=1:Qmax
   P = load([Rep 'Erdos2_Q' num2str(Q)]);
   Lc(Q) = P.LC;
   pBIC(Q) = Lc(Q) - log(n)*(Q-1) - log((n*(n-1)/2))*(Q*(Q+1));
end
[pBICmax QpBIC] = max(pBIC)
figure(1)
Qmin = 2
plot((Qmin:Qmax), Lc(Qmin:Qmax), 'k-', (Qmin:Qmax), pBIC(Qmin:Qmax), 'b--', ...
   [QpBIC QpBIC], [min(pBIC(Qmin:Qmax)) max(Lc(Qmin:Qmax))], 'r:', 'LineWidth', 2), xlabel('Q')
saveas(1, 'FigEColi-LcBIC.eps', 'epsc')

% Afficahge des paramètres
P = load([Rep 'Erdos2_Q' num2str(QpBIC)])
round(10000*P.alpha)/100
round(10000*P.pi)/100
pause

% Ajustement du modèle ERMG
F_GPlot(Rep, QpBIC)
saveas(1, 'FigEColi-Graph.eps', 'epsc')
saveas(1, 'FigEColi-Graph.fig', 'fig')
%saveas(2, 'FigEColiERMG.eps', 'epsc')
saveas(2, 'FigEColi-ERMG.fig', 'fig')

% Connectivité entre les classes
P.alpha, P.pi
Atheo = n*(n-1)*(P.alpha'*P.alpha).*P.pi/2
Atau = P.tau'*G.X*P.tau/2
ZMAP = F_MAP(P.tau);
AMAP = ZMAP'*G.X*ZMAP/2
Comp = round(triu(Atheo)) + round(tril(AMAP) - diag(diag(AMAP)))
diag(AMAP)'


