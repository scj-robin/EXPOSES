%function F_Select(Rep, Q, Init)
clear all, clc, 

%Rep = '..\..\Melange\Karate'
Rep = '..\..\Exemples\Karate'
Q = 4
Init = 'Ward'
LineWidth = 2
MarkerSize = 15
FontSize = 20

if Rep(end) ~= '/',  Rep = [Rep '/']; end
G = load([Rep 'Graphe']);
Id = G.Id;
P = load([Rep 'Erdos2_' Init '_Q' num2str(Q)]);
%P.alpha, P.pi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Effectifs des classes');
ZMAP = F_MAP(P.tau);
n = sum(ZMAP);
N = sum(n);
ncum = cumsum(n);
[x, y, z] = find(G.X);
Znum = ZMAP * (1:Q)';
disp('Définition d''un ordre');
Eq = P.tau*(1:Q)';
[Eq, t] = sort(Eq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'Dot-plot'
figure(1)
set(1, 'PaperOrientation', 'landscape');
set(1, 'PaperType', 'A4');
set(1, 'PaperPosition', [0.63 0.63 28.41 19.72]);
% Avant
subplot(121)
u = rand(1, N);
[u Rnd] = sort(u);
[x, y, z] = find(G.X(Rnd, Rnd));
plot(x, y, '.k', 'MarkerSize', MarkerSize), axis([0 N+1 0 N+1])
% Après
subplot(122)
[x, y, z] = find(G.X(t, t));
subplot(122)
plot(x, y, '.k', 'MarkerSize', MarkerSize), axis([0 N+1 0 N+1]), hold on
for q=1:Q, 
   plot([ncum(q)+0.5 ncum(q)+0.5], [1 N], 'b-', ...
      [1 N], [ncum(q)+0.5 ncum(q)+0.5], 'b-', ...
      'LineWidth', LineWidth);
end
hold 
%saveas(1, '../../Exposes/Figures/Karate-Dotplot.eps', 'epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
set(2, 'PaperOrientation', 'landscape');
set(2, 'PaperType', 'A4');
set(2, 'PaperPosition', [0.63 0.63 28.41 19.72]);
% Avant
subplot(121)
Rho = rand(1, N)'; Theta = rand(1, N)'*2*acos(-1);
Rho = ones(1, N)'; Theta = (1:N)'*2*acos(-1)/N;
[XY] = [Rho(Rnd).*cos(Theta(Rnd)) Rho(Rnd).*sin(Theta(Rnd))];
gplot(G.X, XY), axis off, hold on, 
text(XY(:, 1), XY(:, 2), G.Id(:, :), 'Fontsize', FontSize);
hold off
% Après
subplot(122), axis off
disp('Graphe');
D = log((1-P.pi)./P.pi);
D = D .* (ones(Q) - eye(Q));
Dbar = mean(D);
Dbbar = mean(Dbar);
XpX = -(D - repmat(Dbar, Q, 1) - repmat(Dbar', 1, Q) + ones(Q)*Dbbar);
[V L] = eig(XpX);
[V1 Ordre] = sort(V(:, 1));
XY = zeros(N, 2);
sigma = 0.5 * pi/sqrt(max(n))/Q
for q=1:Q
   V(q, 1) = cos(2*pi*Ordre(q)/Q);
   V(q, 2) = sin(2*pi*Ordre(q)/Q);
   ntmp = 0;
   for i=1:N
      if Znum(i) == q
         ntmp = ntmp+1;
         XY(i, 1) = V(q, 1)+sigma*sqrt(n(q))*cos(2*acos(-1)*ntmp/n(q));
         XY(i, 2) = V(q, 2)+sigma*sqrt(n(q))*sin(2*acos(-1)*ntmp/n(q));
         plot(XY(i, 1), XY(i, 2), '+k'); 
         hold on
      end
   end
end
gplot(G.X, XY), axis off %, hold on, 
text(XY(:, 1), XY(:, 2), G.Id(:, :), 'Fontsize', FontSize);
%axis([-1.5 1.5 -1.5 1.5]), hold on
%for q=1:Q
%   text(1.1*V(q, 1), 1.1*V(q, 2), sprintf('{\\bf %i} (%i)', q, n(q)));
%end
hold off
saveas(2, '../../Exposes/Figures/Karate-Graph.eps', 'epsc')

%% To add circles around groups: does no work
%Theta = 2*acos(-1)*(0:100)/100;
%C = [cos(Theta)' sin(Theta)'];
%Shift = .5
%Coef = 2
%Color = 'brgm'
%for q=1:Q
%   V(q, 1) = cos(2*pi*Ordre(q)/Q)+ Shift;
%   V(q, 2) = sin(2*pi*Ordre(q)/Q);
%   plot(V(q, 1) + Coef*sigma*sqrt(n(q))*C(:,1), ...
%       V(q, 2) + Coef*sigma*sqrt(n(q))*C(:,2), ...
%       [Color(q) '-'], 'LineWidth', LineWidth)   
%end
hold off
saveas(2, '../../Exposes/Figures/Karate-Graph-Circle.eps', 'epsc')

