clear all, clc, clf
% Figure article

cd D:\RECHERCHE\RESEAUX\EXPOSES\figures
figure(1), axis off
set(1, 'PaperOrientation', 'landscape');
set(1, 'PaperType', 'A4');
%set(1, 'Fontsize', 20);
set(1, 'PaperPosition', [0.63 0.63 28.41 19.72]);
Symbol = ['o' '^' '*' 'v']
Color = ['b' 'r' 'g' 'k']
MarkerSize = 20
Fontsize = 20
LineWidth = 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Etoiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = 5 % taille de chaque étoile
R = 0.3 % rayon de chaque étoile
X = zeros(2*C+2, 2*C+2);
Z = [ones(1, C) 2 3 4*ones(1, C)]
Coord = zeros(2*C+2, 2);
Coord(C+1, :) = [0 0];
Coord(C+2, :) = [1 0];
X(C+1, C+2) = 1;
for c=1:C
   X(c, C+1) = 1;
   Coord(c, :) = [R*cos(2*pi*(c+0.5)/C) R*sin(2*pi*(c+0.5)/C)];
   X(C+2, C+2+c) = 1;
   Coord(C+2+c, :) = [1 + R*cos(2*pi*c/C) R*sin(2*pi*c/C)];
end
X = X+X';
figure(1)%, subplot(221)
gplot(X, Coord, 'k-'), hold on, axis off%, axis([-0.5 1.5 -0.4 0.4])
%saveas(1, 'FigNetworks-Star', 'epsc')
%saveas(1, 'FigNetworks-Star', 'fig')
saveas(1, 'FigNetworks-Star', 'png')
for i=1:length(Z)
   plot(Coord(i, 1), Coord(i, 2), ...
      [Symbol(Z(i)) Color(Z(i))], 'MarkerSize', MarkerSize, 'LineWidth', LineWidth)
end
hold off
%saveas(1, '../../Articles/StatComput/Figures/FigNetworks-Star-Symb', 'epsc')
%saveas(1, '../../Articles/StatComput/Figures/FigNetworks-Star-Symb', 'fig')
%saveas(1, 'FigNetworks-Star-Col', 'epsc')
%saveas(1, 'FigNetworks-Star-Col', 'fig')
saveas(1, 'FigNetworks-Star-Col', 'png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clusters = 2 cliques
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = 5 % taille de chaque clique
X = [ones(C, C) zeros(C, C); zeros(C, C) ones(C, C)];
X(C, C+1) = 1;
Z = [ones(1, C) 2*ones(1, C)]
Coord = zeros(2*C, 2);
for c=1:C
   Coord(c, :) = [R*cos(2*pi*(c+0.5)/C) R*sin(2*pi*(c+0.5)/C)];
   Coord(C+c, :) = [1 + R*cos(2*pi*c/C) R*sin(2*pi*c/C)];
end
X = X+X';
%figure(1), subplot(222)
gplot(X, Coord, 'k-'), hold on, axis off %, axis([-0.5 1.5 -0.4 0.4])
%saveas(1, 'FigNetworks-Clusters', 'epsc')
%saveas(1, 'FigNetworks-Clusters', 'fig')
saveas(1, 'FigNetworks-Clusters', 'png')
for i=1:length(Z)
   plot(Coord(i, 1), Coord(i, 2), ...
      [Symbol(Z(i)) Color(Z(i))], 'MarkerSize', MarkerSize, 'LineWidth', LineWidth)
end
hold off
%saveas(1, '../../Articles/StatComput/Figures/FigNetworks-Clusters-Symb', 'epsc')
%saveas(1, '../../Articles/StatComput/Figures/FigNetworks-Clusters-Symb', 'fig')
%saveas(1, 'FigNetworks-Clusters-Col', 'epsc')
%saveas(1, 'FigNetworks-Clusters', 'fig')
saveas(1, 'FigNetworks-Clusters', 'png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Product connection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = 5 % taille des groupes
a = 0.8, b = 0.2
X = zeros(2*C, 2*C);
Z = [ones(1, C) 2*ones(1, C)]
Coord = zeros(2*C, 2);
for c=1:C
   X(c, C+1:2*C) = (rand(1, C) < a*b);
   Coord(c, :) = [R*cos(2*pi*(c+0.5)/C) R*sin(2*pi*(c+0.5)/C)];
   Coord(C+c, :) = [1 + R*cos(2*pi*c/C) R*sin(2*pi*c/C)];
   for d = c+1:C
      X(c, d) = (rand(1) < a^2);      
      X(C+c, C+d) = (rand(1) < b^2);
   end
   if sum(X(c, :))==0, X(c, c) = 1; end
   if sum(X(C+c, :))==0, X(C+c, C+c) = 1; end
end
X11 = [0 1 1 0 0; 0 0 1 0 1; 0 0 0 1 1; 0 0 0 0 1; 0 0 0 0 0]; 
X12 = [0 0 0 1 0; 0 0 1 0 0; 0 1 0 0 1; 0 0 0 1 0; 0 1 0 0 0]; 
X22 = [0 0 0 1 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
X = [X11 X12; X12' X22];
X = X+X';
%figure(1), subplot(223)
gplot(X, Coord, 'k-'), axis off, hold on%, axis([-0.5 1.5 -0.4 0.4])
%saveas(3, 'FigNetworks-Indep', 'epsc')
%saveas(3, 'FigNetworks-Indep', 'fig')
saveas(3, 'FigNetworks-Indep', 'png')
for i=1:length(Z)
   plot(Coord(i, 1), Coord(i, 2), ...
      [Symbol(Z(i)) Color(Z(i))], 'MarkerSize', MarkerSize, 'LineWidth', LineWidth)
end
hold off
%saveas(1, '../../Articles/StatComput/Figures/FigNetworks-Indep-Symb', 'epsc')
%saveas(1, '../../Articles/StatComput/Figures/FigNetworks-Indep-Symb', 'fig')
%saveas(1, 'FigNetworks-Indep-Col', 'epsc')
%saveas(1, 'FigNetworks-Indep-Col', 'fig')
saveas(1, 'FigNetworks-Indep-Col', 'png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erdos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = zeros(2*C, 2*C);
Z = [ones(1, C) ones(1, C)]
Coord = zeros(2*C, 2);
p = 0.2
for c=1:2*C
   X(c, c+1:2*C) = (rand(1, 2*C-c) < p);
   Coord(c, :) = [R*cos(2*pi*c/2/C) R*sin(2*pi*c/2/C)];
   if sum(X(c, :))==0, X(c, c) = 1; end
end
X = X+X';
%figure(1), subplot(224)
gplot(X, Coord, 'k-'), axis off, hold on%, axis([-0.4 0.4 -0.4 0.4])
%saveas(4, 'FigNetworks-Erdos', 'epsc')
%saveas(4, 'FigNetworks-Erdos', 'fig')
saveas(4, 'FigNetworks-Erdos', 'png')
for i=1:length(Z)
   plot(Coord(i, 1), Coord(i, 2), ...
      [Symbol(Z(i)) Color(Z(i))], 'MarkerSize', MarkerSize, 'LineWidth', LineWidth)
end
hold off
%saveas(1, '../../Articles/StatComput/Figures/FigNetworks-Erdos-Symb', 'epsc')
%saveas(1, '../../Articles/StatComput/Figures/FigNetworks-Erdos-Symb', 'fig')
%saveas(1, 'FigNetworks-Erdos-Col', 'epsc')
%saveas(1, 'FigNetworks-Erdos-Col', 'fig')
saveas(1, 'FigNetworks-Erdos-Col', 'png')

% %saveas(1, 'FigNetworks.eps', 'epsc'), 
% %saveas(1, 'FigNetworks.fig', 'fig')
% %saveas(1, 'FigNetworks.fig', 'png')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphe de dépendance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Name = ['Z_1   '; 'X_{12}'; 'Z_2   '; 'X_{23}'; 'Z_3   '; 'X_{13}']
Angle = (2*acos(-1)) * (0:1/6:5/6)'
Coef = 0.03*[1 1.5 1 1.5 1 1.5]'
Coord = [cos(Angle) sin(Angle)]
Coord = Coord .* repmat(Coef, 1, 2)
X = [0 1 1 0 1 1; 1 0 1 0 0 0; 1 1 0 1 1 0; 
   0 0 1 0 1 0; 0 0 0 1 0 1; 1 0 0 0 1 0];
X = X + X';
gplot(X, Coord, 'ko'), axis off, hold on
for i=1:size(Name, 1)
   text(Coord(i, 1), Coord(i, 2)+.003, ...
      Name(i, :), 'Fontsize', Fontsize)
end
%saveas(1, 'FigNetworks-DepGraph', 'fig')
pause
gplot(X, Coord, 'k-o'), axis off, hold on
%plot(Coord(:, 1), Coord(:, 2), 'o')
hold off
%saveas(1, 'FigNetworks-DepGraph-Moral', 'fig')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pref. attachement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = [3 6 9] % effectif de chaque class
R = [0.1 0.25 0.4] % rayon de cercle
Pi = [0.8 0.6 0.4; 0.6 0.35 0.2; 0.4 0.1 0.1]
Q = length(N);
Z = [];
Coord = [];
for q=1:Q
   Z = [Z q*ones(1, N(q))];
   for n=1:N(q);
      Coord = [Coord; R(q)*cos(2*pi*(n+0.5)/N(q)) R(q)*sin(2*pi*(n+0.5)/N(q))];
   end      
end
X = zeros(sum(N), sum(N));
for i=1:sum(N);
   for j=i+1:sum(N)
      X(i, j) = (rand(1, 1) <= Pi(Z(i), Z(j)));
   end
end
X = X+X';
figure(1), %subplot(221)
gplot(X, Coord, 'k-'), axis off, hold on%, axis([-0.5 1.5 -0.4 0.4])
% saveas(1, 'FigNetworks-PrefAtt', 'epsc')
% saveas(1, 'FigNetworks-PrefAtt', 'fig')
for i=1:length(Z)
   plot(Coord(i, 1), Coord(i, 2), ...
      [Symbol(Z(i)) Color(Z(i))], 'MarkerSize', MarkerSize)
end
hold off
% saveas(1, 'FigNetworks-PrefAtt-Col', 'epsc')
% saveas(1, 'FigNetworks-PrefAtt-Col', 'fig')

return

