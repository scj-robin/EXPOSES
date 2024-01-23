clear all, clc, fclose all;

RepFig = 'D:\RECHERCHE\RESEAUX\Exposes\FIGURES\';
% Motif
%E = [0 1 ; 0 0]; Motif = 'I';
E = [0 1 1 ; 0 0 0; 0 0 0]; Motif = 'V';
%E = [0 1 1 ; 0 0 1; 0 0 0]; Motif = 'Triangle';
%E = [0 1 0 1; 0 0 1 0; 0 0 0 1; 0 0 0 0]; Motif = 'Square';
%E = [0 1 0 0; 0 0 1 0; 0 0 0 1; 0 0 0 0]; Motif = 'Chain4';
%E = [0 1 1 1; 0 0 1 1; 0 0 0 1; 0 0 0 0]; Motif = 'Clique4';
%E = [0 1 1 1; 0 0 1 0; 0 0 0 0; 0 0 0 0]; Motif = 'QT';
%E = [0 1 1 1; 0 0 0 0; 0 0 0 0; 0 0 0 0]; Motif = 'Star3';
%E = [0 1 1 1 1; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]; Motif = 'Star4';
%E = [0 1 1 1 1 1; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]; Motif = 'Star5';
% Symétrisation
E = min(1, E + E')
k = size(E, 1)

figure(1), axis off

% Coordonnées
if length(strfind(Motif, 'Star')) > 0
   Coord = [0 0];
   Theta = (0.5:k-0.5)'*2*acos(-1)/(k-1);
   Coord = [Coord; cos(Theta) sin(Theta)];
else
   Theta = (1:k)'*2*acos(-1)/(k);
   Coord = [cos(Theta) sin(Theta)];
end

% Figure et exportation
gplot(E, Coord, 'k-o'), axis off
%saveas(1, [RepFig 'Motif-' Motif '.fig'], 'fig')
%saveas(1, [RepFig 'Motif-' Motif '.eps'], 'eps')
