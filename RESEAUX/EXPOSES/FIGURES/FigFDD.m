clear all, clc;

RepFig = 'D:\RECHERCHE\RESEAUX\Exposes\FIGURES\';

% FDD : rewiring
figure(1)
plot([0 0], [0 1], '.-k', 'MarkerSize', 100, 'LineWidth', 6)
hold on, 
plot([1 1], [0 1], '.-k', 'MarkerSize', 100, 'LineWidth', 6)
text(-0.2, 0, 'i_1', 'FontSize', 40)
text(-0.2, 1, 'j_1', 'FontSize', 40)
text(0.8, 0, 'i_2', 'FontSize', 40)
text(0.8, 1, 'j_2', 'FontSize', 40)
axis([-0.3, 1.1, -0.2, 1.2]), axis off
hold off
saveas(1, [RepFig 'FigFDD-Rewire-1.fig'], 'fig')
saveas(1, [RepFig 'FigFDD-Rewire-1.eps'], 'eps')

figure(2)
plot([0 1], [0 1], 'k.-', 'MarkerSize', 100, 'LineWidth', 6)
hold on, 
plot([1 0], [0 1], 'k.-', 'MarkerSize', 100, 'LineWidth', 6)
text(-0.2, 0, 'i_1', 'FontSize', 40)
text(-0.2, 1, 'j_1', 'FontSize', 40)
text(0.8, 0, 'i_2', 'FontSize', 40)
text(0.8, 1, 'j_2', 'FontSize', 40)
axis([-0.3, 1.1, -0.2, 1.2]), axis off
hold off
saveas(2, [RepFig 'FigFDD-Rewire-2.fig'], 'fig')
saveas(2, [RepFig 'FigFDD-Rewire-2.eps'], 'eps')

% FDD : semi-edges
k = [0 2 1 4 1]
n = 5
m = sum(k)/2
r = 0.4

figure(3)
for i=1:n
    theta = 2*(i+0.5)/n
    plot(cos(theta*pi), sin(theta*pi), 'k.', 'MarkerSize', 50), 
    axis([-1.2, 1.2, -1.2, 1.2]), axis off
    hold on;
    k(i)
    phi = (-(k(i)+1)/2 +(1:k(i))) / (k(i))
    for j=1:k(i)
        phi(j) = theta + (n-2)/n * phi(j);
        plot([cos(theta*pi) cos(theta*pi)-r*cos(phi(j)*pi)], ...
           [sin(theta*pi) sin(theta*pi)-r*sin(phi(j)*pi)], ...
           'k-.', 'LineWidth', 3),         
    end    
end
hold off
saveas(3, [RepFig 'FigFDD-SemiEdge.fig'], 'fig')
saveas(3, [RepFig 'FigFDD-SemiEdge.eps'], 'eps')
