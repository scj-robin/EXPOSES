% Figure segmentation CGH
clear all, clc

LineWidth = 2
FontSize = 14

% Segment
figure(1), orient portrait
plot([0, 10], [0, 0], 'k-', 'LineWidth', LineWidth), ...
   axis([-0.2, 10.2, -0.2, 0.5]), axis off, hold on

% Rupture
T = [0 1.5 4.5 8 10]
mu = [0.05 0.1 0.05 0.07]
for t=1:length(T)
   plot([T(t) T(t)], [-0.05 0.12], 'b-', 'LineWidth', LineWidth)
   text(T(t)-0.05, 0.13, ['t_' num2str(t-1)], ...
      'FontSize', FontSize, 'Color', 'blue')
   if t>1
      text((T(t)+T(t-1))/2-0.1, -0.035, ['I_' num2str(t-1)], ...
         'FontSize', FontSize, 'Color', 'red')
      plot([T(t-1) T(t)], [mu(t-1) mu(t-1)], 'r-','LineWidth', LineWidth) 
      text((T(t)+T(t-1))/2-0.5, mu(t-1)+0.025, ['\lambda_t = \mu_' num2str(t-1)], ...
         'FontSize', FontSize, 'Color', 'red')
   end
end
text(0.05, 0.02, 'Chromosome', 'FontSize', FontSize)
hold off
saveas(1, 'FigSeg_Intro.fig', 'fig')
saveas(1, 'FigSeg_Intro.eps', 'epsc')

% 2ème version (sans lambda, avec sigma)
figure(2), orient portrait
plot([0, 10], [0, 0], 'k-', 'LineWidth', LineWidth), ...
   axis([-0.2, 10.2, -0.2, 0.5]), axis off, hold on
for t=1:length(T)
   plot([T(t) T(t)], [-0.05 0.12], 'b-', 'LineWidth', LineWidth)
   text(T(t)-0.05, 0.13, ['t_' num2str(t-1)], ...
      'FontSize', FontSize, 'Color', 'blue')
   if t>1
      text((T(t)+T(t-1))/2-0.1, -0.035, ['I_' num2str(t-1)], ...
         'FontSize', FontSize, 'Color', 'red')
      plot([T(t-1) T(t)], [mu(t-1) mu(t-1)], 'r-','LineWidth', LineWidth) 
      text((T(t)+T(t-1))/2-0.5, mu(t-1)+0.025, ['\mu_' num2str(t-1) ', \sigma_{(' num2str(t-1) ')}'], ...
         'FontSize', FontSize, 'Color', 'red')
   end
end
text(0.05, 0.02, 'Chromosome', 'FontSize', FontSize)
hold off
saveas(2, 'FigSeg_Intro_bis.fig', 'fig')
saveas(2, 'FigSeg_Intro_bis.eps', 'epsc')
