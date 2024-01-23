Data = load('LoessLamesJaunes.txt');

M = Data(:, 1);
A = Data(:, 2)/2;
Loess = Data(:, 3);
Residu = Data(:, 4);

figure(1)
plot(A, Residu, '.g', [min(A) max(A)], [0 0], 'b-.', ...
   A, M, '.k', A, Loess, 'r-', 'MarkerSize', 10, 'LineWidth', 2)

saveas(1, 'LoessLamesJaunes.fig', 'fig')
orient landscape
saveas(1, 'LoessLamesJaunes.eps', 'epsc')
