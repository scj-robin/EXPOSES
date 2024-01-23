% Raw CGH profile

clear all, clc, close all

Dir = '/media/donnees/RECHERCHE/RUPTURES/F-Picard/Programmes/matlab/data/'
File = 'Bt474.txt'
Data = load([Dir File]);
size(Data)

% Chrom 9
Data = Data((Data(:, 3)==9), :);
%Data(:, 2) = Data(:, 2)/1e6;
Data(:, 1) = Data(:, 1) - 1.35;
size(Data)

% Plot
plot(Data(:, 2), Data(:, 1), 'b.', 'MarkerSize', 15, 'LineWidth', 2), ...
   xlabel('genomic position'), ylabel('log_2 rat'), axis([1.565e6 1.675e6 -3 3])
#saveas(3, 'raw_profile_example.ps', 'ps')
