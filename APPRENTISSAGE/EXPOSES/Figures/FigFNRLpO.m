clear all, clc

% Param
s = [5 10 25 50] % coef beta H1
m = 5 % nb method
Color = 'rbgck'

% Data
cd('D:\RECHERCHE\APPRENTISSAGE\EXPOSES\Figures')
Tab = load('CeR10-JSPI-Tab5.txt')
% s p0	LPO LOO Twil BH Best
S = Tab(:, 1);
p0 = Tab(:, 2);
FNR = Tab(:, 3:end)/100;
FNR = FNR - repmat(FNR(:, m), 1, m)
m = m-1

% plot
for i = 1:4
    subplot(2, 2, i)
    [p0(S==s(i)) FNR((S==s(i)), :)]
    for j=1:m
        plot(p0(S==s(i)), FNR((S==s(i)), j), [Color(j) '-'], ...
            'LineWidth', 2), hold on
    end
    hold off
end

    	
