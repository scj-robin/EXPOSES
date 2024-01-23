clear all, clc

MarkerSize = 5
LineWidth = 2
Color = 'krbgmc'
NbColor = length(Color)

% Lecture des données
Rep = '/RECHERCHE/RUPTURES/CLIMATO/DATA/'
Temp = load([Rep 'temp.min.annuelles.dat']);
Temp = Temp(1:end-2, :);
Time = Temp(:, 1);
Temp = Temp(:, 2:end);
Temp = Temp/10;
[T S] = size(Temp)

% Lecture des stations
FicStat = fopen([Rep 'location.dat'])
Station = fgetl(FicStat)
LgStat = length(Station)
size(Station)
for s=2:S
   Line = fgetl(FicStat)
   if length(Line) < LgStat, 
      Line = [Line ' '*ones(1, LgStat-length(Line))];
   else
      Station = [Station ' '*ones(s-1, length(Line)-LgStat)];
      size(Station)
      LgStat = length(Line)
   end
   Station = [Station; Line];
end
Station
   
% Suppression des stations avec données manquantes
NbMiss = sum(isnan(Temp));
Temp = Temp(:, (NbMiss==0));
Station = Station((NbMiss==0), :)
[T S] = size(Temp)
TempMax = ceil(max(max(Temp)))
TempMin = floor(min(min(Temp)))

% Graphe par station
% for s=1:S
%    plot(Time, Temp(:, s), 'o', 'MarkerSize', MarkerSize), ...
%       axis([Time(1) Time(T) TempMin TempMax]); 
%    pause
% end

% Régression par période
Tau = [0 31 47]
%Tau = [0 6 47]
K = length(Tau)-1
a = zeros(K, S);
b = zeros(K, 1); sigma2 = b;
for k=1:K
   [a(k, :), b(k), sigma2(k)] = F_RegSimul(Temp(Tau(k)+1:Tau(k+1), :), ...
      Time(Tau(k)+1:Tau(k+1)));
end
a, b

% Graphique
figure(1)
Col = 0;
ListStation = [1 2 5 6 10 24]
for s=ListStation
   %if rand(1) <= NbColor/S
      Col = Col+1;
      Color(mod(Col, NbColor)+1)
      plot(Time, Temp(:, s), [Color(mod(Col, NbColor)+1) 'o'], 'MarkerSize', MarkerSize), ...
         axis([Time(1) Time(T) TempMin TempMax]); ...
         hold on
      for k=1:K
         plot(Time(Tau(k)+1:Tau(k+1)), a(k, s) + b(k)*Time(Tau(k)+1:Tau(k+1)), ...
            [Color(mod(Col, NbColor)+1) '-'], 'LineWidth', LineWidth)
      end
      text(Time(round((Tau(2)+Tau(3))/2)), mean(Temp(1:Tau(2), s))+0.2, deblank(Station(s, :)))
   %end
end
plot([Time(Tau(2))+0.5 Time(Tau(2))+0.5], [TempMin TempMax], 'k-.', ...
   'LineWidth', LineWidth)
hold off
%saveas(1, '/RECHERCHE/RUPTURES/Exposes/Figures/SeglinregFloraison.eps', 'epsc')

mu = mean(mean(a))
sigma_U = std(mean(a)-mu)
b
sigma_E = sqrt(mean(sigma2))