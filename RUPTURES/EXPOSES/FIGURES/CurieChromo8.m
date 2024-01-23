clear all, clc

%Rep = '/RECHERCHE/RUPTURES/E-Budinska/FINAL RESULTS/analyses by EM/Chromosome 8/balanced/EM/Franck 0.5 outer/'
Rep = '/RECHERCHE/RUPTURES/E-Budinska/FINAL RESULTS/analyses by EM/Chromosome 8/unbalanced/results EM/Emily/';

% Lecture des données
%Data = load([Rep 'rupt.CGH.chr.8.dat']);
%Data = load([Rep 'rupt.CGH.chr.8.initial.noNA.dat']);
Data = load([Rep 'rupt.CGH.chr.8.initial.NA.Emily.dat']);
Group = Data(:, 1); Patient = Data(:, 2);
Begin = Data(:, 3); End = Data(:, 4);

% Suppression des extémités triviales 1 et T
G = max(Group); T = max(End);
Begin(Begin==1) = NaN;
End(End==T) = NaN;

for g=1:G
   % Sélection du group
   BeginTmp = Begin(Group==g);
   PatientTmp = Patient(Group==g);
   
   % Histo du nombre de ruptures
   figure(1)
   Pg = max(PatientTmp)
   K = zeros(1, Pg);
   for p=1:Pg, K(p) = sum(PatientTmp == p)-1; end
   [Kmax pmax] = max(K)
   hist(K, (0:15));
   saveas(1, ['/RECHERCHE/RUPTURES/Exposes/Figures/CurieChromo8Gp' num2str(g) '_K.eps'], 'epsc')
   
   % Histo de des positions des rupture
   figure(2)
   hist(BeginTmp, (1:T)), axis([1 T 0 10]);;
   saveas(2, ['/RECHERCHE/RUPTURES/Exposes/Figures/CurieChromo8Gp' num2str(g) '_T.eps'], 'epsc')
end
   