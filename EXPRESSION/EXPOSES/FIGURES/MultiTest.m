disp('Tests multiple : données de Golub');
disp('Expression de 1000 gènes, données d''apprentissage');
disp('A = ALL, B = AML');

disp('');
disp('Données');
%Fic = fopen('train-1000.prn', 'r');
Fic = fopen('train-7070.prn', 'r');
Titre = fgetl(Fic)
Data = fscanf(Fic, '%g');
fclose(Fic);
G = length(Data)/39
Data = reshape(Data, 39, G)';
NA = 27, NB = 11
Num = (1:G)';
A = Data(:, 2:28);
B = Data(:, 29:39);
alpha = 0.05

%disp(''); disp('Reduction'); G = 100, A = A(1:G, :); B = B(1:G, :); 
Liste = [(1:10) (G-9:G)]';
Debut = ceil(G/5)

disp('');
disp('Test');
MA = mean(A')'; VA = var(A')';
MB = mean(B')'; VB = var(B')';
figure %subplot(2, 3, 1), 
loglog(sqrt(VA), sqrt(VB), '.k', [1e1 1e4], [1e1 1e4], '-k'), ...
   %%title('std-dev B / A'), ...
   pause(0);
S1 = sqrt(((NA-1)*VA + (NB-1)*VB)/(NA+NB-2) * (1/NA + 1/NB));
S2 = sqrt(VA/NA + VB/NB);
T1 = (MA - MB)./S1;
T2 = (MA - MB)./S2;
ddlT = NA+NB-2
PT1 = zeros(G, 1); PT2 = PT1; PW = PT1;
for g=1:G;
%   PT1(g) = 2*(1-normcdf(abs(T1(g))));
%   PT2(g) = 2*(1-normcdf(abs(T2(g))));
   PT1(g) = 2*(1-Fdr_norm(abs(T1(g))));
   PT2(g) = 2*(1-Fdr_norm(abs(T2(g))));
%   PW(g) = ranksum(A(g, :), B(g, :));
%   if PW(g) < 0.5, 
%      PW(g) = 2*PW(g);
%   else
%      PW(g) = 2*(1-PW(g));
%   end
end
figure %subplot(2, 3, 4), 
plot(T1, T2, '.k', [1e-15 1e0], [1e-15 1e0], '-k', ...
   [1e-15 1e0], [0 0], '-.k', [0 0], [1e-15 1e0], '-.k'), ...
   %title('Welch / student (T)'), ...
   pause(0); 
figure %subplot(2, 3, 2), 
loglog(PT1, PT2, '.k', [1e-15 1e0], [1e-15 1e0], '-k', ...
   [1e-15 1e0], [alpha alpha], '-.k', [alpha alpha], [1e-15 1e0], '-.k'), ...
   %title('Welch / student (P='), ...
   pause(0);  
%%figure %subplot(2, 2, 2), loglog(PT1, PW, '.'), %title('Wilcoxon / student');
%%figure %subplot(2, 2, 3), loglog(PT2, PW, '.'), %title('Wilcoxon / Welch');
%T1T2 = crosstab(1+(PT1<alpha), 1+(PT2<alpha))
%T1W = crosstab(1+(PT1<alpha), 1+(PW<alpha))
%T2W = crosstab(1+(PT2<alpha), 1+(PW<alpha))

disp('');
disp('Palmares');
[X Ordre] = sort(-abs(T2));
T1 = T1(Ordre); PT1 = PT1(Ordre); 
T2 = T2(Ordre); PT2 = PT2(Ordre); 
PW = PW(Ordre);
disp([Liste T1(Liste) PT1(Liste) T2(Liste) PT2(Liste) PW(Liste)])

disp('');
disp('Ajustement multiples : seuils (Welch)');
Liste = (1:10)';
T = T2; P = PT2;
SBon = alpha/G*ones(G, 1);
SSidak = (1 - (1-alpha)^(1/G))*ones(G, 1);
SHolm = alpha./(G:-1:1)';
SSidakAdapt = 1 - ((1-alpha)).^(1./(G:-1:1)');
figure %subplot(2, 3, 5), 
semilogy([1 G], [alpha alpha], '-k', (1:G), P, '.k', ...
   (1:G), SBon, '-r', (1:G), SSidak, '-g', ...
   (1:G), SHolm, '-.r', (1:G), SSidakAdapt, '-.g'), ...
   %title('thresholds'), ...
pause(0);

disp('');
disp('Nombre de genes différentiellement exprimés');
disp(['            P   Bonferroni        Sidak         Holm    Sidak adapt.']);
disp([sum(P < alpha) sum(P < SBon) sum(P < SSidak) ...
      sum(P < SHolm) sum(P < SSidakAdapt)]);

%return
   
disp('');
disp('Probabilités critiques ajustées');
Bon = min([G*P ones(G, 1)]')';
Sidak = min([1-(1-P).^G ones(G, 1)]')';
Holm = min([(G:-1:1)'.*P ones(G, 1)]')';
SidakAdapt = min([1-(1-P).^((G:-1:1)') ones(G, 1)]')';
for g=2:G
   if Holm(g) < Holm(g-1), Holm(g) = Holm(g-1); end
   if SidakAdapt(g) < SidakAdapt(g-1),  SidakAdapt(g) = SidakAdapt(g-1);  end
end
FDR = min([G*P./((1:G)') ones(G, 1)]')';
FDR = flipud(FDR);
for g=2:G
   if FDR(g) > FDR(g-1), FDR(g) = FDR(g-1); end
end
FDR = flipud(FDR);
%Holm = min([Holm ones(G, 1)]')';
%SidakAdapt = min([SidakAdapt ones(G, 1)]')';
%FDR = min([FDR ones(G, 1)]')';
disp(['         Num             T            P   Bonferroni        Sidak         Holm    Sidak adapt.          FDR']);
disp([Liste T(Liste) P(Liste) Bon(Liste) Sidak(Liste) Holm(Liste) SidakAdapt(Liste) FDR(Liste)]);
figure %subplot(2, 3, 3), 
semilogy([1 G], [alpha alpha], '-k', (1:G), P, '.k', ...
   (1:G), Bon, '-r', (1:G), Sidak, '-g', (1:G), Holm, '-.r', (1:G), SidakAdapt, '-.g', ...
   (1:G), FDR, '-.b'), ...
   %title('adjusted p-values'), ...
   pause(0);
figure %subplot(2, 3, 6), 
semilogy([1 Debut], [alpha alpha], '-k', (1:Debut), P(1:Debut), '.k', ...
   (1:Debut), Bon(1:Debut), '-r', (1:Debut), Sidak(1:Debut), '-g', ...
   (1:Debut), Holm(1:Debut), '-.r', (1:Debut), SidakAdapt(1:Debut), '-.g', ...
   (1:Debut), FDR(1:Debut), '-.b'), ...
   %title('adjusted p-values (zoom)'), ...
   pause(0);

disp('');
disp('Nombre de genes différentiellement exprimés');
disp(['            P   Bonferroni        Sidak         Holm    Sidak adapt.          FDR']);
disp([sum(P < alpha) sum(Bon < alpha) sum(Sidak < alpha) ...
      sum(Holm < alpha) sum(SidakAdapt < alpha) sum(FDR < alpha) ]);


return

