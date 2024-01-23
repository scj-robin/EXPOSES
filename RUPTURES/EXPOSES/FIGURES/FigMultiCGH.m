% Figure exposé multi CGH
clear all, clc

% Rep
DataDir = 'D:\RECHERCHE\RUPTURES\ETUDES\Nakao\'
DataFile = 'nakao-mat.txt'

% Data
Data = load([DataDir DataFile]);
Chrom = Data(:, 1);
Position = Data(:, 2);
X = Data(:, 3:end)';
[m n] = size(X)
% plot(X, '.');

% Chrom
NbChrom = max(Chrom)
Length = zeros(1, NbChrom);
for c=1:NbChrom,    
   Length(c) = sum((Chrom==c));
end
Length

% % General plot
C = 20
NbCol = 10;    NbRow = ceil(m/NbCol);
Xmin = floor(min(min(X(:, (Chrom==C)))));
Xmax = ceil(max(max(X(:, (Chrom==C)))));
%NbRow = floor(sqrt(m));    NbCol = ceil(m/NbRow);
% for i=1:m
%    subplot(NbRow, NbCol, i), ...
%       plot(X(i, (Chrom==C)), '.'), ...
%       axis([0 Length(C) Xmin Xmax]);
% end

% Slide plot
figure(2)
Patient = [1 8 14 24 47 65 103]
Break = [1 16 28 97; 1 12 40 97; 1 97 Inf Inf; 1 40 97 Inf; 1 11 97 Inf; 1 40 97 Inf; 1 35 41 97]
PosMark = [11 43 58]
NbPat = length(Patient);
for p=1:NbPat, 
   for b = 1:size(Break, 2)-1
      if Break(p, b+1) <= 97
         Xtmp = X(Patient(p), (Chrom==C));
         Xtmp = Xtmp(Break(p, b):Break(p, b+1));
         Ltmp = length(Xtmp);
         mu = mean(Xtmp((isnan(Xtmp)==0)));
         subplot(NbPat, 1, p), ...
            plot((Break(p, b):Break(p, b+1)), mu*ones(1, Ltmp), 'g-', 'LineWidth', 2), ...
            axis([0 Length(C) Xmin Xmax]), ...
            hold on
      end
   end
   for n=1:length(PosMark)
      subplot(NbPat, 1, p), ...
         plot([PosMark(n)-0.5 PosMark(n)-0.5], [Xmin Xmax], '-r', ...
         [PosMark(n)+0.5 PosMark(n)+0.5], [Xmin Xmax], '-r', 'LineWidth', 1), ...
   end
   subplot(NbPat, 1, p), ...
   plot((1:Length(C)), X(Patient(p), (Chrom==C)), '.', 'MarkerSize', 10);
   axis off,        
   hold off;
end
saveas(2, [DataFile '-MixSeg.eps'], 'epsc')

% Slide plot -V2
figure(3)
Patient = [1 8 14 24 47 65 103]
Break = [1 16 28 97; 1 12 40 97; 1 97 Inf Inf; 1 40 97 Inf; 1 11 97 Inf; 1 40 97 Inf; 1 35 41 97]
PosMark = [43]
NbPat = length(Patient);
for p=1:NbPat, 
   for b = 1:size(Break, 2)-1
      if Break(p, b+1) <= 97
         Xtmp = X(Patient(p), (Chrom==C));
         Xtmp = Xtmp(Break(p, b):Break(p, b+1));
         Ltmp = length(Xtmp);
         mu = mean(Xtmp((isnan(Xtmp)==0)));
         subplot(NbPat, 1, p), ...
            plot((Break(p, b):Break(p, b+1)), mu*ones(1, Ltmp), 'c-', 'LineWidth', 2), ...
            axis([0 Length(C) Xmin Xmax]), ...
            hold on
        if b>1
            plot([Break(p, b) Break(p, b)], [Xmin Xmax], 'b-.', 'LineWidth', 2)
        end        
      end
   end
   for n=1:length(PosMark)
      subplot(NbPat, 1, p), ...
         plot([PosMark(n) PosMark(n)], [Xmin Xmax], '-r', 'LineWidth', 2), ...
   end
   subplot(NbPat, 1, p), ...
       plot((1:Length(C)), X(Patient(p), (Chrom==C)), 'k.', 'MarkerSize', 8);
   Xtmp = X(Patient(p), (Chrom==C));
   plot(43, Xtmp(43), 'k.', 'MarkerSize', 12);
   axis off,        
   hold off;
end
saveas(3, [DataFile '-MixSeg-V2.eps'], 'epsc')