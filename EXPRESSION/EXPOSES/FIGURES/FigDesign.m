clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plan a 2 ESPCI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Individus
XX = [1 1; 1 -1; -1 1; -1 -1];
R = 5
X = [];
for x = 1:size(XX, 1);
   for r = 1:R
      X = [X; XX(x, :)];
   end
end
X = [X X(:, 1).*X(:, 2)];
n = size(X, 1)

% plan proposé
P = zeros(n/2, n/2);
for i=1:n/2
   % Plan ESPCI
   P(i, i) = 1;
   P(i, rem(i, n/2)+1) = -1;
   P(i, rem(i+R-1, n/2)+1) = 1;
   P(i, rem(i+R, n/2)+1) = -1;
   
   % Plan complet
   for j=1:n/2
      P(i, j) = (-1)^(i+j);
   end
end
%P, pause

% matrice des hybridations
[I, J] = find(P ~=0);
[I Rg] = sort(I);
J = J(Rg);
D = zeros(length(I), n);
for i=1:length(I)
   if P(I(i), J(i)) == 1
      D(i, I(i)) = 1;
      D(i, n/2+J(i)) = -1;
   else
      D(i, I(i)) = -1;
      D(i, n/2+J(i)) = 1;
   end
end
%D

% matrices de variances
s2 = 1, g2 = 2
DX = D*X;
VDX = 2*s2*eye(size(D, 1)) + 2*g2*D*D';
IVDX = inv(VDX);
VarT = inv(DX'*IVDX*DX)
DetV = det(VarT)
%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plan a 2 facteurs : modèle mixte
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%orient landscape, 
figure(1)
gamma2 = 10, sigma2 = 1
x = (-2:0.01:2);
Ymax = 1.1/sqrt(2*pi*gamma2);
Coef = sqrt(sigma2/gamma2)/3;
FontSize = 24

plot(x*sqrt(gamma2), exp(-x.^2)/sqrt(2*pi*gamma2), 'b-', ...
   'LineWidth', 3), hold on
text(0, 0.97*Ymax, ['\mu'], 'FontSize', FontSize), 
plot([0 0], [0 Ymax], 'b:', 'LineWidth', 2),


A = sqrt(gamma2)*[1.1 -0.5 0.75]
I = length(A)
for i=1:I
   %A = sqrt(gamma2) * randn(1)
   plot(A(i)+x*sqrt(sigma2), Coef*exp(-x.^2)/sqrt(2*pi*sigma2), 'r-', ...
      'LineWidth', 3), 
   text(A(i)-0.75, (0.3+0.1*i)*Ymax, ['\mu+A_' num2str(i)], 'FontSize', FontSize), 
   plot([A(i) A(i)], [0 0.8*Ymax], 'r:', 'LineWidth', 2),
end
plot([A(2)-0.8 A(2)-0.8], [0 Ymax/4], 'k:', 'LineWidth', 2),
h = text(A(2)-2.75, 0.1*Ymax, ['\mu+A_2+E_{2r}'], 'FontSize', FontSize), 
hold off

saveas(1, 'FigMixedModel.eps', 'epsc')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plan a 2 facteurs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma2 = 5, sigma2 = 1

% Matrice des traitements
X = [1 1; 1 1; 1 -1; 1 -1; -1 1; -1 1; -1 -1; -1 -1];
X = [X X(:, 1).*X(:,2)]
VX = gamma2*eye(size(X, 1));

% Matrice des hybridations et des swaps
D = [1 0 0 0 -1 0 0 0; 0 1 0 0 0 0 -1 0; 0 0 1 0 0 -1 0 0; 0 0 0 1 0 0 0 -1];
DD = D;
D = []; SD = [];
for d=1:size(DD, 1)
   D = [D; DD(d, :); -DD(d, :)];
   SD = [SD; DD(d, :)+DD(d, :)];
end
D, DX = D*X
VDX = D*VX*D' + 2*sigma2*eye(size(D, 1))
SD, SDX = SD*X, 
VSDX = SD*VX*SD' + 4*sigma2*eye(size(SD, 1))
IVSDX = inv(VSDX)

% Estimateur
T = inv(SDX'*IVSDX*SDX)*SDX'*IVSDX
VarT = inv(SDX'*IVSDX*SDX)
%T = inv(SX'*SX)*SX'
%VarT = inv(SX'*SX)

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù
% plan 2^3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù
X = []
for i=-1:2:1
   for j=-1:2:1
      for k=-1:2:1
         X = [X; i j k];
      end
   end
end
X = [ones(8, 1) X X(:,1).*X(:,2) X(:,1).*X(:,3) X(:,2).*X(:,3) X(:,1).*X(:,2).*X(:,3)]

