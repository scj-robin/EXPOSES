function F = Fdr_norm(x)

% fonction de répartition de la loi normale N(0,1)

if x>0, F = 0.5*(1+erf(x/sqrt(2)));
else F = 0.5*(1-erf(-x/sqrt(2)));
end;
