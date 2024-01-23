% ajustement du mélange de lois de poisson
Graphe=textread('N:\Recherche\Réseaux\Graphes\Données_enzyme\degreS.txt','%u');

mx=max(Graphe);

% SR : Graphe = load('Graphe');

%Degre = Graphe.K ;
%hist(Degre,35)
%pause

a1=0.003;a2=0.016;a3=0.039;a4=0.058;a5=0.176;a6=0.22;a7=1-a1-a2-a3-a4-a5-a6;
l1=99.6;l2=70.3;l3=51.5;l4=26.1;l5=15.3;l6=8.6;l7=3.2;
for i=1:150
p(i)=(a1*exp(-l1)*(l1^i)+a2*exp(-l2)*(l2^i)+a3*exp(-l3)*(l3^i)+...
    a4*exp(-l4)*(l4^i)+a5*exp(-l5)*(l5^i)+a6*exp(-l6)*(l6^i)+...
    a7*exp(-l7)*(l7^i))/factorial(i);
end;
m_th=[1:150]*p'
m_emp=mean(Graphe)

%troncature
t=0
% donnees utilisées
n=sum(Graphe>t)

zeta=sum(p(t+1:150));
p=p(t+1:150)./zeta;
cdft=cumsum(p);
cdft(mx-t)=1;
cdft=cdft(1:mx-t);

plot(cdft(1:mx-t));
pause;

e=sum(Graphe>t)*p;

for i=1:mx
freq(i)=sum(Graphe==i);
end

cft=cumsum(freq(t+1:mx))/(sum(freq(t+1:mx)));

plot([t+1:mx],e(1:mx-t),'+',[t+1:mx], freq(t+1:mx),'*');
pause;

plot(cdft,cdft,'-',cdft,cft,'+');
pause;


th=e(1:mx-t);th(mx-t)=sum(e(mx-t:150-t));
emp=freq(t+1:mx);
th(81)=sum(th(81:mx-t));emp(81)=sum(emp(81:mx-t));
th=th(1:81);emp=emp(1:81);
chi2=sum((th-emp).^2./th)
ddl=size(th,2)-13


return


