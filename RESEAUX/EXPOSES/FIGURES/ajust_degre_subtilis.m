% ajustement de la loi de Zipf
%Graphe = load('Graphe.mat') 
%Degre = Graphe.K 
Graphe=textread('N:\Recherche\Réseaux\Graphes\Données_enzyme\degreS.txt','%u');
%hist(Graphe,300);
%pause
mx=max(Graphe);

% estimation chi2 minimum
rho=1.5

m=mean(Graphe((Graphe>1)));
p=[1:1000].^(-rho);

%troncature
t=1
% donnees utilisées
n=sum(Graphe>t)

zeta=sum(p(t+1:1000));
p=p(t+1:1000)./zeta;
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

%en(1:mx-t)=e(1:mx-t);en(23)=sum(e(23:1000-t));
th=e(1:mx-t);th(mx-t)=sum(e(mx-t:1000-t));
emp=freq(t+1:mx);


chi2=sum((th-emp).^2./th)
ddl=size(th,2)-2

f=freq./(sum(freq));
x=log([1:size(freq,2)]);l=log(f+0.001);
plot(x,l,'*')
return

