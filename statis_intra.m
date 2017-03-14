function [ B, Wd, VAPU, VEPU, corrvars, V_pour ] = statis_intra( X, Wn, Wcomp, indnames, varetudes, varnames, p )
%% Fonction de calcul de l'intrastructure pour la methode STATIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input variables
% Wcomp = Matrice avec le Compromis entre les objets W
% X = Matrice des variables
% indnames = variable de type string qui a le nom des individus
% varnames = variable de type string qui a le nom des variables
% p = pourcentage de inertie minimum requis (pour choisir le nb d'axes)
%     par default on garde les 2 premiers axes si cette variable n'est pas
%     definie
%
% Output Variables
% B = Image euclidienne compromis des individus
% Corr = Correlation des variables avec les axes du compromis
%        Matrice N x M x T qui a les correlations o?:
%        - N nb des vecteurs propes
%        - M nb nombre des variables
%        - T nb d'?tudes
%
% Use:
% [ B ] = statis_intrastatis_intra( Wcomp,Xc,p,indnames, varnames)
%
% Author: Rodrigo Andres Rivera Martinez
% Corrections: Larbi Mouchou, Mounir Bendali-Braham, Nafise Gouard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Centrage et reduction des donn?es 
[L,C,n] = size(X);
for i=1:n
    Xc(:,:,i) = centrer(X(:,:,i),mean(mean(X(:,:,i))), std(std(X(:,:,i))));
end
%% Definition de l'image euclidienne
[L,C] = size(Wcomp);
D = 1/L*eye(L);
Wd = Wcomp*D;

[XU,VAPU, VEPU] = ACP(Wd);
n = size(VAPU,1);
V_pour = (VAPU*100)/sum(VAPU);
p_tot = sum(V_pour(1:2));
j=2;
if nargin < 3
    while p > p_tot
        j=j+1;
        if j>=n; break; end;
        p_tot = p_tot + V_pour(j);
    end
end

B = XU(:,1:j); 

disp('Valeur propres');
disp(VAPU);
fprintf('Pourcentage de l''inertie cumule: %.3f %%\n',p_tot);
fprintf('Nb d''axes: %d\n',j);

%% Correlations des variables avec les axes du compromis
[L1,C1,n] = size(Xc);

for l=1:j % axes
    for i = 1:C1 % variables
        for k = 1:n % etudes
            corrvars(k,i,l)=Xc(:,i,k)'*D*VEPU(:,l);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot de l'image euclidienne compromis des individus
figure;
scatter(B(:,1),B(:,2)); grid on;
xlabel(sprintf('Axe 1 (Inertie: %.2f %%)',V_pour(1)));
ylabel(sprintf('Axe 2 (Inertie: %.2f %%)',V_pour(2)));
title('Image euclidienne compromis des individus')

if nargin <=4
    for i=1:L
        indnames{i} = sprintf('Individu %d',i);
    end
end

for i=1:L
    text(B(i,1), B(i,2),indnames(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot des correlations des variables
[nb_inds, nb_vars, nb_etudes] = size(X);
figure; hold on; grid on;
for var= 1:nb_vars
  plot(corrvars(:,var,1), corrvars(:,var,2), '-O');
 for t = 1: nb_etudes
     text(corrvars(t,var,1), corrvars(t,var,2), [varnames(var) num2str(t)]);
 end
end

xlabel(sprintf('Axe 1 (Inertie: %.2f %%)',V_pour(1)));
ylabel(sprintf('Axe 2 (Inertie: %.2f %%)',V_pour(2)));
title('Corr?lations des variables')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end



function [XU, VAPU, VEPU] = ACP(X)
%--------------------------------
% Calcul ACP
%--------------------------------
% Recherche des valeurs et vecteurs propres
[VEPU, VAPU] = eig(X'*X);    
VAPU         = diag(VAPU);        
%
% Ordonnancement des valeurs et vecteurs propres
[VAPU,s] = sort(VAPU);
VAPU     = VAPU(flipud(s)); 
VEPU     = VEPU(:,flipud(s)); 
%
% Nouvelles Coordonn?es (Composantes principales)
XU = X * VEPU; 
end

function [Ac] = centrer(A,mean_A,std_A)
%--------------------------------
% Centrage des donn?es
%--------------------------------
UN = ones(size(A));
Me = UN * mean_A;
Ecart_type = UN * diag(std_A);
Ac  = (A - Me)./Ecart_type;
end