function [ B ] = statis_intra( Wcomp, varnames, p )
%% Fonction de calcul de l'intrastructure pour la methode STATIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input variables
% Wcomp = Matrice avec le Compromis entre les objets W
%
% Output Variables
% B = Image euclidienne compromis des individus
% varnames = variable de type string qui a le nom des variables
% p = pourcentage de inertie minimum requis (pour choisir le nb d'axes)
%     par default on garde les 2 premiers axes si cette variable n'est pas
%     definie
% Use:
% [ B ] = statis_intra( Wcomp, varnames )
%
% Autor: Rodrigo Andres Rivera Martinez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Definition de l'image euclidienne
[L,C] = size(Wcomp);
D = 1/L*eye(L);
Wd = Wcomp*D;

[XU,VAPU, VEPU] = ACP(Wd);
n = size(VAPU,1);
V_pour = (VAPU*100)/sum(VAPU);
p_tot = sum(V_pour(1:2));
j=2;
if nargin > 2
    while p > p_tot
        j=j+1;
        if j>n break; end;
        p_tot = p_tot + V_pour(j);
    end
end

 B = XU(:,1:j); 

disp('Valeur propres');
disp(VAPU);


fprintf('Pourcentage de l''inertie: %3f4\n',p_tot);
fprintf('Nb d''axes: %d\n',j);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot de l'image euclidienne compromis des individus
figure;
scatter(B(:,1),B(:,2)); grid on; xlabel('Axe 1'); ylabel('Axe 2');
title('Image euclidienne compromis des individus')

if nargin <2
    for i=1:L
        varnames{i} = sprintf('Individu %d',i);
    end
end

for i=1:L
    text(B(i,1), B(i,2),varnames(i));
end

end
function [XU,VAPU, VEPU] = ACP(X)
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