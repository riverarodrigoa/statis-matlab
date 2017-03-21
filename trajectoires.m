function trajectoires( X, Wn, D, VEPU, VAPU, V_pour, indnames )
%% Fonction de calcul de de l'interstructure pour la methode STATIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input variables
% X = Tableaux avec les t ?tudes
% Wn = Matrice avec les objets des t etudes normes
% VEPU = Vecteurs propres de WD
% VAPU = Valeurs propres de WD
% V_pour = Pourcentages d'inertie
% indnames = noms des individus
%
% Use:
% trajectoires ( X, Wn, Wcomp, indnames, varnames, p )
%
% Authors: Larbi Mouchou, Mounir Bendali-Braham, Nafise Gouard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin<6)
    error('SYNTAXE ERROR: trajectoires( X, Wn, VEPU, VAPU, V_pour, indnames )');
end

%% Corr?lations individus pour le trac?e des trajectoires
[L,C,n] = size(X);

nb_etudes = n;
nb_inds = L;

for etude= 1:nb_etudes
    WU(etude, :, :) = (Wn(:,:,etude) *D* VEPU) * diag(1./sqrt(VAPU)); 
end;
figure;
hold on;
for ind= 1:nb_inds
    plot(WU(:,ind,1), WU(:,ind,2), '-O');
    for t = 1: nb_etudes
        text(WU(t,ind,1), WU(t,ind,2), [indnames(ind) 'An' num2str(t)]);
    end;
end;

grid on;

xlabel(sprintf('Axe 1 (Inertie: %.2f %%)',V_pour(1)));
ylabel(sprintf('Axe 2 (Inertie: %.2f %%)',V_pour(2)));
title('Corr?lations des variables')
end