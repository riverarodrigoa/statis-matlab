function [ Wcomp ] = compromis(W,S,Delta,VaP,VeP,norm)
%% Fonction de calcul de compromis pour la methode STATIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input variables
% W = Matrice avec les objets des t etudes
% S = Matrice des produits scalaire entre les tableaux
% Delta = Matrice diagonal avec les poids
% norm = Option si on veut faire l'analyse en prendre en compre la norme
% VaP = Valeurs propres du matrice S
% VeP = Vecteurs propres du matrice S
%
% Output Variables
% Wcomp = Matrice avec le Compromis entre les objets W
%
% Use:
% [ Wcomp ] = compromis(W,S,Delta,VaP,VeP,norm)
%
% Autor: Rodrigo Andres Rivera Martinez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calcul de l'expresion definitive du compromis
[L,C,T] = size(W);
pi_t=diag(Delta);
SS = diag(S);
gamma = VeP(:,1);
Wcomp = 0;

if ~norm
    alpha_c = sum(pi_t.*sqrt(SS));
    for i = 1:T
      Wcomp = Wcomp + ((1/sqrt(VaP(i))*alpha_c*pi_t(i)*gamma(i))*W(:,:,i));
    end
else
    for i = 1:T
        Wcomp = Wcomp + ((1/sqrt(VaP(i))*pi_t(i)*gamma(i))*W(:,:,i));
    end
end

end

