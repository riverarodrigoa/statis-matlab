function [Co,S,SS,RV] = statis_interstructure (X,M,Delta,Sup,norm,varnames)
%% Fonction de calcul de de l'interstructure pour la methode STATIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input variables
% X = Tableaux avec les t ?tudes
% M = Metrique (usuelment matrice identit?)
% Delta = Matrice diagonal avec les poids
% Sup = Matrice avec les tableaux supplementaires
% norm = Option si on veut faire l'analyse en prendre en compre la norme
%
% PARAMETRES OPCIONELS
% varnames = variable de type string qui a le nom des variables
%
% Output Variables
% Co = Matrice avec les composantes principales
% S = Matrice des produits scalaire entre les tableaux
% SS = Matrice des produits scalaire entre les tableaux afect? par le poids
%      Delta
% RV = Matrice avec les coefficients RV
%
% Use:
% [Co,S] = statis_interstructure (X,M,Delta,Sup,norm)
%
% Autor: Rodrigo Andres Rivera Martinez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% D?finition des objets repr?sentatifs
% Verification des matrices
if(nargin<5)
    error('SYNTAXE ERROR: statis_interstructure (X,M,Delta,Sup,norm)');
end

[L,C,n] = size(X);
[L1,C1]=size(M);
[L2,C2]=size(Delta);
[L3,C3]=size(Sup);

if L~=L3
    error('ERROR: Nb d''individus doit etre egal dans les matrices X et Sup');
end
if C ~= 3
    error('ERROR: Nb de variables doit etre egal dans les matrices X et Sup');
end
if size(varnames,2) < C || size(varnames,2) > C
    error('ERROR: <varnames> doit etre de la meme taille que le nb de variables');
end
%-------------------------------------------------------------------------------
% Definition des objets
%-------------------------------------------------------------------------------

for i = 1:n
     W(:,:,i) = X(:,:,i)*X(:,:,i)';
end
%-------------------------------------------------------------------------------
% Calcul de la matrice des produits scalaires S
%-------------------------------------------------------------------------------
if ~norm
    for i=1:n
        for j=1:n
            S(i,j)=prod_hs(W(:,:,i),W(:,:,j));
        end
    end
else
    for i=1:n
        for j=1:n
        S(i,j)=(prod_hs(W(:,:,i),W(:,:,j)))/(norme(W(:,:,i))*(norme(W(:,:,j))));
        end
    end
end
%-------------------------------------------------------------------------------
% Calcul de coefficient RV
%-------------------------------------------------------------------------------
for i=1:n
    for j=1:n
        RV(i,j) = S(i,j)/(sqrt(S(i,i))*(sqrt(S(j,j))));
    end
end
%-------------------------------------------------------------------------------
% Image euclidienne des objets
%-------------------------------------------------------------------------------
SS = S*Delta;

Cp = ACP(SS);
% Par le th?oreme de Frobenius on garde seulement les 2 premiers axes
Co = Cp(:,1:2) 

scatter(Co(:,1),Co(:,2)); grid on; xlabel('Axe 1'); ylabel('Axe 2');
title('Image euclidienne des objets')
for i=1:n
    text(Co(i,1), Co(i,2),varnames(i));
end

end

function [r] = prod_hs(A,B)
%--------------------------------
% Definition de produit scalaire
%--------------------------------
r = trace(A*B');
end

function [An]= norme(A)
%--------------------------------
% Definition de norme
%--------------------------------
An= prod_hs(A,A);
end

function [XU] = ACP(X)
% - Recherche des valeurs et vecteurs propres
[VEPU, VAPU] = eig(X'*X);    
VAPU         = diag(VAPU);        
%
% - Ordonnancement des valeurs et vecteurs propres
[VAPU,s] = sort(VAPU);
VAPU     = VAPU(flipud(s)); 
VEPU     = VEPU(:,flipud(s)); 
%
% - Nouvelles Coordonn?es (Composantes principales)
XU = X * VEPU; 
end