function [Co,S,SS,RV,W,Wn,VaP,VeP] = statis_interstructure (X,M,Delta,Sup,norm,D,varnames)
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
% D = M?trique des poids, permettant le calcul des distances entre variables,
%     usuelment 1/n * I, (o? I est la matrice identit?)
%
% Output Variables
% Co = Matrice avec les composantes principales
% S = Matrice des produits scalaire entre les tableaux
% SS = Matrice des produits scalaire entre les tableaux afect? par le poids
%      Delta
% RV = Matrice avec les coefficients RV
% W = Matrice avec les objets des t etudes
% Wn = Matrice avec les objets des t etudes normes
% VaP = Valeurs propres du matrice SS
% VeP = Vecteurs propres du matrice SS
%
% Use:
% [Co,S,SS,RV,W,Wn,VaP,VeP] = statis_interstructure (X,M,Delta,Sup,norm,D,varnames)
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
if C ~= C3
    error('ERROR: Nb de variables doit etre egal dans les matrices X et Sup');
end

if nargin<7
    for i=1:n
        varnames{i} = sprintf('Objet %d',i);
    end
else 
    if size(varnames,2) < n || size(varnames,2) > n
        error('ERROR: <varnames> doit etre de la meme taille que le nb d''?tudes');
    end
end
%-------------------------------------------------------------------------------
% Definition des objets
%-------------------------------------------------------------------------------
% Centrage des tableaux

for i = 1:n
    Xc(:,:,i) = centrer(X(:,:,i),mean(mean(X(:,:,i))), std(std(X(:,:,i))));
end

for i = 1:n
     W(:,:,i) = Xc(:,:,i)*M*Xc(:,:,i)';
     Wn(:,:,i) = W(:,:,i)/sqrt(norme(W(:,:,i)));
end
%-------------------------------------------------------------------------------
% Calcul de la matrice des produits scalaires S
%-------------------------------------------------------------------------------
if nargin == 6
    if ~norm
        for i=1:n
            for j=1:n
                S(i,j)=prod_hs(W(:,:,i),W(:,:,j),D);
            end
        end
    else
        for i=1:n
            for j=1:n
                S(i,j)=(prod_hs(W(:,:,i),W(:,:,j),D))/(norme(W(:,:,i),D)*(norme(W(:,:,j),D)));
            end
        end
    end
else 
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

[Cp,VaP,VeP] = ACP(SS);
% Par le th?oreme de Frobenius on garde seulement les 2 premiers axes
Co = Cp(:,1:2); 

figure;
scatter(Co(:,1),Co(:,2)); grid on; xlabel('Axe 1'); ylabel('Axe 2');
title('Image euclidienne des objets')


for i=1:n
    text(Co(i,1), Co(i,2),varnames(i));
end

end

function [r] = prod_hs(A,B,D)
%--------------------------------
% Definition de produit scalaire
%--------------------------------
if nargin < 3
    n= size(A,2);
    D =1/n * eye(n);
end
r = trace(D*A*D*B);
end

function [An]= norme(A,D)
%--------------------------------
% Definition de norme
%--------------------------------
if nargin < 3
    n= size(A,2);
    D =1/n * eye(n);
end

An= sqrt(prod_hs(A,A,D));
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

function [Ac] = centrer(A,mean_A,std_A)
%--------------------------------
% Centrage des donn?es
%--------------------------------
UN = ones(size(A));
Me = UN * mean_A;
Ecart_type = UN * diag(std_A);
Ac  = (A - Me)./Ecart_type;
end