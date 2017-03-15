function [Co,S,SS,RV,W,Wn,VaP,VeP,p] = statis_inter (X,M,Delta,Sup,reduit,norm,D,etunames)
%% Fonction de calcul de de l'interstructure pour la methode STATIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input variables
% X = Tableaux avec les t ?tudes
% M = Metrique (usuelment matrice identit?)
% Delta = Matrice diagonal avec les poids
% Sup = Matrice avec les tableaux supplementaires
% reduit= option si on veut reduire nos données
% norm = Option si on veut faire l'analyse en prendre en compre la norme
%
% PARAMETRES OPTIONELS
% varnames = variable de type string qui a le nom des variables
% D = M?trique des poids, permettant le calcul des distances entre variables,
%     usuelment 1/n * I (I est la matrice identit?)
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
% Xc = Donn?es centr?es et reduites
%
% Use:
% [Co,S,SS,RV,W,Wn,VaP,VeP,Xc] = statis_inter (X,M,Delta,Sup,reduit,norm,D,etunames)
%
% Author: Rodrigo Andres Rivera Martinez
% Corrections: Larbi Mouchou, Mounir Bendali-Braham, Nafise Gouard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Definition des objets representatifs
% Verification des matrices
if(nargin<6)
    error('SYNTAXE ERROR: statis_interstructure (X,M,Delta,Sup,norm,reduit)');
end

[L,C,n] = size(X);
[L3,C3]=size(Sup);

if L~=L3
    error('ERROR: Nb d''individus doit etre egal dans les matrices X et Sup');
end
if C ~= C3
    error('ERROR: Nb de variables doit etre egal dans les matrices X et Sup');
end

if nargin<8
    for i=1:n
        etunames{i} = sprintf('Objet %d',i);
    end
else 
    if size(etunames,2) < n || size(etunames,2) > n
        error('ERROR: <varnames> doit etre de la meme taille que le nb d''etudes');
    end
end
%-------------------------------------------------------------------------------
% Definition des objets
%-------------------------------------------------------------------------------
% Centrage des tableaux
for i = 1:n
    Xc(:,:,i) = centrer(X(:,:,i));
end


% Reduire des tableaux
if ~reduit
    for i=1:n
        Xcr(:,:,i)=Xc(:,:,i)
    end    
else
    for i=1:n        
        Xcr(:,:,i)=reduire(Xc(:,:,i))
    end
end

for i = 1:n
     W(:,:,i) = Xcr(:,:,i)*M*Xcr(:,:,i)';
     Wn(:,:,i) = W(:,:,i)/sqrt(norme(W(:,:,i)));
end


%-------------------------------------------------------------------------------
% Calcul de la matrice des produits scalaires S
%-------------------------------------------------------------------------------
if nargin > 6
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

% Pourcentage d'inertie
p= (VaP*100)/sum(VaP);

figure;
scatter(Co(:,1),Co(:,2)); grid on; 
xlabel(sprintf('Axe 1 (Inertie: %.2f %%)',p(1)));
ylabel(sprintf('Axe 2 (Inertie: %.2f %%)',p(2)));
title('Image euclidienne des objets')


for i=1:n
    text(Co(i,1), Co(i,2),etunames(i));
end

end

function [r] = prod_hs(A,B,D)
%--------------------------------
% Definition de produit scalaire
%--------------------------------
if nargin < 3
    n= size(A,2);
    D = 1/n * eye(n);
end
r = trace(D*A*D*B);
end

function [Ar]=reduire(A)
%--------------------------------
% Reduire les rableaux
%--------------------------------
[n,p]=size(A);
Ar=A.*repmat(sqrt(n-1./(n.*var(A))),n,1);
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

function [XU, VAPU, VEPU] = ACP(X)
%--------------------------------
% Calcul ACP
%--------------------------------
% Recherche des valeurs et vecteurs propres
[VEPU, VAPU] = eig(X);    
VAPU         = diag(VAPU);        
%
% Ordonnancement des valeurs et vecteurs propres
[VAPU,s] = sort(VAPU, 'descend'); 
VEPU     = VEPU(:,s); 
%
% Nouvelles Coordonn?es (Composantes principales)
XU = VEPU * diag(sqrt(VAPU)); 
end

function [Ac] = centrer(A)
%--------------------------------
% Centrage des donn?es
%--------------------------------
 [n,p]=size(A);
 Ac= A - repmat(mean(A),n,1);
end