function [Co,S,SS,RV,W,VaP,VeP,Xc] = statis_inter (X,M,Delta,Sup,norm,D,varnames)
%% Fonction de calcul de de l'interstructure pour la methode STATIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input variables
% X = Tableaux avec les t études
% M = Metrique (usuelment matrice identit?)
% Delta = Matrice diagonal avec les poids
% Sup = Matrice avec les tableaux supplementaires
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
% [Co,S,SS,RV,W,Wn,VaP,VeP] = statis_inter (X,M,Delta,Sup,norm,D,varnames)
%
% Authors: Larbi Mouchou, Rodrigo Andres Rivera Martinez, Mounir Bendali-Braham, Nafise Gouard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% D?finition des objets repr?sentatifs
% Verification des matrices
if(nargin<5)
    error('SYNTAXE ERROR: statis_interstructure (X,M,Delta,Sup,norm)');
end

[L,C,n] = size(X);
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
    Xc(:,:,i) = centrer(X(:,:,i),mean(X(:,:,i)), std(X(:,:,i)));
end

for i = 1:n
    if ~norm
       W(:,:,i) = Xc(:,:,i)*M*Xc(:,:,i)';
    else
        W(:,:,i) = Xc(:,:,i)*M*Xc(:,:,i)';
        W(:,:,i) = W(:,:,i)./norme(W(:,:,i));
    end     
end
%-------------------------------------------------------------------------------
% Calcul de la matrice des produits scalaires S
%-------------------------------------------------------------------------------
if nargin > 5
    for i=1:n
        for j=1:n
            S(i,j)=prod_hs(W(:,:,i),W(:,:,j),D);
        end
    end
else 
    for i=1:n
        for j=1:n
            S(i,j)=prod_hs(W(:,:,i),W(:,:,j));
        end
    end
end
%-------------------------------------------------------------------------------
% Calcul de coefficient RV
%-------------------------------------------------------------------------------
for i=1:n
    for j=1:n
        RV(i,j) = S(i,j)/(sqrt(S(i,i))*sqrt(S(j,j)));
    end
end
%-------------------------------------------------------------------------------
% Image euclidienne des objets
%-------------------------------------------------------------------------------
SS = Delta'*S*Delta;

[Cp,VaP,VeP] = ACP(SS);
% Par le th?oreme de Frobenius on garde seulement les 2 premiers axes
Co = Cp(:,1:2);

% Pourcentage d'inertie
p= (VaP*100)/sum(VaP);

figure;
scatter(Co(:,1),Co(:,2)); grid on; 
axis([0 (max(Co(:,1))*1.3) (min(Co(:,2))*2) (max(Co(:,2))*2)])
xlabel(sprintf('Axe 1 (Inertie: %.2f %%)',p(1)));
ylabel(sprintf('Axe 2 (Inertie: %.2f %%)',p(2)));
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
    D = 1/n * eye(n);
end
r = trace(A*D*B*D);
end

function [An]= norme(A,D)
%--------------------------------
% Definition de norme
%--------------------------------
if nargin < 2
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
[VEPU, VAPU] = eig(X); 
VAPU         = diag(VAPU);        
%
% Ordonnancement des valeurs et vecteurs propres
[VAPU,s] = sort(VAPU, 'descend');
%VAPU     = VAPU(s); 
VEPU     = VEPU(:,s);
VEPU=sign(VEPU(1,1)).*VEPU;
%
% Nouvelles Coordonn?es (Composantes principales)
XU = VEPU * diag(sqrt(VAPU));
%XU=sign(XU(1,1)).*XU;
end

function [Ac] = centrer(A,mean_A,std_A)
%--------------------------------
% Centrage des donn?es
%--------------------------------
UN = ones(size(A));
Me = UN * diag(mean_A);
Ecart_type = UN * diag(std_A);
ec  = 1./Ecart_type;
Ac  = (A - Me).*ec;
end
