function [Co,S,SS,RV,W,Wn,VaP,VeP] = statis_inter (X,M,D,Delta,norm,r,etunames)
%% Fonction de calcul de de l'interstructure pour la methode STATIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input variables
% X = Tableaux avec les t ?tudes
%
% PARAMETRES OPTIONELS
% M = Metrique (usuelment matrice identit?)
% Delta = Matrice diagonal avec les poids
% Sup = Matrice avec les tableaux supplementaires
% norm = Option si on veut faire l'analyse en prendre en compre la norme
% varnames = variable de type string qui a le nom des variables
% D = M?trique des poids, permettant le calcul des distances entre variables,
%     usuelment 1/n * I (I est la matrice identit?)
% r = Centrage et/ou reduction (1: Centrage et reduction, 0:Centrage)
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
% [Co,S,SS,RV,W,Wn,VaP,VeP] = statis_inter (X,M,D,Delta,norm,r,etunames)
%
% Author: Rodrigo Andres Rivera Martinez
% Corrections: Larbi Mouchou, Mounir Bendali-Braham, Nafise Gouard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('***************************************');
disp('******** STATIS INTERSTRUCTURE ********');
disp('***************************************');
%% Definition des objets repr?sentatifs
% Verification des matrices et definition des parametres par default
if ~exist('X','var') || isempty(X)
    error('[STATIS] You must provide a matrix X');
else
    disp('****************************');
    disp('***Parametres par default***');
    disp('****************************');
    disp('');
    [L,C,n] = size(X);
    if ~exist('M','var') || isempty(M)
        M = eye(C);
        fprintf('[DEFAULT] M = I [%d x %d]\n',size(M));
    end
    if ~exist('D','var') || isempty(D)
        D = (1/L).* eye(L);
        fprintf('[DEFAULT] D = 1/%d *I [%d x %d]\n',L,size(D));
    end
    if ~exist('Delta','var') || isempty(Delta)
        Delta = (1/n).*eye(n);
        fprintf('[DEFAULT] Delta = 1/%d *I [%d x %d]\n',n,size(Delta));
    end
    if ~exist('norm','var') || isempty(norm) % norm? par default
        norm = 1;
        fprintf('[DEFAULT] Norme = %d\n',norm);
    end
    if ~exist('r','var') || isempty(r) % reduit par default
        r = 1;
        fprintf('[DEFAULT] Centrage et reduction = %d\n',r);
    end
    if ~exist('etunames','var') || isempty(etunames)
        for i=1:n
            etunames{i} = sprintf('Objet %d',i);
        end
        fprintf('[DEFAULT] etunames [%d x %d]\n',size(etunames));
    else
        if size(etunames,2) < n || size(etunames,2) > n
        error('ERROR: <etunames> doit etre de la meme taille que le nb d''?tudes');
        end
    end
end
disp('****************************');
%-------------------------------------------------------------------------------
% Definition des objets
%-------------------------------------------------------------------------------
% Centrage/reduction des tableaux
for i = 1:n
    Xc(:,:,i) = centrer(X(:,:,i),r);
end

% Calcul des objets
W = zeros(L,L,n);
if norm
    for i = 1:n
        W(:,:,i) = Xc(:,:,i)*M*Xc(:,:,i)';
        Wn(:,:,i) = W(:,:,i)./sqrt(norme(W(:,:,i),D));
    end
else
    for i = 1:n
        W(:,:,i) = Xc(:,:,i)*M*Xc(:,:,i)';
        Wn(:,:,i)= Xc(:,:,i)*M*Xc(:,:,i)';
    end
   
end
%-------------------------------------------------------------------------------
% Calcul de la matrice des produits scalaires S
%-------------------------------------------------------------------------------
S = zeros(n);
if ~norm
    for i=1:n
        for j=1:n
            S(i,j)=prod_hs(W(:,:,i),W(:,:,j),D);
        end
    end
else
    for i=1:n
        for j=1:n
            S(i,j)=(prod_hs(W(:,:,i),W(:,:,j),D))./(norme(W(:,:,i),D)*(norme(W(:,:,j),D)));
        end
    end
end
%-------------------------------------------------------------------------------
% Calcul de coefficient RV
%-------------------------------------------------------------------------------
RV = zeros(n);
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
disp(Co)
% Pourcentage d'inertie
p = (VaP*100)/sum(VaP);

figure;
scatter(Co(:,1),Co(:,2)); grid on;
if norm; xlim([0,1]); ylim([-1,1]); else xlim([0,Inf]);  end
xlabel(sprintf('Axe 1 (Inertie: %.2f %%)',p(1)));
ylabel(sprintf('Axe 2 (Inertie: %.2f %%)',p(2)));
title('Image euclidienne des objets','FontSize',16);

for i=1:n
    text(Co(i,1), Co(i,2),etunames(i));
end

end

function [r] = prod_hs(A,B,D)
%--------------------------------
% Definition de produit scalaire
%--------------------------------
r = trace(D*A*D*B);
end

function [An]= norme(A,D)
%--------------------------------
% Definition de norme
%--------------------------------
An= sqrt(prod_hs(A,A,D));
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
[VAPU,s] = sort(VAPU, 'descend'); 
VEPU     = VEPU(:,s); 
%
% Nouvelles Coordonn?es (Composantes principales)
XU = VEPU * diag(sqrt(VAPU)); 
end

function [Acr] = centrer(A,r)
%--------------------------------
% Centrage des donn?es
%--------------------------------

 [n,p]=size(A);
 
 if r
     %--------------------------------
     % Reduire les rableaux
     %--------------------------------
     Ac= A - repmat(mean(A),n,1);
     Acr=Ac.*repmat(sqrt(n-1./(n.*var(Ac))),n,1);
 else
     Acr= A - repmat(mean(A),n,1); 
 end
end
