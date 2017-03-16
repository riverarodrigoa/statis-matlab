function [Co,SS,RV,W,VaP,VeP,Xcr] = statis_inter (X,M,D,Delta,norm,etunames)
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
% SS = Matrice des produits scalaire entre les tableaux afect? par le poids Delta
% RV = Matrice avec les coefficients RV
% W = Matrice avec les objets des t etudes
% Wn = Matrice avec les objets des t etudes normes
% VaP = Valeurs propres du matrice SS
% VeP = Vecteurs propres du matrice SS
% Xcr = Tableaux des etudes centrees et reduites
%
% Use:
% [Co,SS,RV,W,VaP,VeP,Xcr] = statis_inter (X,M,D,Delta,norm,etunames)
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
    disp('----------------------------');
    disp('-- Parametres par default --');
    disp('----------------------------');
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
        disp('[DEFAULT] Norme');
    else
        if norm
            disp('[USER] Norme');
        else
            disp('[USER] Non norme');
        end
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
disp('***************************************');
%-------------------------------------------------------------------------------
% Definition des objets
%-------------------------------------------------------------------------------
% Centrage/reduction des tableaux
Xcr = zeros(size(X));
for i = 1:n
    Xcr(:,:,i) = centrer(X(:,:,i));
end

% Calcul des objets
W = zeros(L,L,n);
if norm
    for i = 1:n
        W(:,:,i) = Xcr(:,:,i)*M*Xcr(:,:,i)';
        W(:,:,i) = W(:,:,i)./norme(W(:,:,i));
    end
else
    for i = 1:n
        W(:,:,i) = Xcr(:,:,i)*M*Xcr(:,:,i)';
    end
end
%-------------------------------------------------------------------------------
% Calcul de la matrice des produits scalaires S
%-------------------------------------------------------------------------------
S = zeros(n);
for i=1:n
    for j=1:n
        S(i,j)=prod_hs(W(:,:,i),W(:,:,j),D);
    end
end
disp('S');
disp(S)
%-------------------------------------------------------------------------------
% Calcul de coefficient RV
%-------------------------------------------------------------------------------
RV = zeros(n);
if norm
    RV = S;
else
    for i=1:n
        for j=1:n
            RV(i,j) = S(i,j)/(sqrt(S(i,i))*(sqrt(S(j,j))));
        end
    end
end
disp('RV');
disp(RV)
%-------------------------------------------------------------------------------
% Image euclidienne des objets
%-------------------------------------------------------------------------------
SS = S*Delta;

[Cp,VaP,VeP] = ACP(SS);

% Par le th?oreme de Frobenius on garde seulement les 2 premiers axes
Co = Cp(:,1:2);
% Pourcentage d'inertie
p = (VaP*100)/sum(VaP);

% Plot de l'image eucidenne
disp('Plot de l''image euclidienne');

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

function [An]= norme(A)
%--------------------------------
% Definition de norme
%--------------------------------
VaP = eig(A);
An = sqrt(sum(VaP.^2));
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

function [Acr] = centrer(A)
%----------------------------------
% Centrage et reduction des donnees
%----------------------------------
[n,p]=size(A);
Ac= A - repmat(mean(A),n,1);
Acr=Ac.*repmat(sqrt(n-1./(n.*var(Ac))),n,1);
end