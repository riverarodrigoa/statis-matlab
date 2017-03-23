function [ Wcomp, alpha_t ] = compromis(W,S,Delta,VaP,VeP,norm)
%% Fonction de calcul de compromis pour la methode STATIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input variables
% W = Matrice avec les objets des t etudes
% S = Matrice des produits scalaire entre les Objets (W) affectes par le poids
%     Delta
% Delta = Matrice diagonal avec les poids (pi_i)
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
% Author: Rodrigo Andres Rivera Martinez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('***************************************');
disp('********** STATIS  COMPROMIS **********');
disp('***************************************');
% Verification des matrices et definition des parametres par default
if ~exist('W','var') || isempty(W)
    error('[STATIS] You must provide a matrix W');
end
if ~exist('S','var') || isempty(S)
    error('[STATIS] You must provide a matrix S');
end
if ~exist('Delta','var') || isempty(Delta)
    error('[STATIS] You must provide a matrix Delta');
end
if ~exist('VaP','var') || isempty(VaP)
    error('[STATIS] You must provide the eigen values');
end
if ~exist('VeP','var') || isempty(VeP)
    error('[STATIS] You must provide the eigen vectors');
end
if ~exist('norm','var') || isempty(norm)
    norm=1; 
    disp('[DEFAULT] Norme = 1');
else
    if norm
        disp('[USER] Norme');
    else
        disp('[USER] Non norme');
    end
end
%% Expresion definitive du compromis
[~,~,T] = size(W);
pi_t=diag(Delta);
SS = diag(S);
gamma = VeP(:,1);
Wcomp = 0;
alpha_t = zeros(1,T);
if ~norm
    alpha_c = (sqrt(pi_t')*sqrt(SS));
    for i = 1:T
        alpha_t(i) = (1./sqrt(VaP(1))*alpha_c*sqrt(pi_t(i))*gamma(i));
      Wcomp = Wcomp + (alpha_t(i)*W(:,:,i));
    end
else
    for i = 1:T
        alpha_t(i) = (1/sqrt(VaP(1))*sqrt(pi_t(i))*gamma(i));
        Wcomp = Wcomp + (alpha_t(i)*W(:,:,i));
    end
end
fprintf('Compromis [W] [%d x %d]\n',size(Wcomp));
disp(Wcomp);
disp('***************************************');
end
