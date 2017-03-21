function [ B, Wd, VAPU, VEPU, corrvars, V_pour ] = statis_intra( X, Wn, Wcomp, alpha_t, indnames, varetudes, varnames, p )
%% Fonction de calcul de l'intrastructure pour la methode STATIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input variables
% Wcomp = Matrice avec le Compromis entre les objets W
% X = Matrice des variables
% indnames = variable de type string qui a le nom des individus
% varnames = variable de type string qui a le nom des variables
% p = pourcentage de inertie minimum requis (pour choisir le nb d'axes)
%     par default on garde les 2 premiers axes si cette variable n'est pas
%     definie
%
% Output Variables
% B = Image euclidienne compromis des individus
% Corr = Correlation des variables avec les axes du compromis
%        Matrice N x M x T qui a les correlations o?:
%        - N nb des vecteurs propes
%        - M nb nombre des variables
%        - T nb d'?tudes
%
% Use:
% [ B ] = statis_intrastatis_intra( Wcomp,Xc,p,indnames, varnames)
%
% Author: Rodrigo Andres Rivera Martinez
% Corrections: Larbi Mouchou, Mounir Bendali-Braham, Nafise Gouard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Centrage et reduction des donn?es 
[L,C,T] = size(X);
for i=1:T
    Xc(:,:,i) = centrer(X(:,:,i),mean(mean(X(:,:,i))), std(std(X(:,:,i))),0);
end
%% Definition de l'image euclidienne
[L,C] = size(Wcomp);
D = 1/L*eye(L);
Wd = Wcomp*D;

[XU,VAPU, VEPU] = ACP(Wd);
n = size(VAPU,1);

for i = 1:n
    CP_Wd(:,i)= (1./sqrt(VAPU(i)))*Wd*VEPU(:,i);
end

V_pour = (VAPU*100)/sum(VAPU);
p_tot = sum(V_pour(1:2));

j=2;
if nargin < 3
    while p > p_tot
        j=j+1;
        if j>=n; break; end;
        p_tot = p_tot + V_pour(j);
    end
end

B = CP_Wd(:,1:j);

disp('Valeur propres');
disp(VAPU);
fprintf('Pourcentage de l''inertie cumule: %.3f %%\n',p_tot);
fprintf('Nb d''axes: %d\n',j);
%% Image eucclidienne
% B = VEPU * diag(sqrt(VAPU));
% B = B(:,1:j);
%% Validation de l'image euclidienne compromis
alpha = sqrt(alpha_t);
B_val = alpha(1)*Xc(:,:,1);
for i=2:T
    B_val =[B_val [alpha(i)*Xc(:,:,i)]];
end
size(B_val)
% ACP du B_val
[XU_v, ~, ~] = ACP4(B_val);

B_val_c = XU_v(:,1:j);

Decision = isequal(B,B_val_c);
if Decision
    disp('[Validation test] Image euclidienne compromis correcte');
else
    disp('[Validation test] Image euclidienne compromis incorrecte');
    disp('Image euclidienne');
    disp(B);
    disp('Validation');
    disp(B_val_c);
end
%% Correlations des variables avec les axes du compromis
[L1,C1,n] = size(Xc);

for l=1:j % axes
    for i = 1:C1 % variables
        for k = 1:n % etudes
            corrvars(k,i,l)=Xc(:,i,k)'*D*VEPU(:,l);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot de l'image euclidienne compromis des individus
if nargin <=4
    for i=1:L
        indnames{i} = sprintf('Individu %d',i);
    end
end

for m=1:j
    for k=2:j
        if m==k
        else
        figure;
        scatter(B(:,m),B(:,k)); grid on;
        xlabel(sprintf('Axe %d (Inertie: %.2f %%)',m,V_pour(m)));
        ylabel(sprintf('Axe %d (Inertie: %.2f %%)',k,V_pour(k)));
        title(sprintf('Image euclidienne compromis des individus (Plan %d - %d)',m,k));
            for i=1:L
                text(B(i,1), B(i,2),indnames(i));
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot des correlations des variables
[nb_inds, nb_vars, nb_etudes] = size(X);
figure; hold on; grid on;
for var= 1:nb_vars
  plot(corrvars(:,var,1), corrvars(:,var,2), '-O');
 for t = 1: nb_etudes
     text(corrvars(t,var,1), corrvars(t,var,2), [varnames(var) num2str(t)]);
 end
end

xlabel(sprintf('Axe 1 (Inertie: %.2f %%)',V_pour(1)));
ylabel(sprintf('Axe 2 (Inertie: %.2f %%)',V_pour(2)));
title('Corr?lations des variables')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function [XU, VAPU, VEPU] = ACP(X)
%--------------------------------
% Calcul ACP modifie
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
XU = (X * VEPU) * diag(1./sqrt(VAPU)); 
end

function [XU, VAPU, VEPU] = ACP2(X)
%--------------------------------
% Calcul ACP classique
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

function [XU, VAPU, VEPU] = ACP3(X)
%--------------------------------
% Calcul ACP classique mod
%--------------------------------
% Recherche des valeurs et vecteurs propres
[VEPU, VAPU] = eig(X'*X);    
VAPU         = diag(VAPU);        
%
% Ordonnancement des valeurs et vecteurs propres
[VAPU,s] = sort(VAPU, 'descend');
%VAPU     = VAPU(s); 
VEPU     = VEPU(:,s);
%
% Nouvelles Coordonn?es (Composantes principales)
XU = (X * VEPU) * diag(1./sqrt(VAPU)); 
end

function [XU, VAPU, VEPU] = ACP4(X)
%--------------------------------
% Calcul ACP classique mod2
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


function [Ac] = centrer(A,mean_A,std_A,r)
%--------------------------------
% Centrage des donn?es
%--------------------------------
UN = ones(size(A));
Me = UN * mean_A;
    if r
        Ecart_type = UN * diag(std_A);
        Ac  = (A - Me)./Ecart_type;
    else
        Ac  = (A - Me);
    end
end