function [ B, B_val_c, Wd, VAPU, VEPU, corrvars, V_pour ] = statis_intra( X, M, W, Wcomp, alpha_t, indnames, varetudes, varnames, norm )
%% Fonction de calcul de l'intrastructure pour la methode STATIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input variables
% X = Tableaux avec les t études
% Wcomp = Matrice avec le Compromis entre les objets W
% indnames = variable de type string qui a le nom des individus
% varetudes = variable de type string qui a les noms des variables
% p = pourcentage de inertie minimum requis (pour choisir le nb d'axes)
%     par default on garde les 2 premiers axes si cette variable n'est pas
%     definie
%
% Output Variables
% B = Image euclidienne compromis des individus
%
% Use:
% [ B ] = statis_intra( Wcomp, indnames,p )
%
% Authors: Larbi Mouchou, Rodrigo Andres Rivera Martinez, Mounir Bendali-Braham, Nafise Gouard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Definition de l'image euclidienne
[L,C] = size(Wcomp);
D = 1/L*eye(L);
Wd = Wcomp*D;

[XU,VAPU, VEPU] = ACP(Wd);
n = size(VAPU, 1);
V_pour = (VAPU*100)/sum(VAPU); %Valeur propre_pourcentage (Inertie) ?
p_tot = sum(V_pour(1:2)); %Pourcentage total ?
j=2;

B = XU(:,1:j); 

disp('Valeur propres');
disp(VAPU);
fprintf('Pourcentage de l''inertie cumule: %.3f %%\n',p_tot);
fprintf('Nb d''axes: %d\n',j);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot de l'image euclidienne compromis des individus
figure;
scatter(B(:,1),B(:,2)); grid on;
xlabel(sprintf('Axe 1 (Inertie: %.2f %%)',V_pour(1)));
ylabel(sprintf('Axe 2 (Inertie: %.2f %%)',V_pour(2)));
title('Image euclidienne compromis des individus')

disp('V_pour');
disp(V_pour);
if nargin <5
    for i=1:L
        indnames{i} = sprintf('Individu %d',i);
    end
end

for i=1:L
    text(B(i,1), B(i,2),indnames(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Corrélations des variables avec les axes du compromis

[nb_inds, nb_vars, nb_etudes] = size(X);
nb_axes = 2;

corrvars = ones(nb_etudes, nb_vars, nb_axes);
[L,C,n] = size(X);
for i = 1:n
    Xc(:,:,i) = X(:,:,i);
end

for axe = 1:nb_axes
    for var= 1:nb_vars
        for etude= 1:nb_etudes
            corrvars(etude, var, axe) = Xc(:,var,etude)'*D*VEPU(:,axe);
        end;
    end;
end;
% Plot des corrélations des variables

figure;
hold on;
for var= 1:nb_vars
    plot(corrvars(:,var,1), corrvars(:,var,2), '-O');
    for t = 1: nb_etudes
        text(corrvars(t,var,1), corrvars(t,var,2), [varnames(var) num2str(t)]);
    end;
end;

xlabel(sprintf('Axe 1 (Inertie: %.2f %%)',V_pour(1)));
ylabel(sprintf('Axe 2 (Inertie: %.2f %%)',V_pour(2)));
title('Corrélations des variables')

grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Validation de l'image euclidienne compromis

if norm
    B_val=(sqrt(abs(alpha_t(1)/norme(Xc(:,:,1)*M*Xc(:,:,1)'))))*Xc(:,:,1);
    for i=2:nb_etudes
        B_val =[B_val (sqrt(abs(alpha_t(i)/norme(Xc(:,:,i)*M*Xc(:,:,i)'))))*Xc(:,:,i)];
    end
else
    B_val=(sqrt(abs(alpha_t(1))))*Xc(:,:,1);
    for i=2:nb_etudes
        B_val =[B_val (sqrt(abs(alpha_t(i))))*Xc(:,:,i)];
    end
end

%B_val = centrer(B_val,mean(B_val), std(B_val));
% ACP du B_val
[XU_v, Vp, VE] = ACP2(B_val,D);
B_val_c = XU_v(:,1:2);

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


disp('norme Wcomp');
disp(norme(Wcomp));

disp('produit Delta * normes(Wt)')
Pi = 1/4;
disp(Pi*norme(W(:,:,1)));
disp(Pi*norme(W(:,:,2)));
disp(Pi*norme(W(:,:,3)));
disp(Pi*norme(W(:,:,4)));



end
function [XU,VAPU, VEPU] = ACP(X)
%--------------------------------
% Calcul ACP
%--------------------------------
% Recherche des valeurs et vecteurs propres
[VEPU, VAPU] = eig(X);
VAPU         = diag(VAPU); 
%VAPU=sign(VAPU(1))*VAPU;
%
% Ordonnancement des valeurs et vecteurs propres
[VAPU,s] = sort(VAPU, 'descend');
%VAPU     = VAPU(s); 
VEPU     = VEPU(:,s);
%
% Nouvelles Coordonnées (Composantes principales)
XU = VEPU *diag(sqrt(VAPU)); 
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

function [XU,VAPU, VEPU] = ACP2(X,D)
%--------------------------------
% Calcul ACP
%--------------------------------
% Recherche des valeurs et vecteurs propres
[VEPU, VAPU] = eig(X*((1/3)*eye(size(X,2)))*X'*D);
VAPU         = diag(VAPU); 
%VAPU=sign(VAPU(1))*VAPU;
VEPU=sign(VEPU(1,1)).*VEPU;
%
% Ordonnancement des valeurs et vecteurs propres
[VAPU,s] = sort(VAPU, 'descend');
%VAPU     = VAPU(s); 
VEPU     = VEPU(:,s);
%
% Nouvelles Coordonnées (Composantes principales)
XU = VEPU *diag(sqrt(VAPU)); 
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
