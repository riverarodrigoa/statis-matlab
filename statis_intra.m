function [ B, Wd, VAPU, VEPU, corrvars, V_pour ] = statis_intra( X, Wn, Wcomp, indnames, varetudes, varnames, p )
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
if nargin > 6
    while p > p_tot
        j=j+1;
        if j>n break; end;
        p_tot = p_tot + V_pour(j);
    end
end

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
if nargin <2
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
    Xc(:,:,i) = centrer_reduire(X(:,:,i),mean(mean(X(:,:,i))), std(std(X(:,:,i))));
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
[VAPU,s] = sort(VAPU);
VAPU     = VAPU(flipud(s)); 
VEPU     = VEPU(:,flipud(s)); 
%
% Nouvelles Coordonnées (Composantes principales)
XU = (X * VEPU) * diag(1./sqrt(VAPU)); 
end


function [Ac] = centrer_reduire(A,mean_A,std_A)
%--------------------------------
% Centrage des donn?es
%--------------------------------
UN = ones(size(A));
Me = UN * mean_A;
Ecart_type = UN * diag(std_A);
Ac  = (A - Me)./Ecart_type;
end