% %% STATIS DB 1
% clear; close all; clc;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Load Data 1
% path_data = 'Data/'; 
% filename=[path_data,'nnotes_FAT.xls'];
% Data=xlsread(filename);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Parametres Interstructure
% X = zeros(6,3);
% Data;
% j=1;
% for i = 1:3:11
%     X(:,:,j) = Data(:,i:i+2);
%     j=j+1;
% end
% 
% M = 1/3*eye(size(X,2));
% Sup = X(:,:,4);
% Delta = 1/4*eye(size(X,3));
% norm=1;
% D =(1/size(X,1))*eye(size(X,1));
% varetude = {'Année 1','Année 2','Année 3','Année 4'};
% varnames = {'Francais', 'Maths', 'Histoire'};
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Co,S,SS,RV,W,VaP,VeP,Xc] = statis_inter(X,M,Delta,Sup,norm,D, varetude);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ Wcomp, alpha_t ] = compromis(W,S,Delta,VaP,VeP,norm);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% indnames = {'Eleve 1','Eleve 2','Eleve 3','Eleve 4','Eleve 5','Eleve 6'};
% [ B, B_val_c, Wd, VAPU, VEPU, corrvars, V_pour ] = statis_intra( Xc, M, W, Wcomp, alpha_t, indnames, varetude, varnames, norm, Delta);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trajectoires( X, W, D, VEPU, VAPU, V_pour, indnames)

%% STATIS DB 2
clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data 2
path_data = 'Data/'; 
filename=[path_data,'croiss_tall.xls'];
Data=xlsread(filename);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametres Interstructure
X = zeros(12,9);
for i=1:30
    X(:,:,i) = Data(i:i+11,:);
end

M = (1/size(X,2))*eye(size(X,2));
Sup = X(:,:,4);
Delta = (1/size(X,3))*eye(size(X,3));
norm=0;
D =(1/size(X,1))*eye(size(X,1));

for i=1:size(X,1) varetude{i} = sprintf('Etude %d',i); end
for i=1:size(X,1) indnames{i} = sprintf('Individu %d',i); end
varnames{1}=sprintf('annee');
varnames{2}=sprintf('weight');
varnames{3}=sprintf('taille');
varnames{4}=sprintf('coccys');
varnames{5}=sprintf('poitrine');
varnames{6}=sprintf('bras');
varnames{7}=sprintf('mollet');
varnames{8}=sprintf('tete');
varnames{9}=sprintf('pelvis');
for t=1:size(X,3) varetude{t} = sprintf('Annee %d', t); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Co,S,SS,RV,W,VaP,VeP,Xc] = statis_inter(X,M,Delta,Sup,norm,D, varetude);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ Wcomp, alpha_t ] = compromis(W,S,Delta,VaP,VeP,norm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ B, B_val_c, Wd, VAPU, VEPU, corrvars, V_pour ] = statis_intra( Xc, M, W, Wcomp, alpha_t, indnames, varetude, varnames, norm, Delta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trajectoires( X, W, D, VEPU, VAPU, V_pour, indnames)
