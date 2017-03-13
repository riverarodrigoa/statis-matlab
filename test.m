%% STATIS DB 1
clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data 1
path_data = 'Data/'; 
filename=[path_data,'nnotes_FAT.xls'];
Data=xlsread(filename);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametres Interstructure
X = zeros(6,3);
Data;
j=1;
for i = 1:3:11
    X(:,:,j) = Data(:,i:i+2);
    j=j+1;
end

M = eye(size(X,2));
Sup = X(:,:,4);
Delta = [1/4 0 0 0; 0 1/4 0 0; 0 0 1/4 0; 0 0 0 1/4];
norm=1;
D =1/6 * eye(6);
varetude = {'Année 1','Année 2','Année 3','Année 4'};
varnames = {'Francais', 'Maths', 'Histoire'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Co,S,SS,RV,W,Wn,VaP,VeP] = statis_inter(X,M,Delta,Sup,norm,D, varetude);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ Wcomp ] = compromis(Wn,SS,Delta,VaP,VeP,norm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indnames = {'Eleve 1','Eleve 2','Eleve 3','Eleve 4','Eleve 5','Eleve 6'};
[ B, Wd, VAPU, VEPU, corrvars ] = statis_intra( X, Wn, Wcomp, indnames, varetude, varnames, 99 )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% STATIS DB 2
% clear; close all; clc;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Load Data 2
% path_data = 'Data/'; 
% filename=[path_data,'croiss_tall.xls'];
% Data=xlsread(filename);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Parametres Interstructure
% X = zeros(12,9);
% for i=1:30
%     X(:,:,i) = Data(i:i+11,:);
% end
% 
% M = eye(size(X,2));
% Sup = X(:,:,4);
% Delta = eye(size(X,3));
% norm=1;
% D =1/size(X,2) * eye(size(X,2));
% %varnames = {'Ann?e 1','Ann?e 2','Ann?e 3','Ann?e 4'};
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Co,S,SS,RV,W,Wn,VaP,VeP] = statis_inter(X,M,Delta,Sup,norm);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ Wcomp ] = compromis(Wn,SS,Delta,VaP,VeP,norm);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:size(X,1) indnames{i} = sprintf('Individu %d',i); end
% varnames{1}=sprintf('annee');
% varnames{2}=sprintf('weight');
% varnames{3}=sprintf('taille');
% varnames{4}=sprintf('coccys');
% varnames{5}=sprintf('poitrine');
% varnames{6}=sprintf('bras');
% varnames{7}=sprintf('mollet');
% varnames{8}=sprintf('tete');
% varnames{9}=sprintf('pelvis');
% [ B, Wd, VAPU, VEPU, corrvars ] = statis_intra( X, Wn, Wcomp, indnames, varnames, 99);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%