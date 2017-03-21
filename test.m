%% STATIS DB 1
clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data 1
path_data = '~/Documents/UNIVERSITE_PARIS_SACLAY/M2_TRIED/Projet_long/Data/'; 
filename=[path_data,'nnotes_FAT.xls'];
Data=xlsread(filename);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametres
X = zeros(6,3);
j=1;
for i = 1:3:11 
    X(:,:,j) = Data(:,i:i+2);
    j=j+1;
end

norm=1;
Delta = eye(size(X,3));
D =1/size(X,1) * eye(size(X,1));

etunames{size(X,3)}=[];
for i=1:size(X,3)
    etunames{i} = sprintf('Ann?e %d',i);
end
varnames = {'Francais', 'Math', 'Histoire'};
indnames = {'Eleve 1','Eleve 2','Eleve 3','Eleve 4','Eleve 5','Eleve 6'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Co,SS,RV,W,VaP,VeP,Xcr]=statis_inter(X,[],[],[],norm,etunames);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ Wcomp, alpha_t ] = compromis(W,SS,Delta,VaP,VeP,norm);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ B, Wd, VAPU, VEPU, corrvars, V_pour ] = statis_intra( X, W, Wcomp, alpha_t, indnames, etunames, varnames, 80 );
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trajectoires( X, W, D, VEPU, VAPU, V_pour, indnames )
%%
%% STATIS DB 2
clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data 2
path_data = '~/Documents/UNIVERSITE_PARIS_SACLAY/M2_TRIED/Projet_long/Data/'; 
filename=[path_data,'croiss_tall.xls'];
Data=xlsread(filename);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametres Interstructure

for(i=1:size(Data,1)/12)
     X(:,:,i)=Data(12*(i-1)+1:12*(i),:);
end

X=X(:,2:end,:);
M = eye(size(X,2));
Sup = X(:,:,4);
Delta = eye(size(X,3));
norm=1;
D =1/size(X,1) * eye(size(X,1));

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
% for t=1:size(X,3) varetude{t} = sprintf('Annee %d', t); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Co,SS,RV,W,VaP,VeP,Xcr] = statis_inter (X);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ Wcomp ] = compromis(W,SS,Delta,VaP,VeP,norm);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ B, Wd, VAPU, VEPU, corrvars, V_pour ] = statis_intra( X, Wn, Wcomp, indnames, varetude, varnames, 99 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%7
trajectoires( X, Wn, D, VEPU, VAPU, V_pour, indnames )