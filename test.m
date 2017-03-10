%% STATIS DB 1
clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data 1
path_data = '~/Documents/UNIVERSITE_PARIS_SACLAY/M2_TRIED/Projet long/Data/'; 
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
Delta = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
norm=1;
D =1/3 * eye(3);
varnames = {'Ann?e 1','Ann?e 2','Ann?e 3','Ann?e 4'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Co,S,SS,RV,W,Wn,VaP,VeP] = statis_inter(X,M,Delta,Sup,norm,D,varnames);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ Wcomp ] = compromis(Wn,SS,Delta,VaP,VeP,norm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varnames = {'Eleve 1','Eleve 2','Eleve 3','Eleve 4','Eleve 5','Eleve 6'};
[ B ] = statis_intra( Wcomp,varnames,99);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% STATIS DB 2
clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data 2
path_data = '~/Documents/UNIVERSITE_PARIS_SACLAY/M2_TRIED/Projet long/Data/'; 
filename=[path_data,'croiss_tall.xls'];
Data=xlsread(filename);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametres Interstructure
X = zeros(12,9);
for i=1:30
    X(:,:,i) = Data(i:i+11,:);
end

M = eye(size(X,2));
Sup = X(:,:,4);
Delta = eye(size(X,3));
norm=1;
D =1/size(X,2) * eye(size(X,2));
%varnames = {'Ann?e 1','Ann?e 2','Ann?e 3','Ann?e 4'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Co,S,SS,RV,W,Wn,VaP,VeP] = statis_inter(X,M,Delta,Sup,norm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ Wcomp ] = compromis(Wn,SS,Delta,VaP,VeP,norm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size(X,1) varnames{i} = sprintf('Individu %d',i); end
[ B ] = statis_intra( Wcomp,varnames,99);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%