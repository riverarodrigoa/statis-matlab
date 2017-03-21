%% STATIS DB 1
clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
Delta = 1/4*eye(size(X,3));
norm=1;
D =(1/size(X,1))*eye(size(X,1));
varetude = {'Année 1','Année 2','Année 3','Année 4'};
varnames = {'Francais', 'Maths', 'Histoire'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Co,S,SS,RV,W,VaP,VeP,Xc] = statis_inter(X,M,Delta,Sup,norm,D, varetude);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ Wcomp, alpha_t ] = compromis(W,S,Delta,VaP,VeP,norm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indnames = {'Eleve 1','Eleve 2','Eleve 3','Eleve 4','Eleve 5','Eleve 6'};
[ B, B_val_c, Wd, VAPU, VEPU, corrvars, V_pour ] = statis_intra( Xc, M, W, Wcomp, alpha_t, indnames, varetude, varnames, norm, Delta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trajectoires( X, W, D, VEPU, VAPU, V_pour, indnames)
