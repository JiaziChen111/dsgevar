clear all
close all
clc
%Elaborando um MCS para determinar qual conjunto de modelos, M*, consiste 
%nos melhores modelos da coleção de modelos M^0. 
addpath('C:\Users\Fernanda\Google Drive\Programa\DSGE-VAR');
load forecasting.mat
addpath('C:\Users\Fernanda\Google Drive\Programa\BVAR-Minn');
load bvar.mat
load rw.mat
n=4;
for j=4

L_1 = L1(:,j,:);
L_1 = L_1(:,:);
L1_m = L1m(:,j);
L1c = [L_1, L1_m]; %losses, TxK
L1_c(:,:,j) = L1c;

L_3 = L3(:,j,:);
L_3 = L_3(:,:);
L3_m = L3m(:,j);
L3c = [L_3, L3_m]; %losses, TxK
L3_c(:,:,j) = L3c;

L_6 = L6(:,j,:);
L_6 = L_6(:,:);
L6_m = L6m(:,j);
L6c = [L_6, L6_m]; %losses, TxK
L6_c(:,:,j) = L6c;
end
a = .05; %o MCS contém o melhor modelo com probabilidade 1-a, onde a=0.05 
%é o nível de significância

B = 999; %Hansen
w = 12; %Hansen 
for j=4
[includedR1,pvalsR1,excludedR1]=mcs(L1_c(:,:,j),a,B,w);
[includedR3,pvalsR3,excludedR3]=mcs(L3_c(:,:,j),a,B,w);
[includedR6,pvalsR6,excludedR6]=mcs(L6_c(:,:,j),a,B,w);
end