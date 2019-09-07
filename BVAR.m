clear all
close all
clc

prompt = 'Defina a ordem de defasagem (1, 2, ou 3): ';
p = input(prompt); %Number of lags (1, 2 or 3) 

%Importa os dados
load dados_dissertacao.mat
y_obs = diff(log(IBCBr));
pi_obs = IPCA(2:end,:);
q_obs = diff(log(CAMBIO));
i_obs = JUROS(2:end,:);

%y_{t}
Y = [y_obs(p+1:end), q_obs(p+1:end), pi_obs(p+1:end), i_obs(p+1:end)];
t = length(Y);
[~,n] = size(Y);

kk = t - 130;
result = bvar(Y,p,.2,.2,1);
C = [1 2 3 4 5 6];
AA = [2 4 5];
for FTm = 1:length(C)
    if (C(FTm)~=AA)
        for j = 1:kk-FTm
            i = kk;
B = Y(j:end-i+j,:);
iEnd=length(B)+1;
%gera previsões usando a função toolbox LeSage
yfor = bvarf(B,p,FTm,iEnd,0.2,0.2,1,[],0);
Yactual = Y(end-i+j+FTm,:);
L = (Yactual - yfor(end,:)).^2;
    if FTm == 1 
    L1m(j,:) = L;
        elseif FTm == 3
        L3m(j,:) = L;
    else 
    L6m(j,:) = L;  
    end
        end
    end
end


