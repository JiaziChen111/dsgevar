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
X = [y_obs(p:end-1,:), q_obs(p:end-1,:), pi_obs(p:end-1,:), i_obs(p:end-1,:), y_obs(p-1:end-2,:), q_obs(p-1:end-2,:), pi_obs(p-1:end-2,:), i_obs(p-1:end-2,:), ones(t,1)];

kk = t - 130;

C = [1 2 3 4 5 6];
AA = [2 4 5];
for FTv = 1:length(C)
    if (C(FTv)~=AA)
        for j = 1:kk-FTv
            i = kk;
YY = Y(j:end-i+j,:); %rolling window
XX = X(j:end-i+j,:);
phi = (XX'*XX)\XX'*YY; %MQO
X_s = [YY(end,:) YY(end-1,:) ones(1,1)];
yfor = X_s*(phi.^FTv);
Yactual = Y(end-i+j+FTv,:);
L = (Yactual - yfor).^2;
    if FTv == 1 
    L1v(j,:) = L;
        elseif FTv == 3
        L3v(j,:) = L;
    else 
    L6v(j,:) = L;  
    end
        end
    end
end