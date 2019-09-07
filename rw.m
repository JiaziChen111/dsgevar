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

C = [1 2 3 4 5 6];
AA = [2 4 5];
for h = 1:length(C)
    if (C(h)~=AA)
        for j = 1:kk-h
            i = kk;
yfor = Y(end-i+j,:);
Yactual = Y(end-i+j+h,:);
L = (Yactual - yfor).^2;
    if h == 1 
    L1v(j,:) = L;
        elseif h == 3
        L3v(j,:) = L;
    else 
    L6v(j,:) = L;  
    end
        end
    end
end