clear all
close all
clc

tic;
prompt = 'Defina a ordem de defasagem p: ';
p = input(prompt);

prompt = 'Defina o grid de \lambda. Formato: [\lambda_1, \lambda_2, ... ,\lambda_q]: ';
lambda_s = input(prompt); %Define o grid com diferentes \lambdas

%Importa os dados
load dados_dissertacao.mat
y_obs = diff(log(IBCBr));
pi_obs = IPCA(2:end,:);
q_obs = diff(log(CAMBIO));
i_obs = JUROS(2:end,:);

C = [1 2 3 4 5 6]; %Horizontes de previsão
AA = [2 4 5]; %Horizontes de previsão que não serão avaliados
n = 4; %Número de variáveis observadas
k = n*p+1; %Número de parâmetros + constante  
pp = 16; %Número de parâmetros a serem estimados no DSGE (preencher manualmente)
t = size(y_obs,1)-p; %Número de observações da amostra
kk = t - 130; %134 %Invervalo inicial da janela expansiva
[~,lsize] = size(lambda_s);

%Define a priori
prior.lower = [0;0;0;0;0;0;0;0;0;0;0;0; -inf; -inf; -inf;0];
prior.upper = [inf; inf; inf; inf; inf; inf; 1; inf; inf;1;1;1;inf;inf;inf;inf];
prior.pshape = [4;4;4;4;2;2;1;2;2;1;1;1;3;3;3;2]; %4: inv_gamma, 2: gamma, 1: beta, 3: normal
prior.mean = [0.01; 0.01; 0.01; 0.01;1;0.25; 0.5; 1.75; 0.5; 0.5; 0.5; 0.5; 0.1; .5;-0.1; 1];
prior.var = [4; 4; 4; 4; 0.1; 0.1; 0.15; 0.15; 0.1; 0.15; 0.15; 0.15; 0.15; 0.15; 0.15; 0.5];

for index_lambda = 1:lsize
%Razão entre as variáveis simuladas do DSGE e amostradas
lambda = lambda_s(1,index_lambda);
%Matriz de dados amostrados
YY = [y_obs(p+1:end), q_obs(p+1:end), pi_obs(p+1:end), i_obs(p+1:end)];
%y_{t-h}
if p == 1;
XX = [y_obs(p:end-1,:), q_obs(p:end-1,:), pi_obs(p:end-1,:), i_obs(p:end-1,:), ones(t,1)];
    elseif p == 2;
XX = [y_obs(p:end-1,:), q_obs(p:end-1,:), pi_obs(p:end-1,:), i_obs(p:end-1,:), y_obs(p-1:end-2,:), q_obs(p-1:end-2,:), pi_obs(p-1:end-2,:), i_obs(p-1:end-2,:), ones(t,1)];    
else
XX = [y_obs(p:end-1,:), q_obs(p:end-1,:), pi_obs(p:end-1,:), i_obs(p:end-1,:), y_obs(p-1:end-2,:), q_obs(p-1:end-2,:), pi_obs(p-1:end-2,:), i_obs(p-1:end-2,:), y_obs(p-2:end-3,:), q_obs(p-2:end-3,:), pi_obs(p-2:end-3,:), i_obs(p-2:end-3,:), ones(t,1)];        
end;


    for h= 1:length(C) %Define os h diferentes horizontes de previsão (1, 3 e 6)
        if (C(h)~=AA)
            
            for j = 1:kk-h %Janela expansiva
                i = kk;
            Y = YY(j:end-i+j,:);
            X = XX(j:end-i+j,:);
            
                if lambda == 0; %\lambda = 0, usa-se somente os dados e a priori será degenerada
                npd = 0;
                %Posterior
                %Parâmetros da distribuição Wishart-invertida
                s0_post = Y'*Y - (Y'*X)*inv(X'*X)*(X'*Y); %matriz de covariância
                v0_post = (lambda+1)*t-k; %graus de liberdade      
                [~,g] = cholcov(s0_post);
                    if g ~=0
                    s0_post = nearestSPD(s0_post);
                    npd = npd+1;
                    end;
                Sigma_u_post = iwishrnd(s0_post,v0_post); %sigma_u | theta

                %Parâmetros da distribuição Normal
                Phi = inv(X'*X)*(X'*Y);
                mu_post = vec(inv(X'*X)*(X'*Y)); %média
                Sigma_post = kron(Sigma_u_post, inv(X'*X)); %matriz de covariância
                [~,g] = cholcov(Sigma_post);
                    if g ~=0
                    Sigma_post = nearestSPD(Sigma_post);
                    npd = npd+1;
                    end;
                phi_post = mvnrnd(mu_post, Sigma_post); %phi| sigma_u, theta
                phi_post = reshape(phi_post,k,n);
%--------------------------------Previsões--------------------------------%                
                    %Para h=1, tem-se a preditiva com fórmula fechada.
                    if h==1
                    [Y_s] = forecast_one_step(Y, Phi, p); %(Phi_tilde)
                    Y_for(:,:,j)=Y_s;
                    else
                            if p ==1
                            X_s = [Y(end,:) ones(1,1)];   
                            elseif p ==2
                            X_s = [Y(end,:) Y(end-1,:) ones(1,1)];    
                            else
                            X_s = [Y(end,:) Y(end-1,:) Y(end-2,:) ones(1,1)];    
                            end
                            
                            for jj=1:h
                            Yh = X_s*phi_post+mvnrnd(zeros(1,n),Sigma_u_post);
                                if p == 1
                                X_s = [Yh, ones(1,1)];
                                elseif p ==2
                                V = X_s(:,1:4);
                                X_s = [Yh, V, ones(1,1)];
                                else
                                V = X_s(:,1:8);
                                X_s = [Yh, V, ones(1,1)];    
                                end
                            
                            Y_forh(:,jj) = Yh;
                            end
                                if h==3
                                Y_h = Y_forh(:,end);
                                    elseif h==6
                                    Y_h = Y_forh(:,end);
                                end
                            if h==3
                            Y_for_3(:,:,j) = Y_h;
                                elseif h ==6
                                Y_for_6(:,:,j) = Y_h;
                            end
                    end
                end                
clear ('Y_forh', 'Y_h', 'Yh', 'X_s', 'Y_s');
         if lambda~=0 
%--------------------------------Calibração-------------------------------%
        sigma = 1; %Silveira (2008)
        etta = 2; % Galí and Monacelli (2005)
        gama = 1; % Galí and Monacelli (2005)
        varphi = 0.25; % Silveira (2008)
        theta = 0.5; 
        %betta  = 0.989; % Silveira (2008)
        alppha = 0.12; % Importações/PIB
        % Measurement Equation
        gamma_s = 0.1;
        q_s = -0.1;
        pi_s = 0.5;
        rho = 1;
        % Monetary Policy Parameters 
        phi_pi = 1.75; % Silveira (2008)
        phi_y = 0.5; % Silveira (2008)
        % AR(1) Process
        rho_a = 0.5; % Galí and Monacelli (2005)                                                               
        rho_y = 0.5; % Galí and Monacelli (2005)
        rho_z = 0.5; % Silveira (2008)
        % Shock 
        sig_eps_y = 0.01; % Silveira (2008)
        sig_eps_nu = 0.01; % Silveira (2008)
        sig_eps_z = 0.01; % Silveira (2008)
        sig_eps_a = 0.01; % Silveira (2008)
        % Relations
        betta           = 1/(1+rho);
        omega           = sigma*gama+(1-alppha)*(sigma*etta-1);
        sigma_alppha    = sigma/((1-alppha)+alppha*omega);
        Theta           =(sigma*gama-1)+(1-alppha)*(sigma*etta-1);
        lambbda          = (1-(betta*theta))*(1-theta)/theta;
        kappa           = lambbda*(sigma_alppha+varphi);
        Gamma           = (1+varphi)/(sigma_alppha+varphi);
        Psi             = -Theta*sigma_alppha/(sigma_alppha+varphi);

%-----------------------------Parâmetros DSGE-----------------------------%
        theta_est(1,1) = sig_eps_a;
        theta_est(2,1) = sig_eps_y;
        theta_est(3,1) = sig_eps_z;
        theta_est(4,1) = sig_eps_nu;
        theta_est(5,1)  = sigma;
        theta_est(6,1)  = varphi;
        theta_est(7,1)  = theta;
        theta_est(8,1)  = phi_pi;
        theta_est(9,1)  = phi_y;
        theta_est(10,1)  = rho_a;
        theta_est(11,1)  = rho_y;
        theta_est(12,1)  = rho_z;
        theta_est(13,1)  = gamma_s;
        theta_est(14,1)  = pi_s;
        theta_est(15,1)  = q_s;
        theta_est(16,1)  = rho;

        pp = size(theta_est,1); %Número de parâmetros a serem estimados 

%-----------------------------------MCMC----------------------------------%
        S_total = 2000;
        burnin = 0.5;
        S0 = burnin*S_total; % burnin
        S1 = S_total-S0; %Número de iterações
        nm = 0; %contador

%----------------Matriz de Covariância da Densidade Proposta--------------%
        objective_function_penalty_base = 100000000;    
        [bounds] = priordensity(prior.lower, prior.upper, prior.pshape, prior.mean, prior.var);
        LB = bounds(:,1);
        UB = bounds(:,2);

        opt = optimset('display','iter','MaxFunEvals',1000000,'MaxIter',60000,'TolFun',1e-8,'TolX',1e-6);
        c = 0.05;    
        rng(0);
        [fval,theta_real,grad,hessian_csminwel,itct,fcount,retcodeh] = csminwel1(@(x)dsgevar_posterior(x, alppha, gama, etta, Y, p,lambda,k,X,n,prior.mean,prior.var),prior.mean,1e-4*eye(pp,pp),[],1e-7,1000,2,1e-6);
        hh = reshape(hessian(@(x)dsgevar_posterior(x, alppha, gama, etta, Y, p,lambda,k,X,n,prior.mean,prior.var),theta_real, [0.01;1]),pp,pp);  
        invHess = inv(hh);
        V = c*(nearestSPD(invHess));    

        %Pré alocando 
        Y_one = nan(n,h,S1);
        theta_s = nan(S1, pp); 
        Phi_star_s = nan(k,n,S1); %Equation (22) 
        Sigma_star_s = nan(n,n,S1); %Equation (23) 
        Phi_tilde_s = nan(k,n,S1); %Equation (28)
        Sigma_tilde_s = nan(n,n,S1); %Equation (29)
        Sigma_u_s = nan(n,n,S1); %Equation (24) - prior DSGE-VAR
        phi_s = nan(k,n,S1); %Equation (25) - prior DSGE-VAR
        Sigma_u_post_s = nan(n,n,S1); %Equation (30) - posterior DSGE-VAR
        phi_post_s = nan(k,n,S1); %Equation (31) - posterior DSGE-VAR
        logPost_s= nan(S1,1);% log posterior

        npd = 0;
        it = 0; %Contador do número de iterações
       
while it<=S0+S1
     
%---------------------------Densidade Proposta----------------------------%

theta_prop = theta_real + mvnrnd(zeros(1,pp),V)';  

if all(theta_prop > LB) && all(theta_prop < UB) || theta_prop(8,1)>1
    [npost_p, ~, ~, exit_flag, logLikeValue_p, logPrior_p, sigma_u_star_p, phi_star_p, Sigma_tilde_p, Phi_tilde_p, GXX_p, mXX_p, Z_p, R_p, T_p] = dsgevar_posterior(theta_prop, alppha, gama, etta, Y, p, lambda,k,X,n, prior.mean, prior.var);
    if exit_flag == 0
    clear('theta_prop');
    theta_prop = theta_real;
    continue;
    end   
else
clear('theta_prop');
theta_prop = theta_real;
continue;
end

%-------------------------------------------------------------------------%

[npost,  ~, ~, exit_flag, logLikeValue, logPrior, sigma_u_star, phi_star, Sigma_tilde, Phi_tilde, GXX, mXX, Z, R, T] = dsgevar_posterior(theta_real, alppha, gama, etta, Y, p, lambda,k,X,n, prior.mean, prior.var);

%---------------------------Probability of move---------------------------%

MHprob = exp((-npost_p)-(-npost));

%----------------------------Taxa de Aceitação----------------------------%
if rand<MHprob % atualiza o estado com probabilidade MHprop
npost = npost_p;
theta_real = theta_prop;
sigma_u_star = sigma_u_star_p;
phi_star = phi_star_p;
Sigma_tilde = Sigma_tilde_p;
Phi_tilde = Phi_tilde_p;
GXX = GXX_p;     
nm = nm+1; %atualiza o contador 
end

if it>S0 %burn-in
theta_s(it-S0,:) = theta_real'; %deep parameters 
Phi_star_s(:,:,it-S0) = phi_star; %Equation (22)
Sigma_star_s(:,:,it-S0) = sigma_u_star; %Equation (23)
Phi_tilde_s(:,:,it-S0) = Phi_tilde; %Equation (28)    
Sigma_tilde_s(:,:,it-S0) = Sigma_tilde; %Equation (29)
logPost_s(it-S0,:) = -npost;
end

%--------------------------------DSGE-VAR---------------------------------%
%---------------------------------Priori----------------------------------%
%Parâmetros da distribuição Wishart Invertida
s0 = lambda*t*sigma_u_star; %matriz de covariância 
v0 = lambda*t-k; %graus de liberdade
[~,g] = cholcov(s0);
    if g ~=0
    s0 = nearestSPD(s0);
    npd = npd+1;
    end;
Sigma_u = iwishrnd(s0,v0); %sigma_u | theta

%Parâmetros da distribuição Normal
mu = vec(phi_star); %média
Sigma = kron(Sigma_u, inv(lambda*t*GXX)); %matriz de covariância
[~,g] = cholcov(Sigma);
    if g ~=0
    Sigma = nearestSPD(Sigma);
    npd = npd+1;
    end;
phi = mvnrnd(mu, Sigma); %phi| Sigma_u, theta
phi = reshape(phi,k,n);

%-------------------------------Posteriori--------------------------------%
%Parâmetros da distribuição Wishart Invertida
s0_post = (lambda+1)*t*Sigma_tilde; %matriz de covariância
v0_post = (lambda+1)*t-k; %graus de liberdade     
[~,g] = cholcov(s0_post);
    if g ~=0
    s0_post = nearestSPD(s0_post);
    npd = npd+1;
    end;
Sigma_u_post = iwishrnd(s0_post,v0_post); %sigma_u | theta

%Parâmetros da distribuição Normal
mu_post = vec(Phi_tilde); %média
Sigma_post = kron(Sigma_u_post, inv(lambda*t*GXX+mXX)); %matriz de covariância
[~,g] = cholcov(Sigma_post);
    if g ~=0
    Sigma_post = nearestSPD(Sigma_post);
    npd = npd+1;
    end;
phi_post = mvnrnd(mu_post, Sigma_post); %phi| sigma_u, theta
phi_post = reshape(phi_post,k,n);

%--------------------------------Previsões--------------------------------%
if h ==1
[Y_s] = forecast_one_step(Y, Phi_tilde, p);
else
        if p ==1
        X_s = [Y(end,:) ones(1,1)];   
        elseif p ==2
        X_s = [Y(end,:) Y(end-1,:) ones(1,1)];    
        else
        X_s = [Y(end,:) Y(end-1,:) Y(end-2,:) ones(1,1)];    
        end
        
    for jj=1:h
    Yh = X_s*phi_post+mvnrnd(zeros(1,n),Sigma_u_post);

        if p == 1
        X_s = [Yh, ones(1,1)];
        elseif p ==2
        S = X_s(:,1:4);
        X_s = [Yh, S, ones(1,1)];
        else
        S = X_s(:,1:8);
        X_s = [Yh, S, ones(1,1)];    
        end
    Y_for_h(:,jj) = Yh;
    end
        if h==3
        Y_for3 = Y_for_h(:,end);
            elseif h==6
            Y_for6 = Y_for_h(:,end);
        end
        clear('Y_for_h');
end
    if it>S0
        if h==3
            Yhh3(:,:,it-S0) = Y_for3;
            Yhh_3 = mean(Yhh3, 3);
        elseif h==6
            Yhh6(:,:,it-S0) = Y_for6;
            Yhh_6 = mean(Yhh6, 3);
        else 
            Yone(:,:,it-S0) = Y_s;
            Y_one_m = mean(Yone, 3);
        end
phi_s(:,:,it-S0) = phi;   %Equation (24)
Sigma_u_s(:,:,it-S0) = Sigma_u; %Equation (25)
phi_post_s(:,:,it-S0) = phi_post; %Equation (29)
Sigma_u_post_s(:,:,it-S0) = Sigma_u_post; %Equation (30)
Sigma = mean(Sigma_u_post_s,3);
Phi = mean(phi_post_s,3);
thetas = mean(theta_s);
thetavar = var(theta_s);
    end  
    it = it+1;
end
ar_t = nm/(S1+S0); %acceptance rate

[pvalor] = convergence(theta_s, S1);
fprintf('P-valor para a hipótese nula de convergência da cadeia de Markov%1.2f\n',pvalor)
pvalor_lambda(:,j,index_lambda) = pvalor;
    
    if h == 1
    Y_for(:,:,j) = Y_one_m;
        elseif h ==3
        Y_for_3(:,:,j) = Yhh_3;
        elseif h==6
        Y_for_6(:,:,j) = Yhh_6;    
    end
        end
            end        
        end
    end
yforlambda(:,:,:,index_lambda)=Y_for; %(n,h-steps-ahead,t,lambda)
yforlambda3(:,:,:,index_lambda)=Y_for_3; %(n,h-steps-ahead,t,lambda)
yforlambda6(:,:,:,index_lambda)=Y_for_6; %(h-steps-ahead,n,t,lambda)
clear ('Y_for', 'Y_for_3', 'Y_for_6'); 
end
%%
%--------------------------------Gráficos---------------------------------%
%h=1 
% for i=1:index_lambda
% figure(i),    
%     for j=1:n
% subplot(2,2,j)
% yforecast1 = yforlambda(j,:,:,i);
% yforecast1 = yforecast1(:,:)';
% plot(Data(end-kk+p:end),[YY(end-kk+p:end,j) yforecast1]) 
%     end
% end
% %h=3
% for i=1:index_lambda
% figure(i),    
%     for j=1:n
% subplot(2,2,j)
% yforecast3 = yforlambda3(j,:,:,i);
% yforecast3 = yforecast3(:,:)';
% plot(Data(end-kk+p:end),[YY(end-kk+p:end,j) yforecast3]) 
%     end
% end
% %h=6
% for i=1:index_lambda
% figure(i),    
%     for j=1:n
% subplot(2,2,j)
% yforecast6 = yforlambda6(j,:,:,i);
% yforecast6 = yforecast6(:,:)';
% plot(Data(end-kk+p:end),[YY(end-kk+p:end,j) yforecast6]) 
%     end
% end
%%
%------------------------------Função Perda-------------------------------%
h = 1;
for i=1:index_lambda
    for j=1:n
yforecast1 = yforlambda(j,:,:,i);
yforecast1 = yforecast1(:,:)';
L = (YY(end-kk+h+1:end,j) - yforecast1).^2;
L1(:,j,i) = L;
    end
end

h = 3;
for i=1:index_lambda
    for j=1:n
yforecast3 = yforlambda3(j,:,:,i);
yforecast3 = yforecast3(:,:)';             
L_3 = (YY(end-kk+h+1:end,j) - yforecast3).^2;
L3(:,j,i) = L_3;
    end
end
h = 6;
for i=1:index_lambda
    for j=1:n
yforecast6 = yforlambda6(j,:,:,i);
forecast6 = yforecast6(:,:)';
L_6 = (YY(end-kk+h+1:end,j) - forecast6).^2;
L6(:,j,i) = L_6;
    end
end
%-------------------------------------------------------------------------%

