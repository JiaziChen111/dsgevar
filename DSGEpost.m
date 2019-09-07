function [nlog_post] = DSGEpost(theta_prop, betta, alppha, gama, etta, Y, p)

pp = size(theta_prop,1);
THETA(1:pp,1)= theta_prop;

sigma = THETA(1,1);
varphi = THETA(2,1);
theta = THETA(3,1);
phi_pi = THETA(4,1);
phi_y = THETA(5,1);
phi_i = THETA(6,1);
rho_a = THETA(7,1);
rho_y = THETA(8,1);
rho_z = THETA(9,1);
gamma_s = THETA(10,1);
pi_s = THETA(11,1);
q_s = THETA(12,1);
sig_eps_a = THETA(13,1);
sig_eps_y = THETA(14,1);
sig_eps_nu = THETA(15,1);
sig_eps_z = THETA(16,1);
% Relations
rho            = betta^(-1)-1;
lambbda         = ((1-betta*theta)*(1-theta))/theta; %p. 717
omega          = sigma*gama+(1-alppha)*(sigma*etta-1); %p.716
sigma_alppha   = sigma/((1-alppha)+alppha*omega); %p. 717
Theta          = omega-1; %p.717
kappa          = lambbda*(sigma_alppha+varphi); %p.718
Gamma          = (1+varphi)/(sigma_alppha+varphi); %p.718
Psi            = -(Theta*sigma_alppha)/(sigma_alppha+varphi); %p.718
THETA(17,1) = betta;
THETA(18,1) = alppha;
THETA(19,1) = gama;
THETA(20,1) = etta;
THETA(21,1) = rho;
THETA(22,1) = lambbda;
THETA(23,1) = omega;
THETA(24,1) = sigma_alppha;
THETA(25,1) = Theta;
THETA(26,1) = kappa;
THETA(27,1) = Gamma;
THETA(28,1) = Psi;

%%%%%%%%%% Variable ID assignments %%%%%%%%%%%
    var_y	    = 10;    % Output
    var_pi	    = 11;	% Inflation
    var_i	    = 5;	% Interest Rate
    var_q       = 6;    % Exchange Rate
    var_s_y 	= 7;	% AR(1) process: World Production Shock
    var_s_a  	= 8;    % AR(1) process: Technology Shock
    var_s_z     = 9;    % AR(1) process: Preference Shock  
    var_y_obs    = 1;
    var_pi_obs  = 2;
    var_i_obs   = 4;
    var_q_obs   = 3;

    % Shock ID assignment
    e_i         = 1;	% Monetary Policy Shock
    e_y         = 3;    % World Output Shock
    e_a         = 2;	% Technology Shock 
    e_z         = 4;	% Preference Shock
    % Expectation Erros
    eta_e       = 1;
    eta_ee      = 2;
    
    %%%%%%%%%% MATRIX DEFINITIONS %%%%%%%%%%
    % Number of endogenous variables
    nx = 11; 
    % Number of shocks
    nu = 4; 
    % Declaring matrix
    AA = zeros(nx,nx); %endogenous variables
    CO = zeros(nx,1); %
    BB = zeros(nx,nx); %lagged endogenous variables
    CC = zeros(nx,nu); %exogenous errors
    VC = zeros(nu,nu); %variance-covariance
    DD = zeros(nx,2); %expectation errors 

    %%%%%%%%%% VARIANCE-COVARIANCE MATRIX %%%%%%%%%%
    VC(1,1) = sig_eps_nu*sig_eps_nu; %  productivity shock
    VC(2,2) = sig_eps_a*sig_eps_a; %  world output shock
    VC(3,3) = sig_eps_y*sig_eps_y; %  monetary policy shock
    VC(4,4) = sig_eps_z*sig_eps_z; %  preference shock

    %%%%%%%%%% MODEL EQUATIONS %%%%%%%%%%
    EqNb 				= 0; % starting the counter;

    % Output
    EqNb                    = EqNb+1; % equation number
    AA(EqNb, var_y)         = -1;
    AA(EqNb, var_pi)        = -1/sigma_alppha;
    BB(EqNb, var_i)         = 1/sigma_alppha;
    BB(EqNb, var_s_y)       = -alppha*Theta*(rho_y-1);
    BB(EqNb, var_s_z)       = -(1/sigma_alppha)*(1-rho_z);
    BB(EqNb, var_y)         = 1;
    DD(EqNb, eta_e)         = 1;
    DD(EqNb, eta_ee)        = 1;
    % NK Phillips Curve 
    EqNb                    = EqNb+1; % equation number
    AA(EqNb, var_pi)        = -betta;
    BB(EqNb, var_s_a)       = kappa*Gamma;
    BB(EqNb, var_s_y)       = kappa*alppha*Psi;
    BB(EqNb, var_y)         = -kappa;
    BB(EqNb, var_pi)        = 1; 
    DD(EqNb, eta_ee)        = 1;
    % Open Economy
    EqNb                    = EqNb+1; % equation number
    AA(EqNb, var_q)         = -1;
    AA(EqNb, var_s_y)       = -sigma_alppha*(1-alppha);
    AA(EqNb, var_y)         = sigma_alppha*(1-alppha);
    % Monetary Policy
    EqNb                    = EqNb+1; % equation number
    AA(EqNb, var_i)         = -1;
    AA(EqNb, var_s_y)       = -(1-phi_i)*(phi_y*alppha*Psi);
    AA(EqNb, var_s_a)       = -(1-phi_i)*(phi_y*Gamma);
    BB(EqNb, var_i)         = phi_i;
    AA(EqNb, var_pi)        = (1-phi_i)*(phi_pi);
    AA(EqNb, var_y)         = (1-phi_i)*(phi_y);
    CC(EqNb, e_i)           = 1;
    % AR(1) process: World Output shock
    EqNb                    = EqNb+1; % equation number
    AA(EqNb, var_s_y)       = -1;
    BB(EqNb, var_s_y)       = rho_y;
    CC(EqNb, e_y)           = 1;
    % AR(1) process: Technology shock
    EqNb                    = EqNb+1; % equation number
    AA(EqNb, var_s_a)       = -1;
    BB(EqNb, var_s_a)       = rho_a;
    CC(EqNb, e_a)           = 1;
    % AR(1) process: Preference Shock 
    EqNb                    = EqNb+1; % equation number
    AA(EqNb, var_s_z)       = -1;
    BB(EqNb, var_s_z)       = rho_z;
    CC(EqNb, e_z)           = 1; 
    
    EqNb                    = EqNb+1; % equation number
    AA(EqNb, var_y_obs)     = 1;
    AA(EqNb, var_y)         = -1;
    BB(EqNb, var_y)         = 1;
    AA(EqNb, var_s_a)       = -1;
    CO(EqNb,1)              = gamma_s;
    
    EqNb                    = EqNb+1; % equation number
    AA(EqNb, var_pi_obs)    = 1;
    AA(EqNb, var_pi)        = -1;
    CO(EqNb,1)              = pi_s;
    
    EqNb                    = EqNb+1; % equation number
    AA(EqNb, var_i_obs)     = 1;
    AA(EqNb, var_i)        = -12;
    CO(EqNb,1)              = 12*(pi_s+rho);

    EqNb                    = EqNb+1; % equation number
    AA(EqNb, var_q_obs)     = 1;
    AA(EqNb, var_q)         = -1;
    BB(EqNb, var_q)         = 1;
    CO(EqNb,1)              = q_s;
    
    [T,~,R,~,~,~,~,~,~]=gensys(-AA,BB,CO,CC,DD); %Solve the DSGE model
    
    % Measurement Equation 
    EqNb    				= 0; % starting the counter;
    Z = zeros(nu, nx);
    D = zeros(nu, 1); 
    % Output
    EqNb                    = EqNb+1; % equation number
    Z(EqNb, var_y)          = 1;
    Z(EqNb, var_s_a)        = 1;
    % Inflation
    EqNb                    = EqNb+1; % equation number
    Z(EqNb, var_pi)         = 1;
    % Exchange Rate
    EqNb                    = EqNb+1; % equation number
    Z(EqNb, var_q)          = 1;
    % Interest Rate 
    EqNb                    = EqNb+1; % equation number
    Z(EqNb, var_i)          = 12;

    D(1, 1)             = gamma_s;
    D(2, 1)             = pi_s;
    D(4, 1)             = 12*(pi_s+rho);
    D(3, 1)             = q_s;
    
    [T1, ~] = size(T);
    KsiLast = zeros(T1, 1);
    t = size(Y,1);
    XX = ones(1,t);
    Q = R*VC*transpose(R);
    prior_SS = [0, 0, 0, 0.1335];
    trend = repmat(prior_SS,168-(p+1),1);
    Y = Y-trend;
    %Kalman Filter
    PLast = lyapunov_symm(T,Q,1.000,1.000e-15,0);
    last = length(Y);
    a = zeros(11,1);
    ZZ = [1;2;3;4];
    [lnL] = kalman_filter(Y',1, last, a,PLast, 1.0000e-6, 1.0000e-6, 0, T, VC, R, 0, ZZ, 11, 4, 4, 0, 0); 
    x = [1; 1; 1; 1; 1; 0.25; 0.5; 1.01; 0.25; 0.5; 0.8; 0.8; 0.5; 0;0;0];
    lower = [0;0;0;0;0;0;0;0;0;0;0;0;0; -inf; -inf; -inf];
    upper = [inf; inf; inf; inf; inf; inf; 1; inf; inf;1;1;1;1;inf;inf;inf];
    prior_1 = [0.6366; 0.6366; 0.6366; 0.6366; 11.1111; 0.6944; 6; 11.3344; 0.6944; 12; 12; 12; 12; 0.5; 0.5; 0.5];
    prior_2 = [2;2;2;2;0.09;0.36;14;0.0891; 0.36; 12;12;12; 12; 0.2; 0.2; 0.2];
    pshape = [4;4;4;4;2;2;1;2;2;1;1;1;1;3;3;3];
    logPrior = priordens(x,pshape, prior_1, prior_2, lower, upper,1);
     
    nlog_post = -(-lnL-logPrior); 
end

