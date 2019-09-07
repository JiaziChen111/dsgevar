close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Calibration %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
betta          = 1/(1+rho);
omega = sigma*gama+(1-alppha)*(sigma*etta-1);
sigma_alppha =sigma/((1-alppha)+alppha*omega);
Theta=(sigma*gama-1)+(1-alppha)*(sigma*etta-1);
lambda = (1-(betta*theta))*(1-theta)/theta;
kappa =lambda*(sigma_alppha+varphi);
Gamma = (1+varphi)/(sigma_alppha+varphi);
Psi = -Theta*sigma_alppha/(sigma_alppha+varphi);

%%%%%%%%%% Variable ID assignments %%%%%%%%%%%
    var_y	    = 10;    % Output
    var_pi	    = 11;	% Inflation
    var_i	    = 1;	% Interest Rate
    var_q       = 6;    % Exchange Rate
    var_s_y 	= 7;	% AR(1) process: World Production Shock
    var_s_a  	= 8;    % AR(1) process: Technology Shock
    var_s_z     = 9;    % AR(1) process: Preference Shock  
    var_y_obs   = 2;
    var_pi_obs  = 3;
    var_i_obs   = 4;
    var_q_obs   = 5;

    % Shock ID assignment
    e_i         = 4;	% Monetary Policy Shock
    e_y         = 2;    % World Output Shock
    e_a         = 1;	% Technology Shock 
    e_z         = 3;	% Preference Shock
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
    VC(4,4) = sig_eps_nu*sig_eps_nu; %  productivity shock
    VC(1,1) = sig_eps_a*sig_eps_a; %  world output shock
    VC(2,2) = sig_eps_y*sig_eps_y; %  monetary policy shock
    VC(3,3) = sig_eps_z*sig_eps_z; %  preference shock

    %%%%%%%%%% MODEL EQUATIONS %%%%%%%%%%
    EqNb 				= 0; % starting the counter;

    % Output
    EqNb                    = EqNb+1; % equation number
    AA(EqNb, var_y)         = -1;
    AA(EqNb, var_pi)        = -1/sigma_alppha;
    BB(EqNb, var_i)         = 1/sigma_alppha;
    BB(EqNb, var_s_y)       = -alppha*(omega-1)*(rho_y-1);
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
    AA(EqNb, var_s_y)       = -phi_y*alppha*Psi;
    AA(EqNb, var_s_a)       = -phi_y*Gamma;
    AA(EqNb, var_pi)        = phi_pi;
    AA(EqNb, var_y)         = phi_y;
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
   
    [T,~,R,~,~,~,~,eu,~]=gensys(-AA,BB,CO,CC,DD); %Solve the DSGE model
   