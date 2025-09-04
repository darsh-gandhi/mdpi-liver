function [lL,lC,kL,kC,gL,gC,nu,l0,c0,t0,tf,Nt,t,dt] = parameter_control()
%parameters := stores the parameters and initial conditions for the system of ODEs
%   Inputs: no args
%   Outputs:
%       lL <- lambda_L, healthy proliferation rate
%       lC <- lambda_C, cancer proliferation rate
%       kL <- liver carrying capacity
%       kC <- cancer carrying capacity
%       gL <- gamma_L
%       gC <- gamma_C
%       v  <- nu
%       l0 <- initial liver
%       c0 <- initial cancer
%       t0 <- initial time
%       tf <- final time
%       Nt <- # of time points in discretization
%       t  <- time grid
%       dt <- step size (in time grid)


%params from data-fitting (new)
lL = 2.16918177423846;
lC = 1.620612057098;
kL=1500;
kC=7500;
gL=0.142569012229;
gC=0;
nu = 0;

%initial conditions
l0 = 500;  %mm^3
c0 = 50;   %mm^3

% c_crit > l_crit <-> 1-1/gC > 1-gL <-> gLgC>1
% c_crit < l_crit <-> 1-1/gC < 1-gL <-> gLgC<1

% 0 <= nu <= lambda_C


% for coexistence: lC*(1-1/gC) < nu < lC*(1-gL), and 
% if nu>lC*((gC-gL)/(1+gC)), L^*>C^* (in coexistence)

t0 = 0; % initial time (weeks)
tf = 25; % final time (weeks)

Nt = 1001; % Number of time grid points
t = linspace(t0,tf,Nt)'; % Time grid
dt = t(2)-t(1); % Mesh size.
end