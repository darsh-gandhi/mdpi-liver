function [lambda_l,lambda_c,k_l,k_c,gamma_l,gamma_c,nu,errPercentage] = liver_MultiStart(data,initialvalues,tspan,lb,ub,w,alpha,NoStartPoints,dataChoice)
% liver_MultiStart
% Authored by Emma and Darsh, adopted from code for a simple SIR model by Christina Edholm.
% 
% This code fits the OROV epidemic model formulated by the Biomath class to
% weekly infectious data. The code calls orovModel.m, where the model is
% "initialized" in MATLAB, and orov_run_odeSolver.m, where a numerical solver
% is implemented. It will also plot results.
% 
% Inputs: 
%   data            <- reported weekly new infectious cases in Amazonas, Brazil
%   initialvalues   <- intial values for each state variable 
%   tspan           <- timespan of experiment
%   lb              <- lower bounds for the parameters to be fit
%   ub              <- upper bounds for the parameters to be fit
%   w               <- omega (data-fitting weight)
%   alpha           <- alpha (regularization weight)
%   NoStartPoints   <- number of runs for the MultiStart algorithm
%   dataChoice      <- choice of dataset for data fitting 
% Outputs:
%   []              <- this vector contains the "optimal" values of unknown
%                       parameters and the relative l2 error percentage

%% Upper and Lower bounds for parameters you are fitting.

%Parameter vector we are approximating:
% z = [(1) lambda_l, (2) lambda_c, (3) k_l, (4) k_c, (5) gamma_l, (6) gamma_c, (7) nu]

%What initial parameter values you want to start with
if dataChoice == 0
    xstart=.5*(lb+ub);   % midpoint is the initial guess for fitting to liver only data
    % xstart = rand(size(UpperBounds)).*(UpperBounds-LowerBounds)+LowerBounds;        % random values w/n the bounds for every parameter
else                    % randomizing initial guess for cancer dataset
    xstart=.5*(lb+ub);
    fixed_param_indices = [1 3 5 6 7];
    xstart(fixed_param_indices) = .5*(lb(fixed_param_indices)+ub(fixed_param_indices));
    xstart([2 4]) = rand([1 2]).*(ub([2 4])-lb([2 4]))+lb([2 4]);
    xstart
end

%% MultiStart and fmincon - Fitting Part - Parallelization

%This section of the code develops and solves the optimization problem. We
% minimize the cost functional over the parameters (z) in the space between
% vectors lb and ub, with the value of the functional measured in the
% liver_Functional function. We use given initial conditions and bounds.

% setting up the optimization problem
problem = createOptimProblem('fmincon','objective',@(z) liver_Functional(z,data,initialvalues,tspan,w,alpha,dataChoice)...
    ,'x0',xstart,'lb',lb,'ub',ub);%,'Aineq',A,'bineq',b)%,'Aeq',Aeq,'beq',beq);

problem.options = optimoptions(problem.options,'MaxFunEvals',9999,'MaxIter',9999); %,'TolCon',0)

ms=MultiStart('UseParallel',false,'Display','iter');        % defines a multistart problem
[b,fval,exitflag,output,manymins]=run(ms,problem,NoStartPoints);  % runs the multistart

% the following takes solutions from manymins and makes a matrix out of them

for i=1:length(manymins)
    parameters(i,:)=manymins(i).X;  %what are the parameter values
end

for i=1:length(manymins)
    fvalues(i)=manymins(i).Fval;    %the minimization error
end

for i=1:length(manymins)
    ExitFlags(i)=manymins(i).Exitflag; %how "good" is the solution, we want 1 or 2.
end


%% Plot the "best" solution

%Outputs state variables for "best" fit
[~,y] = ode23s(@(t,y) dimensional_odes(t,y,parameters(1,:)), tspan, initialvalues);

if dataChoice == 0    % fitting to liver only data
    diff = y(:,1) - reshape(data,size(y(:,1)));
elseif dataChoice == 1 % fitting to cancer only data
    diff = y(:,2) - reshape(data,size(y(:,2)));
else                                       % fitting to liver only data
    diff = y(:,1) - reshape(data,size(y(:,1)));
end

errPercentage = sqrt(sum(diff.^2)./sum((data.^2)))*100; % relative l2 error percentage

%% define parameters to send back

lambda_l=parameters(1,1);
lambda_c=parameters(1,2);
k_l = parameters(1,3);
k_c = parameters(1,4);
gamma_l=parameters(1,5);
gamma_c=parameters(1,6);
nu = parameters(1,7);