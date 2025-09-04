% runner.m := this file contains data and runs the parameter estimation

clear; close all; format LONG
%%

% we are estimating parameters in a DIMENSIONAL model, where the healthy
% liver volume is 1500 mm^3 and cancer carrying capacity (k) is estimated.

% params:= [(1) lambda_l, (2) lambda_c, (3) k_l, (4) k_c, (5) gamma_l, (6) gamma_c, (7) nu]
lb=[0.1,0.1,1000,7500,0,1,0];             % lower bounds for parameter estimates
ub=[10,5,2000,7500,1,2,1]; % upper bounds for parameter estimates

dataChoice = input('Which dataset are you using?\n (0) Liver only (1) Tumor (Anton M1) \n');

% dataChoice == 0 <- Koniaris liver regeneration data (already nondimensional)
% dataChoice == 1 <- Anton M1 mouse tumor growth (exponential)

if dataChoice == 0 %fitting lambda_l only (k_l = 1)
    timespan = 0:2:12; %(days)
    timespan = timespan/7; %(weeks)
    data = [0.33 0.41 0.54 0.726 0.995 1.105 1]*1500; % mm^3
    error_bars = [0	0.02 0.02 0.026 0.05 0.105 0.09]*1500; % mm^3
    l_0 = data(1);
    c_0 = 0;
    % ub(2)=0; %lambda_c = 0
    lb(2:4)=1;
    ub(2:4)=1;
    lb(3)=1500;
    ub(3)=1500;
    lb(5:7)=0;
    ub(5:7)=0;
elseif dataChoice == 1 % fitting lambda_c, gamma_l
    timespan = [0 3 17 21 24 28 31]; % timespan of the experiment (days)
    timespan = timespan/7; %(weeks)
    data = 1000*[3.3/1000 0.2356114635 0.2520757125 0.5022632175 0.9649201295 2.011033185 3.435962131]; % mm^3
    l_0 = 1500; % inital liver population
    c_0 = data(1); % inital liver population
    lb(1)=2.16918177423846;ub(1)=2.16918177423846; % known (estimated) lambda_l values (wks)
    lb(5:7)=0;ub(5:7)=0; % not fitting gamma_c or nu
end

plot(timespan,data,'r-o','LineWidth',1.5)
if dataChoice == 0
    hold on
    errorbar(timespan,data,error_bars,'.r','LineWidth',1.25);
end
xlabel('Days since PH')
ylabel('Fraction of Original Liver')
set(gca,'FontSize',14)

init=[l_0, c_0, l_0+c_0]; % initial conditions

n=15;        %trials of the multistart
runs=100;  %runs of fmincon per multistart trial
paramsAndOmega=-1.0*ones(n,length(lb)+1); % matrix containing "best" parameter values for each trial of Multistart, relative l2 error, and the omega value 

% paramsAndOmega(:,end) = 10.^randi([-3 4],1,n); % makes omega (data-fitting weight) value random 10^k, where k is an integer between -3 and 4
w=100; % data-fitting weight
alpha=1;% regularization weight

%% running multistart parameter estimation

inittime = cputime; % start time of trials
for i=0:n-1 % n trials
    i+1 % trial number
    %the next line runs the multistart and saves the parameter outputs and error to the paramsAndOmega matrix
    [paramsAndOmega(i+1,1), paramsAndOmega(i+1,2), paramsAndOmega(i+1,3), paramsAndOmega(i+1,4), paramsAndOmega(i+1,5), paramsAndOmega(i+1,6), paramsAndOmega(i+1,7), paramsAndOmega(i+1,8)] = liver_MultiStart(data,init,timespan,lb,ub,w,alpha,runs(mod(i,length(runs))+1),dataChoice);
    paramsAndOmega(i+1,5)=log(2)/(3*paramsAndOmega(i+1,2));
end
eval_time = cputime-inittime; % total time of trials

%%

[minErr, minRow] = min(paramsAndOmega(:,length(paramsAndOmega(1,:)))) % finding the ms trial with least error

params=paramsAndOmega(minRow,1:length(paramsAndOmega(minRow,:))) % params corresponding to least error
% params(1)=2.16918177423846;

colors = 1/255*[0 0 255; 255 0 0; 160 32 240]; % blue red purple

fineMesh = linspace(0,timespan(end),1001); % finer mesh for pretty pictures
[t,y]=ode23s(@(t,y) dimensional_odes(t,y,params), fineMesh, init); % solving the ode for plotting

%%

figure(2) % plots liver cell population vs real data (if any)
plot(fineMesh, y(:,1),'LineWidth',1.5,'Color',colors(1,:))
xlim([0 fineMesh(end)])
xlabel('Time (weeks)')
ylabel('Relative Liver Volume')
if dataChoice == 0
    hold on
    scatter(timespan,data,'filled','r')
    hold on
    errorbar(timespan,data,error_bars,'.r','LineWidth',1.25);
    xlabel('Time (weeks)')
    ylabel('Fraction of Original Liver')
    legend('Simulated', 'Observed',Location='northwest')
    hold off
    ylim([0 1.25])
end
set(gca,'FontSize',14)

%%

figure(3)
plot(fineMesh, y(:,2),'LineWidth',1.5,'Color',colors(2,:))
xlim([0 fineMesh(end)])
xlabel('Time (weeks)')
ylabel('Volume (mm^3)')
if dataChoice == 1
    hold on
    scatter(timespan,data,'filled','r')
    legend('Simulated', 'Observed',Location='northwest')
    hold off
end
set(gca,'FontSize',14)

%%

[lL,lC,kL,kC,gL,gC,nu,l0,c0,t0,tf,Nt,t,~] = parameter_control();
new_params = [lL,lC,kL,kC,gL,gC,nu];
init = [l0 c0 l0+c0];

longMesh = fineMesh*3;
[t2,y2] = ode23s(@(t,y) odes_after_fitting(t,y,new_params), longMesh, init);

figure(5)
plot(longMesh, y2(:,1),'LineWidth',1.5,'Color',colors(1,:))
hold on
plot(longMesh, y2(:,2),'LineWidth',1.5,'Color',colors(2,:))
xlim([0 longMesh(end)])
xlabel('Time (weeks)')
ylabel('Volume (mm^3)')
legend('Healthy','Cancer','Location','northwest')
set(gca,'FontSize',14)