% make_figs.m := this file creates figures***

clear; close all; format LONG
%% gen figs 1 (a-d)

lambdaL = 1;
lambdaC = 1;
Kl = 1;
Kc = 2;
mu = 0;
fig_index = 1;

for gammaL = [0.45, 5]
    for gammaC = [0.45, 5]
        % Initialize Model
        odefun = @(L, C) [lambdaL.*(1-(L+gammaC*C)/Kl).*L; lambdaC.*(1-mu)*(1-(gammaL*L+C)/((1-mu)*Kc)).*C];
        centers = [0, 0; Kl, 0; 0, Kc*(1-mu); (Kl-gammaC*Kc*(1-mu))/(1-gammaC*gammaL), (-Kl*gammaL+Kc*(1-mu))/(1-gammaC*gammaL)];
        odefun2 = @(t,y) [lambdaL.*(1-(y(1)+gammaC*y(2))/Kl).*y(1); lambdaC.*(1-mu)*(1-(gammaL*y(1)+y(2))/((1-mu)*Kc)).*y(2)];
        
        % Boundary of Plots
        max_Ecol1 = 1.1*max(centers(:,1));
        max_Ecol2 = 1.1*max(centers(:,2));
        min_Ecol1 = 0;
        min_Ecol2 = 0;
        
        % Calculate Directional Field
        index=1;
        for i=linspace(min_Ecol1,max_Ecol1,10)
            for j=linspace(min_Ecol2,max_Ecol2,10)
                k=odefun(i,j);
                X(index)=i;
                Y(index)=j;
                U(index)=k(1);
                V(index)=k(2);
                index = index+1;
            end
        end
        U = U./sqrt(U.^2+V.^2);
        V = V./sqrt(U.^2+V.^2);
        
        % Begin Figure
        figure(fig_index)
        quiver(X,Y,U,V, 'k', 'AutoScaleFactor', 0.5) % Directional Field
        hold on
        fontsize(12, "points")
        colrs = 1/255*[0 0 255; 160 32 240; 255 0 0; 255 165 0; 0 255 255; 0 0 0]; % blue purple red orange cyan black
        
        % Plot Axes
        plot(linspace(0,20,40), zeros(1, 40), 'k', 'LineWidth',1)
        plot(zeros(1, 40),linspace(0,20,40), 'k', 'LineWidth',1)
        
        % Plot Nullclines
        plot(linspace(0,20,40), zeros(1, 40), 'LineStyle','--', 'Color', colrs(1,:), 'LineWidth',2)
        plot(zeros(1, 40), linspace(0,20,40), 'LineStyle','--', 'Color', colrs(3,:), 'LineWidth',2)
        plot(linspace(0,20,40), (Kl-linspace(0,20,40))/gammaC, 'LineStyle','--', 'Color', colrs(3,:), 'LineWidth',2)
        plot(linspace(0,20,40), (Kc*(1-mu)-gammaL*linspace(0,20,40)), 'LineStyle','--', 'Color', colrs(1,:), 'LineWidth',2)
        
        % Plot Equilibria
        scatter(centers(1,1), centers(1,2), 100, "k", 'LineWidth',2) %E_0
        
        if gammaL > Kc*(1-mu)/Kl % E_L
            scatter(centers(2,1), centers(2,2), 100, "k", 'filled', 'LineWidth',2)
        else
            scatter(centers(2,1), centers(2,2), 100, "k", 'LineWidth',2)
        end
        
        
        if gammaC > Kl/(Kc*(1-mu)) % E_C
            scatter(centers(3,1), centers(3,2), 100, "k", "filled", 'LineWidth',2)
        else
            scatter(centers(3,1), centers(3,2), 100, "k", 'LineWidth',2)
        end
        
        if (gammaC > Kl/(Kc*(1-mu)) && gammaL > Kc*(1-mu)/Kl) || (gammaC < Kl/(Kc*(1-mu)) && gammaL < Kc*(1-mu)/Kl) %COE
            if gammaL*gammaC < 1 % COE Stable
                scatter(centers(4,1), centers(4,2), 100, "k", "filled", 'LineWidth',2)
                legend("Direction Field Arrows", "", "", "$$L$$ nullclines", "$$C$$ nullclines", "", "", "$$E_0, E_L, E_C$$ Unstable", "","", "$$COE$$ Stable", 'interpreter','latex')
            
            else % COE Unstable
                scatter(centers(4,1), centers(4,2), 100, "k", 'LineWidth',2)

                % Separatrix
                % sep_index = 1;
                % separatrix = [];
                % for i=linspace(min_Ecol1,max_Ecol1,1)
                %     for j=linspace(min_Ecol2,max_Ecol2,100000)
                %         tspan = linspace(0,5000,5);
                %         [t,y] = ode45(odefun2, tspan, [i,j]);
                %         if y(end,1)==centers(4,1)
                %             separatrix(sep_index, 1:2) = [i,j]
                %             sep_index = sep_index + 1;
                %         end
                %     end
                % end
                % separatrix
                % semilogy(separatrix, 'Color', colrs(2,:), 'LineWidth',2)
                legend("Direction Field Arrows", "", "","", "", "$$L$$ nullclines", "$$C$$ nullclines", "$$E_0, COE$$ Unstable", "$$E_L, E_C$$ Stable", "", "", "Separatrix", 'interpreter','latex')
            end
         
        % Other legends
        elseif (gammaL > Kc*(1-mu)/Kl && gammaC < Kl/(Kc*(1-mu))) %E_L Stable
            legend("Direction Field Arrows", "", "", "$$L$$ nullclines", "$$C$$ nullclines", "", "", "$$E_0, E_C$$ Unstable", "$$E_L$$ Stable", "",'interpreter','latex')
        elseif (gammaL < Kc*(1-mu)/Kl && gammaC > Kl/(Kc*(1-mu))) %E_C Stable
            legend("Direction Field Arrows", "", "", "$$L$$ nullclines", "$$C$$ nullclines", "", "", "$$E_0, E_L$$ Unstable", "", "$$E_C$$ Stable", "",'interpreter','latex')
        else %All Unstable
            legend("Direction Field Arrows", "", "", "$$L$$ nullclines", "$$C$$ nullclines", "", "", "$$E_0, E_L, E_C$$ Unstable", "", "", "",'interpreter','latex')
        end
        
        
        % Titles, Labels, and File Saving
        if gammaL < 1 && gammaC < 1
            plot_title = "Phase Portrait: $$\gamma_L, \gamma_C < 1, \mu=$$ " + string(mu);
        elseif gammaL > 1 && gammaC > 1
            plot_title = "Phase Portrait: $$\gamma_L, \gamma_C > 1, \mu=$$ " + string(mu);
        elseif gammaL > 1 && gammaC < 1
            plot_title = "Phase Portrait: $$\gamma_L > 1, \gamma_C < 1, \mu=$$ " + string(mu);
        elseif gammaL < 1 && gammaC > 1
            plot_title = "Phase Portrait: $$\gamma_L < 1, \gamma_C > 1, \mu=$$ " + string(mu);
        end
        % title(plot_title,'interpreter','latex')
        xlabel('Liver cells $$L$$', 'Interpreter','latex')
        xlim([0-(max(centers(:,1)))*0.025, max_Ecol1])
        ylabel('Cancer cells $$C$$', 'Interpreter','latex')
        ylim([0-(max(centers(:,2)))*0.025, max_Ecol2])
        file_title = "pp_gL" + string(gammaL) + "_gC" + string(gammaC) + "_mu" + string(mu);
        saveas(gcf,file_title+".fig",'fig')
        saveas(gcf,file_title+".eps",'epsc')
        hold off
    fig_index = fig_index + 1;
    end
end

%% gen figs 2 (a-d)

lambdaL = 1;
lambdaC = 1;
Kl = 1;
Kc = 2;
mu = 0.8;

for gammaL = [0.45, 5]
    for gammaC = [0.45, 5]
        % Initialize Model
        odefun = @(L, C) [lambdaL.*(1-(L+gammaC*C)/Kl).*L; lambdaC.*(1-mu)*(1-(gammaL*L+C)/((1-mu)*Kc)).*C];
        centers = [0, 0; Kl, 0; 0, Kc*(1-mu); (Kl-gammaC*Kc*(1-mu))/(1-gammaC*gammaL), (-Kl*gammaL+Kc*(1-mu))/(1-gammaC*gammaL)];
        odefun2 = @(t,y) [lambdaL.*(1-(y(1)+gammaC*y(2))/Kl).*y(1); lambdaC.*(1-mu)*(1-(gammaL*y(1)+y(2))/((1-mu)*Kc)).*y(2)];
        
        % Boundary of Plots
        max_Ecol1 = 1.1*max(centers(:,1));
        max_Ecol2 = 1.1*max(centers(:,2));
        min_Ecol1 = 0;
        min_Ecol2 = 0;
        
        % Calculate Directional Field
        index=1;
        for i=linspace(min_Ecol1,max_Ecol1,10)
            for j=linspace(min_Ecol2,max_Ecol2,10)
                k=odefun(i,j);
                X(index)=i;
                Y(index)=j;
                U(index)=k(1);
                V(index)=k(2);
                index = index+1;
            end
        end
        U = U./sqrt(U.^2+V.^2);
        V = V./sqrt(U.^2+V.^2);
        
        % Begin Figure
        figure(fig_index)
        quiver(X,Y,U,V, 'k', 'AutoScaleFactor', 0.5) % Directional Field
        hold on
        fontsize(12, "points")
        colrs = 1/255*[0 0 255; 160 32 240; 255 0 0; 255 165 0; 0 255 255; 0 0 0]; % blue purple red orange cyan black
        
        % Plot Axes
        plot(linspace(0,20,40), zeros(1, 40), 'k', 'LineWidth',1)
        plot(zeros(1, 40),linspace(0,20,40), 'k', 'LineWidth',1)
        
        % Plot Nullclines
        plot(linspace(0,20,40), zeros(1, 40), 'LineStyle','--', 'Color', colrs(1,:), 'LineWidth',2)
        plot(zeros(1, 40), linspace(0,20,40), 'LineStyle','--', 'Color', colrs(3,:), 'LineWidth',2)
        plot(linspace(0,20,40), (Kl-linspace(0,20,40))/gammaC, 'LineStyle','--', 'Color', colrs(3,:), 'LineWidth',2)
        plot(linspace(0,20,40), (Kc*(1-mu)-gammaL*linspace(0,20,40)), 'LineStyle','--', 'Color', colrs(1,:), 'LineWidth',2)
        
        % Plot Equilibria
        scatter(centers(1,1), centers(1,2), 100, "k", 'LineWidth',2) %E_0
        
        if gammaL > Kc*(1-mu)/Kl % E_L
            scatter(centers(2,1), centers(2,2), 100, "k", 'filled', 'LineWidth',2)
        else
            scatter(centers(2,1), centers(2,2), 100, "k", 'LineWidth',2)
        end
        
        
        if gammaC > Kl/(Kc*(1-mu)) % E_C
            scatter(centers(3,1), centers(3,2), 100, "k", "filled", 'LineWidth',2)
        else
            scatter(centers(3,1), centers(3,2), 100, "k", 'LineWidth',2)
        end
        
        if (gammaC > Kl/(Kc*(1-mu)) && gammaL > Kc*(1-mu)/Kl) || (gammaC < Kl/(Kc*(1-mu)) && gammaL < Kc*(1-mu)/Kl) %COE
            if gammaL*gammaC < 1 % COE Stable
                scatter(centers(4,1), centers(4,2), 100, "k", "filled", 'LineWidth',2)
                legend("Direction Field Arrows", "", "", "$$L$$ nullclines", "$$C$$ nullclines", "", "", "$$E_0, E_L, E_C$$ Unstable", "","", "$$COE$$ Stable", 'interpreter','latex')
            
            else % COE Unstable
                scatter(centers(4,1), centers(4,2), 100, "k", 'LineWidth',2)

                % Separatrix
                % sep_index = 1;
                % separatrix = [];
                % for i=linspace(min_Ecol1,max_Ecol1,1)
                %     for j=linspace(min_Ecol2,max_Ecol2,100000)
                %         tspan = linspace(0,5000,5);
                %         [t,y] = ode45(odefun2, tspan, [i,j]);
                %         if y(end,1)==centers(4,1)
                %             separatrix(sep_index, 1:2) = [i,j]
                %             sep_index = sep_index + 1;
                %         end
                %     end
                % end
                % separatrix
                % semilogy(separatrix, 'Color', colrs(2,:), 'LineWidth',2)
                legend("Direction Field Arrows", "", "","", "", "$$L$$ nullclines", "$$C$$ nullclines", "$$E_0, COE$$ Unstable", "$$E_L, E_C$$ Stable", "", "", "Separatrix", 'interpreter','latex')
            end
         
        % Other legends
        elseif (gammaL > Kc*(1-mu)/Kl && gammaC < Kl/(Kc*(1-mu))) %E_L Stable
            legend("Direction Field Arrows", "", "", "$$L$$ nullclines", "$$C$$ nullclines", "", "", "$$E_0, E_C$$ Unstable", "$$E_L$$ Stable", "",'interpreter','latex')
        elseif (gammaL < Kc*(1-mu)/Kl && gammaC > Kl/(Kc*(1-mu))) %E_C Stable
            legend("Direction Field Arrows", "", "", "$$L$$ nullclines", "$$C$$ nullclines", "", "", "$$E_0, E_L$$ Unstable", "", "$$E_C$$ Stable", "",'interpreter','latex')
        else %All Unstable
            legend("Direction Field Arrows", "", "", "$$L$$ nullclines", "$$C$$ nullclines", "", "", "$$E_0, E_L, E_C$$ Unstable", "", "", "",'interpreter','latex')
        end
        
        
        % Titles, Labels, and File Saving
        if gammaL < 1 && gammaC < 1
            plot_title = "Phase Portrait: $$\gamma_L, \gamma_C < 1, \mu=$$ " + string(mu);
        elseif gammaL > 1 && gammaC > 1
            plot_title = "Phase Portrait: $$\gamma_L, \gamma_C > 1, \mu=$$ " + string(mu);
        elseif gammaL > 1 && gammaC < 1
            plot_title = "Phase Portrait: $$\gamma_L > 1, \gamma_C < 1, \mu=$$ " + string(mu);
        elseif gammaL < 1 && gammaC > 1
            plot_title = "Phase Portrait: $$\gamma_L < 1, \gamma_C > 1, \mu=$$ " + string(mu);
        end
        % title(plot_title,'interpreter','latex')
        xlabel('Liver cells $$L$$', 'Interpreter','latex')
        xlim([0-(max(centers(:,1)))*0.025, max_Ecol1])
        ylabel('Cancer cells $$C$$', 'Interpreter','latex')
        ylim([0-(max(centers(:,2)))*0.025, max_Ecol2])
        file_title = "pp_gL" + string(gammaL) + "_gC" + string(gammaC) + "_mu" + string(mu);
        saveas(gcf,file_title+".fig",'fig')
        saveas(gcf,file_title+".eps",'epsc')
        hold off
    fig_index = fig_index + 1;
    end
end

%% initialize general parameters for figure 3

colrs = 1/255*[0 0 255; 160 32 240; 255 0 0; 255 165 0; 0 255 255; 0 0 0]; % blue purple red orange cyan black
% E_0 <- orange
% E_l <- blue
% E_c <- red
% coe <- purple

l0 = 1/3; % initial healthy population
c0 = 0; % initial tumor population
init = [l0 c0];

lC = 1; % lambda_c
v = linspace(0,1.5,1001*1.5); % nu

k_l = ones(1,length(v));
k_c = 1.5*ones(1,length(v));

%% initialize parameters for figs 3 (a) and (b) 

% a = 1.5>1
gL = 1/1.5; % gamma_l
gC = 3; % gamma_c
a=gL*gC;

l_crit = lC*(1-gL*k_l(1)/k_c(1));
c_crit = lC*(1-k_l(1)/(gC*k_c(1)));

b=1-v/lC; % b=1-\mu, used here for convenience

% equilibrium values
E_0 = 0*ones(1,length(v));
% L_Ec = 0*ones(1,length(v));
C_Ec = b.*k_c;
% L_El = k_l;
C_El = 0*ones(1,length(v));
C_coe = (b.*k_c-gL.*k_l)/(1-a);
L_coe = (k_l-(gC*b).*k_c)/(1-a);

%% generate fig 3 (a) and (b)

% v(556) = l_crit
% v(778) = c_crit
% v(1001) = lC
%'FontSize',14

figure(fig_index)
plot(v(1:end),E_0(1:end),'k--','LineWidth',2) % E_0
hold on
plot(v(1:778),C_El(1:778),'r-','LineWidth',2) % E_C stable
hold on
plot(v(778:end),C_El(778:end),'r--','LineWidth',2) % E_C unstable
hold on
plot(v(1:556),k_l(1:556),'b--','LineWidth',2) % E_L unstable
hold on
plot(v(556:end),k_l(556:end),'b-','LineWidth',2) % E_L stable
hold on
plot(v(556:778),L_coe(556:778),'--','color',colrs(2,:),'LineWidth',2) % L^* COE unstable
ylim([0 1.6])
xlim([0 lC])
xlabel('$\mu$','Interpreter','latex')
ylabel('$L^*$','Interpreter','latex','Rotation',0)
xticks([0 l_crit c_crit lC])
xticklabels({'0','L_{crit}','C_{crit}','1'})
yticks([0 k_l(1) k_c(1)])
yticklabels({'0', 'K_L', 'K_C'})
hAxes.TickLabelInterpreter = 'latex';
vAxes.TickLabelInterpreter = 'latex';
set(gca,'fontsize',16)
legend('E_0','E_C','','','E_L','COE')
fig_index = fig_index + 1;

figure(fig_index)
plot(v(1:end),E_0(1:end),'k--','LineWidth',2) % E_0
hold on
plot(v(1:556),E_0(1:556),'b--','LineWidth',2) % E_L unstable
hold on
plot(v(556:end),E_0(556:end),'b-','LineWidth',2) % E_L stable
hold on
plot(v(1:778),C_Ec(1:778),'r-','LineWidth',2) % E_C stable
hold on
plot(v(778:end),C_Ec(778:end),'r--','LineWidth',2) % E_C unstable
hold on
plot(v(556:778),C_coe(556:778),'--','LineWidth',2,'color',colrs(2,:)) % C^* COE unstable
ylim([0 1.6])
xlim([0 lC])
xlabel('$\mu$','Interpreter','latex')
ylabel('$C^*$','Interpreter','latex','Rotation',0)
xticks([0 l_crit c_crit lC])
xticklabels({'0','L_{crit}','C_{crit}','1'})
yticks([0 k_l(1) k_c(1)])
yticklabels({'0', 'K_L', 'K_C'})
hAxes.TickLabelInterpreter = 'latex';
vAxes.TickLabelInterpreter = 'latex';
set(gca,'fontsize',16)
legend('E_0','','E_L','E_C','','COE')
fig_index = fig_index + 1;

%% initialize parameters for figs 3 (c) and (d) 

% a=0.4<1
gL=1/10.5; % gamma_l
gC=2; % gamma_c
a=gL*gC;

l_crit = lC*(1-gL*k_l(1)/k_c(1));
c_crit = lC*(1-k_l(1)/(gC*k_c(1)));


b=1-v/lC; % b=1-\mu, used here for convenience

% equilibrium values
E_0 = 0*ones(1,length(v));
% L_Ec = 0*ones(1,length(v));
C_Ec = b.*k_c;
% L_El = k_l;
C_El = 0*ones(1,length(v));
C_coe = (b.*k_c-gL.*k_l)/(1-a);
L_coe = (k_l-(gC*b).*k_c)/(1-a);


%% generate fig 3 (c) and (d)

%v(667) = c_crit
%v(937) = l_crit
%v(1001) = lambda_c

figure(fig_index)
plot(v(1:end),E_0(1:end),'k--','LineWidth',2) % unstable E_0
hold on
plot(v(1:667),E_0(1:667),'r-','LineWidth',2) % stable E_C
hold on
plot(v(667:end),E_0(667:end),'r--','LineWidth',2) % unstable E_C
hold on
plot(v(667:937),L_coe(667:937),'color',colrs(2,:),'LineWidth',2) % stable COE
hold on
plot(v(1:937),k_l(1:937),'b--','LineWidth',2) % unstable E_L
hold on
plot(v(937:end),k_l(937:end),'b-','LineWidth',2) % stable E_L
ylim([0 1.6])
xlim([0 lC])
xlabel('$\mu$','Interpreter','latex')
ylabel('$L^*$','Interpreter','latex','Rotation',0)
xticks([0 c_crit l_crit lC])
xticklabels({'0','C_{crit}','L_{crit}','1'})
yticks([0 k_l(1) k_c(1)])
yticklabels({'0', 'K_L', 'K_C'})
hAxes.TickLabelInterpreter = 'latex';
vAxes.TickLabelInterpreter = 'latex';
set(gca,'fontsize',16)
legend('E_0','E_C','','COE','','E_L')
fig_index = fig_index + 1;

figure(fig_index)
plot(v(1:end),E_0(1:end),'k--','LineWidth',2) % unstable E_0
hold on
plot(v(1:667),C_Ec(1:667),'r-','LineWidth',2) % stable E_C
hold on
plot(v(667:end),C_Ec(667:end),'r--','LineWidth',2) % unstable E_C
hold on
plot(v(667:937),C_coe(667:937),'color',colrs(2,:),'LineWidth',2) % stable COE
hold on
plot(v(1:937),E_0(1:937),'b--','LineWidth',2) % unstable E_L
hold on
plot(v(937:end),E_0(937:end),'b-','LineWidth',2) % stable E_L
ylim([0 1.6])
xlim([0 lC])
xlabel('$\mu$','Interpreter','latex')
ylabel('$C^*$','Interpreter','latex','Rotation',0)
xticks([0 c_crit l_crit lC])
xticklabels({'0','C_{crit}','L_{crit}','1'})
yticks([0 k_l(1) k_c(1)])
yticklabels({'0', 'K_L', 'K_C'})
hAxes.TickLabelInterpreter = 'latex';
vAxes.TickLabelInterpreter = 'latex';
set(gca,'fontsize',16)
legend('E_0','E_C','','COE','','E_L')
fig_index = fig_index + 1;

%% initialize and plot raw data

dataChoice = input('Which dataset are you using?\n (0) Liver only (1) Tumor (Anton M1) \n');

% dataChoice == 0 <- Koniaris liver regeneration data (already nondimensional)
% dataChoice == 1 <- Anton M1 mouse tumor growth (exponential)

if dataChoice == 0 %fitting lambda_l only (k_l = 1)
    timespan = 0:2:12; %(days)
    timespan = timespan/7; %(weeks)
    data = 1500*[0.33 0.41 0.54 0.726 0.995 1.105 1]; % mm^3
    error_bars = 1500*[0 0.02 0.02 0.026 0.05 0.105 0.09]; % mm^3
    l0 = data(1);
    c0 = 0;
elseif dataChoice == 1 % fitting lambda_c, gamma_l
    timespan = [0 3 17 21 24 28 31]; % timespan of the experiment (days)
    timespan = timespan/7; % convert to weeks
    data = 1000*[3.3/1000 0.2356114635 0.2520757125 0.5022632175 0.9649201295 2.011033185 3.435962131]; % mm^3
    l0 = 1500; % inital liver population
    c0 = data(1); % inital liver population
end

% plotting data
figure(fig_index)
plot(timespan,data,'r-o','LineWidth',1.5)
if dataChoice == 0
    hold on
    errorbar(timespan,data,error_bars,'.r','LineWidth',1.25);
    xlabel('Days since PH')
    ylabel('Fraction of Original Liver')
else
    xlabel('Weeks post inocualtion')
    ylabel('Tumor (mm^3)')
end
set(gca,'FontSize',14)
fig_index = fig_index + 1;

%% initialize params from data-fitting

% parameter/IC values from data-fitting
lL = 2.16918177423846;
lC = 1.80545101338815;
kL=1500;
kC=7500;
gL=0.639865215044126;
gC=0;
nu = 0;
params = [lL,lC,kL,kC,gL,gC,nu];
init = [1500 3.3 1503.3];
fineMesh = linspace(0,timespan(end),1001); % Time grid

colors = 1/255*[0 0 255; 255 0 0; 160 32 240]; % blue red purple

%solve ode
if dataChoice == 0
    [t,y] = ode23s(@(t,y) odes_after_fitting(t,y,params), fineMesh, init);
elseif dataChoice == 1
    [t,y] = ode23s(@(t,y) dimensional_odes(t,y,params), fineMesh, init);
end

%% generating figs 4 (a) and (b)

figure(fig_index) % plots liver cell population vs real data (if any)
plot(fineMesh, y(:,1),'LineWidth',1.5,'Color',colors(1,:))
xlim([0 fineMesh(end)])
xlabel('Time (weeks)')
ylabel('Relative Liver Volume')
if dataChoice ~= 1 && dataChoice ~= 2 && dataChoice ~= 3 && dataChoice ~= 4
    hold on
    scatter(timespan,data,'filled','r')
    hold on
    errorbar(timespan,data,error_bars,'.r','LineWidth',1.25);
    xlabel('Time (weeks)')
    ylabel('Fraction of Original Liver')
    legend('Simulated', 'Observed',Location='northwest')
    hold off
    ylim(1500*[0 1.25])
end
set(gca,'FontSize',15)
fig_index = fig_index + 1;

figure(fig_index) % plots cancer cell population vs real data (if any)
plot(fineMesh, y(:,2),'LineWidth',1.5,'Color',colors(2,:))
xlim([0 timespan(end)])
xlabel('Time (weeks)')
ylabel('Tumor (mm^3)')
if dataChoice == 1 || dataChoice == 2 || dataChoice == 3 || dataChoice == 4
    hold on
    scatter(timespan,data,'filled','r')
    legend('Simulated', 'Observed',Location='northwest')
    hold off
end
set(gca,'FontSize',15)
fig_index = fig_index + 1;

%% initialize params from parameter_control()
clear params lL lC kL kC gL gC nu l0 c0

% parameter/IC values stored in parameters() function
[lL,lC,kL,kC,gL,gC,nu,~,~,~,~,~,fineMesh,~] = parameter_control();
params = [lL,lC,kL,kC,gL,7.1,nu];
l0=500;
c0=50;
init = [l0 c0 l0+c0];

colors = 1/255*[0 0 255; 255 0 0; 160 32 240]; % blue red purple
mylinestyles = ["-"; "--"; "-.";":"]; % linestyles for figs 5 and 6

L_crit = (1-(gL*kL)/kC)

%% generating figs 5 (a) and (b)

% varying_gammaC = [gL 1.25 5 10];
varying_gammaC = [0.1 gL 1.25 5];

for i = 1:length(varying_gammaC)
    params(6) = varying_gammaC(i);
    [~,y2] = ode23s(@(t,y) odes_after_fitting(t,y,params), fineMesh, init);

    figure(fig_index) % liver
    plot(fineMesh, y2(:,1),'Color',colors(1,:),'LineStyle',mylinestyles(i,:),'LineWidth',2)
    hold on
    if i==length(varying_gammaC)
        legendCellH = cellstr(num2str(varying_gammaC', '$\\gamma_C$=%-.1'));
        legend(legendCellH,'Interpreter','latex','Location','northeast')
    end
    xlim([0 fineMesh(end)])
    xlabel('Time (weeks)')
    ylabel('Volume (mm^3)')
    set(gca,'FontSize',14)
    ylim([-50 1550])

    figure(fig_index+1) % tumor
    plot(fineMesh, y2(:,2),'Color',colors(2,:),'LineStyle',mylinestyles(i,:),'LineWidth',2)
    hold on
    if i==length(varying_gammaC)
        legendCellH = cellstr(num2str(varying_gammaC', '$\\gamma_C$=%-.1d'));
        legend(legendCellH,'Interpreter','latex','Location','east')
    end
    xlim([0 fineMesh(end)])
    xlabel('Time (weeks)')
    ylabel('Volume (mm^3)')
    set(gca,'FontSize',14)
end
fig_index = fig_index + 2;

%% generating figs 6 (a) and (b)

clear t2 y2

params(6) = 1.5; % value of gammaC
varying_mu = [0.1 0.6 0.85 0.975]; % desired values of mu that we are varying
C_crit = (1-Kl/(params(6)*Kc))

% fineMesh2 = linspace(0,4*fineMesh(end),1001);

for i = 1:length(varying_mu)
    params(7) = varying_mu(i)*lC; % converting mu to value of nu for the ode solver
    [t2,y2] = ode23s(@(t,y) odes_after_fitting(t,y,params), fineMesh, init);

    figure(fig_index)
    plot(fineMesh, y2(:,1),'Color',colors(1,:),'LineStyle',mylinestyles(i,:),'LineWidth',2)
    hold on
    if i==length(varying_mu)
        legendCell = cellstr(num2str(varying_mu', '$\\mu$=%-.1d'));
        legend(legendCell,'Interpreter','latex','Location','east')
    end 
    xlim([0 fineMesh(end)])
    ylim([-50 1550])
    xlabel('Time (weeks)')
    ylabel('Volume (mm^3)')
    set(gca,'FontSize',15) 

    figure(fig_index+1)
    plot(fineMesh, y2(:,2),'Color',colors(2,:),'LineStyle',mylinestyles(i,:),'LineWidth',2)
    hold on
    if i==length(varying_mu)
        legendCell = cellstr(num2str(varying_mu', '$\\mu$=%-.1d'));
        legend(legendCell,'Interpreter','latex','Location','east')
    end
    xlim([0 fineMesh(end)])
    ylim([-50 8000])
    xlabel('Time (weeks)')
    ylabel('Volume (mm^3)')
    set(gca,'FontSize',15)
end
fig_index = fig_index + 2;

%% generating figs 6 (c) and (d)

clear t2 y2

params(6) = 5; % value of gammaC
varying_mu = [0.1 0.6 0.85 0.975]; % desired values of mu that we are varying
C_crit = (1-Kl/(params(6)*Kc))

for i = 1:length(varying_mu)
    params(7) = varying_mu(i)*lC; % converting mu to value of nu for the ode solver    
    [t2,y2] = ode23s(@(t,y) odes_after_fitting(t,y,params), fineMesh, init);

    figure(fig_index)
    plot(fineMesh, y2(:,1),'Color',colors(1,:),'LineStyle',mylinestyles(i,:),'LineWidth',2)
    hold on
    if i==length(varying_mu)
        legendCell = cellstr(num2str(varying_mu', '$\\mu$=%-.1d'));
        legend(legendCell,'Interpreter','latex','Location','east')
    end 
    xlim([0 fineMesh(end)])
    ylim([-50 1550])
    xlabel('Time (weeks)')
    ylabel('Volume (mm^3)')
    set(gca,'FontSize',15) 

    figure(fig_index+1)
    plot(fineMesh, y2(:,2),'Color',colors(2,:),'LineStyle',mylinestyles(i,:),'LineWidth',2)
    hold on
    if i==length(varying_mu)
        legendCell = cellstr(num2str(varying_mu', '$\\mu$=%-.1d'));
        legend(legendCell,'Interpreter','latex','Location','east')
    end
    xlim([0 fineMesh(end)])
    ylim([-50 8000])
    xlabel('Time (weeks)')
    ylabel('Volume (mm^3)')
    set(gca,'FontSize',15)
end
fig_index = fig_index + 2;

%% generating figure 7(a)

fig_index=1;

fineMesh2=linspace(0,fineMesh(end)*2,1001);
% fineMesh2=fineMesh;

params(6) = 1.5; % value of gamma_c
params(7) = 0.875*lC; % value of mu
% fineMesh = linspace(0,25,1001);
[~,y2] = ode23s(@(t,y) odes_after_fitting(t,y,params), fineMesh2, init);

figure(fig_index)
plot(fineMesh2, y2(:,1),'Color',colors(1,:),'LineWidth',2)
hold on
plot(fineMesh2, y2(:,2),'Color',colors(2,:),'LineWidth',2)
legend('Healthy','Cancer',Location='east')
xlim([0 fineMesh2(end)])
xlabel('Time (weeks)')
ylabel('Volume (mm^3)')
set(gca,'FontSize',14)
ylim([-50 1550])
fig_index = fig_index + 1;

% %% generating figure 7(b)

params(6) = 2.5; % value of gamma_c
% params(7) = 0.9*lC; % value of mu
% fineMesh = linspace(0,25,1001);
[~,y2] = ode23s(@(t,y) odes_after_fitting(t,y,params), fineMesh2, init);

figure(fig_index)
plot(fineMesh2, y2(:,1),'Color',colors(1,:),'LineWidth',2)
hold on
plot(fineMesh2, y2(:,2),'Color',colors(2,:),'LineWidth',2)
legend('Healthy','Cancer')
xlim([0 fineMesh2(end)])
xlabel('Time (weeks)')
ylabel('Volume (mm^3)')
set(gca,'FontSize',14)
ylim([-50 1550])