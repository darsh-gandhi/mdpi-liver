clear; close all

%% define + initialize

% Parameter Values
lambdaL = 2.16918177423846;
lambdaC = 1.80545133635960;
Kl=1500;
Kc=7500;
gammaL=0.639865100580;
gammaC=5;
mu = 0;

colors = 1/255*[0 0 255; 255 0 0]; % blue red, vector for colors/shading

% Initialize vector for graphing
n = 5000; % Number of steps
Y = linspace(0, 1, n).';
Z = zeros(n, n);
tspan = linspace(0,5000,1001); % timespan of simulations
row_index = 1;

%% validation from qual. analysis

% Calculate critical values
L_crit = lambdaC*(1-(gammaL*Kl)/Kc);
C_crit = lambdaC*(1-Kl/(gammaC*Kc));

if L_crit < C_crit
    fprintf('expecting unstable COE\n')
    % X = linspace(L_crit, C_crit, n).';
    X = linspace(0, lambdaC, n).';
else
    fprintf('expecting stable COE\n')
    X = linspace(C_crit, L_crit, n).';
end

end_behavior = -1.0*ones(n^2,4);

%% finding end behaviors

% Initialize parameters vector
z = [lambdaL, lambdaC, gammaL, gammaC, Kl, Kc, mu];

% X <- nu
% Y <- cancer percent
% Z <- C(t_f)

% system('caffeinate')

for i = 1:n
    i
    nu = X(i);
    z(7) = nu/lambdaC;

    for j = 1:n
        C_percent = Y(j);
        initialvalues = ((1/3)*Kl).*[(1-C_percent), C_percent];
        [t, y] = ode15s(@(t, y) twoD_odes(t, y, z), tspan, initialvalues);

        Z(j, i) = y(end,2);
        end_behavior(row_index, :) = [nu, C_percent, y(end,2), y(end, 1)];
        row_index = row_index + 1;
    end
end

% system('killall caffeinate')
%% generate fig 8

figure(1)
surf(X, Y, Z)
colormap(colors)
% colorbar
clim([0 1])
% colormap(colors)
hold on
shading interp

red1 = patch(X, Y,'r');
blue1 = patch(X, Y,'b');
legend([red1 blue1], {'Cancer-Win', 'Liver-Win'})

xlabel('$\mu$','Interpreter','latex')
ylabel('Fraction $C$ After PH','Interpreter','latex')
zlabel("Final Amount of Cancer")

xticks([L_crit (L_crit + (C_crit-L_crit)/5) (L_crit + 2*(C_crit-L_crit)/5) (L_crit + 3*(C_crit-L_crit)/5) ...
    (L_crit + 4*(C_crit-L_crit)/5) C_crit])
xticklabels({'L_{crit}','','','','','C_{crit}'})

set(gca, 'FontSize', 13)
view(2)

xlim([L_crit C_crit])