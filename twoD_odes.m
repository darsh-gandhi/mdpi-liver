function dydt = twoD_odes(t, y, z)
%twoD_odes := liver ode system
%   Inputs: 
%       t <- time (weeks)
%       y <- state vars (mm^3)
%       z <- parameters (depends on parameter)
%           z = [(1) lambda_l, (2) lambda_c, (3) k_l, (4) k_c, (5) gamma_l, (6) gamma_c, (7) mu]
%   Outputs:
%       dydt <- output of ode (mm^3/weeks)

% Initialize derivative vector
dydt = zeros(2,1);

% Differential equations
dydt(1) = z(1)*y(1)*(1-(y(1)+z(4)*y(2))/z(5));
dydt(2) = z(2)*y(2)*(1-z(7))*(1-(z(3)*y(1)+y(2))/(z(6)*(1-z(7))));
% dydt(3) = z(1)*y(1)*(1-(y(1)+z(4)*y(2))/z(5)) + z(2)*y(2)*(1-z(7))*(1-(z(3)*y(1)+y(2))/(z(6)*(1-z(7))));

end
