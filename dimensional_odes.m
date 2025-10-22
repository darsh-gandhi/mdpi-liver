function dydt = dimensional_odes(t,y,p)
%dimensional_odes := liver ode system without exponential term in dC/dt
%   Inputs: 
%       t <- time (weeks)
%       y <- state vars (mm^3)
%       p <- parameters (depends on parameter)
%           p = [(1) lambda_l, (2) lambda_c, (3) k_l, (4) k_c, (5) gamma_l, (6) gamma_c, (7) nu]
%   Outputs:
%       dydt <- output of ode (mm^3/weeks)

dydt = -1.0*ones(3,1);

p(5)=5*log(2)/(3*p(2));

dydt(1)=p(1)*y(1)*(1-(y(1)+p(6)*y(2))/p(3));        % dL/dt
dydt(2)=p(2)*(1-p(7)/p(2))*y(2)*(1-(p(5)*y(1))/(p(4)*(1-p(7)/p(2))));   % dC/dt
dydt(3)=(p(1)*y(1)*(1-(y(1)+p(6)*y(2))/p(3))) + ((p(2)-p(7)/p(2))*y(2)*(1-(p(5)*y(1)+y(2))/(p(4)*(1-p(7)/p(2))))); % d(L+C)/dt

end