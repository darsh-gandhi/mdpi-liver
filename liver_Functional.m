function jFunctional=liver_Functional(z,data,initialvalues,tspan,w,alpha,dataChoice)     

[~,y] = ode23s(@(t,y) dimensional_odes(t,y,z),tspan,initialvalues);

if dataChoice == 0    % fitting to liver only data
    diff = y(:,1) - reshape(data,size(y(:,1)));
elseif dataChoice == 1 % fitting to cancer only data
    diff = y(:,2) - reshape(data,size(y(:,2)));
else                                       % fitting to liver only data
    diff = y(:,1) - reshape(data,size(y(:,1)));
end

%jFunction is the actual functional we're minimizing
jFunctional = w*0.5*sum(diff.^2) + alpha*0.5*sum(z.^2);

end