function y = lakemodel(tspan, y0, Q0, Aq,Aj,V,J_b0,delta_in,epsilon,tau,f1,f2,a1,b1)

function dydt = odefun(t,y)
dydt = zeros(2,1);
Q = Q0 *( 1 + Aq*sin(f1*t));
Si_in = a1*Q^b1;

% Si_in = 3.6E-9*Q^0.9;
% Si_in = 0.2;


J_b = J_b0 * (1+Aj*sin(f2*(t-tau)));
dydt(1) = Q*(Si_in-y(1)/V) - J_b;
dydt(2) = (Q*Si_in*(delta_in-y(2)) - J_b*epsilon)/y(1);

% dydt(2) = (Q*Si_in*(delta_in-y(2)) + J_b0 * (1+Aj*sin(2*pi*(t-tau)))*epsilon)/y(1);
% dydt(1) = Q_in*(1+0.5*sin(2*pi*t)) * Si_in - Q_in*(1+0.5*sin(t))*(y(1)/V)- J_b*(1+0.1*sin(2*pi*(t)));
% dydt(2) = (Q_in*(1+0.5*sin(2*pi*t)) * Si_in *(delta_in-y(2))+ J_b*(1+0.1*sin(2*pi*(t)))*epsilon)/y(1);

end

sol = ode15s(@odefun,tspan,y0);
y=sol;
end

