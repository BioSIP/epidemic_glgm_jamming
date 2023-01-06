function y = generalizedGrowthFunc(z,timevect,I0)

r = z (1);
p = z (2);

IC = I0; % initial condition of attack nodes

[t,x] = ode45(@generalizedGrowth,timevect,IC,[],r,p);

incidence1 = [x(1,1);diff(x(:,1))];

y = incidence1;
