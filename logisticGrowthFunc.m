function y = logisticGrowthFunc(z,timevect,Io)

r = z (1);
p = z (2);
K = z(3);

IC = Io; % initial number of cases reported

[t,C] = ode45(@logisticGrowth,timevect,IC,[],r,p,K);

incidence = [C(1,1);diff(C(:,1))];

y = incidence;
