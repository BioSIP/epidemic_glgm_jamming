% Forecasting the incidence of jamming attacks in Wireless Sensor Networks
% using Epidemic Logistic Growth Model
% Miguel Lopez, PhD  {m.lopez@uma.es - code.io@icloud.com}
% BioSiP Research Group - University of Malaga - Spain 
% Revision. 3.1  -  Date.  2022/14/12

clear;
close all;

 % Scenario 1: Reached nodes by the coordinator when the jammer node is
 %             near to the network's coordinator
 % Scenario 2: Reached nodes by the coordinator when the jammer node is 
 %             in the middle of the topology
 % Scenario 3: Reached nodes by the coordinator when the jammer node is
 %             far to the network's coordinator
 % Type of jamming: Random & Reactive
 % Protocols involved: IEEE 802.15.4, AODV, DSR & MPH
 % Cases of study from 1 to 3 are 50 p/s random jamming
 %                from 4 to 6 are 80 p/s random jamming
 %                from 7 to 9 are reactive jamming attack
 
 escenario = 1; % load the empirical data for the chosen scanario
 jammingType = 1; % load the empirical data for the chosen jamming type
 N = 49; % total nodes
 So = 48; % total susceptibles nodes 

 % --- fitting the Generalized Grouth Model --- %
% ---   from early jamming attack phase    --- %
% name of the data file containing the time series and the attack data reported
dataFileName = sprintf('rawjammingdata%d%d.txt',escenario,jammingType);
data = load(dataFileName);

% temporal resolution of the data in seconds
% define the length of early attack phase
startTime = 1; % by default 1
endTime = 6; % by default 1 minute
early_attack_phase = startTime:endTime;

% select the early growth phase of the attack
data1 = data(early_attack_phase,:); 

% create time vector in seconds
timevect = (data1(:,1));

% initial guesses for parameters
r = 0.1; % initial 0.1
p = 0.6; % initial 0.6

% initial condition
Co = data1(startTime,2);
z(1) = r;
z(2) = p;

% lower bound for r and p
LowerBound = [0.01 0.01];

% upper bound for r and p
UpperBound = [1 1.2];

% fit parameters from data series least-square method
[param,resnorm,residual,~,~,~,jacobian]=lsqcurvefit(@generalizedGrowthFunc,z,timevect,data1(:,2),...
    LowerBound,UpperBound,[],Co);

% estimated parameters 
r_est = param(1);
p_est = param(2);

% initial number of cases
Ic(1) = Co; 

% obtain the best model fit. F contains the cumulative curve and 
% f contains curve of incidence 
[~,F]=ode45(@generalizedGrowth,timevect,Ic,[],r_est,p_est);

% obtaian the curve of incidence from the derivative of F
incidence1 = [F(1,1);diff(F(:,1))];
f = incidence1;

% --- determining parameter uncertainty --- %

% number of bootstrap realizations
bootstraps = 300;

% parameter estimates
Ptrue = [r_est p_est];

% F contains the cumulative number of cases derived from fitting the model
% to the data
% f contains incidence curve
% vector that will store parameter estimates from each realization
paramEst = zeros(bootstraps,2);

% SSE for each realization
SSEs = zeros(bootstraps,1);
curves = [];

% generate each of the realizations
for count=1:bootstraps
    
    % bootstrap realization is generated
    yirData = zeros(length(F),1);
    yirData(1) = F(1);
    
    % poisson dist. error structure
        for t=2:length(F)
            newcases_t = F(t)-F(t-1);
            yirData(t,1) = poissrnd(newcases_t,1,1);
        end

    % store realization of the epidemic curve
    curves = [curves yirData];

    % estimate parameters from bootstrap realization
    % by default employ the trust-region reflective algorithm
    options = [];
    [param,resnorm,residual] = lsqcurvefit(@generalizedGrowthFunc,z,timevect,yirData,...
        LowerBound,UpperBound,options,Co);
        
    % P is the vector with the estimated parameters from the realization
    r_est = param(1);
    p_est = param(2);
    
    % stores estimated parameters for each realization
    paramEst(count,:) = param; 
    
end

% estimate mean and 95% confidence interval from distribution of parameter
param_r = [mean(paramEst(:,1)) plims(paramEst(:,1),0.025) plims(paramEst(:,1),0.975)];
param_p = [mean(paramEst(:,2)) plims(paramEst(:,2),0.025) plims(paramEst(:,2),0.975)];
par_r(jammingType) = param_r(end,1);  % mean of r
par_rl(jammingType) = param_r(end,2); % lower bound of the 95% confidence interval
par_ru(jammingType) = param_r(end,3); % upper bound of the 95% confidence interval
par_p(jammingType) = param_p(end,1);  % mean of p
par_pl(jammingType) = param_p(end,2); % lower bound of the 95% confidence interval
par_pu(jammingType) = param_p(end,3); % upper bound of the 95% confidence interval

% --- generating short-term forecasts from early phase attack --- %

% number of data points ahead
forecastingperiod = 4;  % forecast horizon in time slot
forecastingtime = startTime:endTime+forecastingperiod;

% create time slot vector in seconds
timevect2 = data(forecastingtime,1);

% vector to store forecast curves
curvesforecasts1 = [];
cummulative1 = [];

% generate forecast curves from each bootstrap realization
for count = 1:bootstraps
    r_est = paramEst(count,1);
    p_est = paramEst(count,2);
    [~,x] = ode45(@generalizedGrowth,timevect2,Ic,[],r_est,p_est);
    incidence1=[x(1,1);diff(x(:,1))];
    curvesforecasts1 = [curvesforecasts1 incidence1];
    cummulative1 = [cummulative1 x(:,1)];
end

% save parameter estimations for r and p
filename1 = sprintf('param_ggm_%d%d.mat',escenario,jammingType);
save(filename1, 'param_r','param_p')


