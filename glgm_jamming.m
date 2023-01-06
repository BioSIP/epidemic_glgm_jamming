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
 
 escenario = 3; % load the empirical data for the chosen scanario
 jammingType = 9; % load the empirical data for the chosen jamming type
 N = 49; % total nodes
 So = 48; % total susceptibles nodes 
 jammerPosition = {' Near',' Middle',' Far'};
 protocol = {' AODV',' DSR',' MPH',' AODV',' DSR',' MPH',' AODV',' DSR',' MPH'};
 jammingAttack = {' 50 p/s random ',' 50 p/s random ',' 50 p/s random ',...
     ' 80 p/s random ',' 80 p/s random ',' 80 p/s random ',' reactive ',...
     ' reactive ',' reactive '};

% --- fitting the Generalized Grouth Model --- %
% ---   from early jamming attack phase    --- %

% name of the data file containing the time series and the attack data reported
dataFileName = sprintf('rawjammingdata%d%d.txt',escenario,jammingType);
data = load(dataFileName);

% parameters obtained from the GGM model
dataFileName = sprintf('param_ggm_%d%d.mat',escenario,jammingType);
load(dataFileName, 'param_p','param_r');

% temporal resolution of the data in seconds
% define the length of early attack phase
startTime = 1; % by default 1
endTime = 6; % by default 1 minute
early_attack_phase = startTime:endTime;

% select the early growth phase of the attack
data1 = data(early_attack_phase,:); 

% create time vector in seconds
timevect = (data1(:,1));

% initial guesses for parameters from GLGM model obtained from GGM model
r = param_r(1);
p = param_p(1);
rLow = param_r(2);
rHigh = param_r(3);
pLow = param_p(2);
pHigh = param_p(3);
k = So*0.5;
kLow = 1;
kHigh = So*0.9;

% initial condition
Co = data1(startTime,2);
z(1) = r;
z(2) = p;
z(3) = k;

% lower bound for r, p and K
LowerBound = [rLow pLow kLow];

% upper bound for r, p and K
UpperBound = [rHigh pHigh kHigh];%1.2; best 1.6 and So*0.9

% fit parameters from data series least-square method
[param,resnorm,Residualssss,~,~,~,jacobian]=lsqcurvefit(@logisticGrowthFunc,z,...
    timevect,data1(:,2),LowerBound,UpperBound,[],Co);

% estimated parameters 
r_est = param(1);
p_est = param(2);
k_est = param(3);

% initial number of cases
Ic(1) = Co; 

% obtain the best model fit. F contains the cumulative curve and 
% f contains curve of incidence 
[~,F]=ode45(@logisticGrowth,timevect,Ic,[],r_est,p_est,k_est);

% obtaian the curve of incidence from the derivative of F
incidence1 = [F(1,1);diff(F(:,1))];
f = incidence1;

% plot the model fit and data
figure(jammingType)
subplot(1,3,[1 2])
plot(timevect,incidence1,'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5);
hold on
plot(data(:,1),data(:,2),'bo','Color',[0 0.4470 0.7410],'LineWidth',1.5);
limitsx = [0 55];
xlim(limitsx);
limitsy = [0 (max(data(:,2))+2)];
ylim(limitsy);
xlabel('Time (s)');
ylabel('Attack incidence')
text(-0.12,1.05,'(a)','Units','normalized','fontsize',15,'fontweight','bold',...
    'FontName','Times New Roman');
legend('Model','Data','Box', 'off','Location','northwest')
set(gca,'FontSize', 12,'FontName','Times New Roman')
set(gcf,'color','white')

% plot the Residuals
subplot(1,3,3)
resid = (incidence1-data1(:,2));
scaledresid = resid./std(resid);
stem(timevect,resid)
limitsx = [0 50];
xlim(limitsx);
xlabel('Time (s)');
ylabel('Residuals')
text(-0.26,1.05,'(b)','Units','normalized','fontsize',15,'fontweight','bold',...
    'FontName','Times New Roman');
set(gca,'FontSize', 12,'FontName','Times New Roman')
set(gcf,'color','white')

x0 = 200;
y0 = 200;
width = 620;
height = 320;
set(gcf,'position',[x0,y0,width,height])
hold off
fig1name = sprintf('figure1_glgm_%d%d.tiff',escenario,jammingType);
saveas(gcf,fig1name)

% --- determining parameter uncertainty --- %

% number of bootstrap realizations
bootstraps = 300;

% parameter estimates
Ptrue = [r_est p_est k_est];

% F contains the cumulative number of cases derived from fitting the model to the data
% f contains incidence curve
% vector that will store parameter estimates from each realization
Phatss = zeros(bootstraps,3);

% SSE for each realization
SSEs = zeros(bootstraps,1);
curves = [];

% generate each of the realizations
for count = 1:bootstraps
    
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
    [param,resnorm,Residualsss] = lsqcurvefit(@logisticGrowthFunc,z,timevect,yirData,...
        LowerBound,UpperBound,options,Co);
        
    % P is the vector with the estimated parameters from the realization
    r_est = param(1);
    p_est = param(2);
    k_est = param(3);
    
    % stores estimated parameters for each realization
    Phatss(count,:) = param; 
    
end

% --- generating short-term forecasts from early phase attack --- %

% number of data points ahead
forecastingperiod = 7;  % forecast horizon to 150 sec

% create time slot vector in seconds for the whole data
% temporal resolution of the data in seconds
timevect2 = data((startTime:early_attack_phase(end)+forecastingperiod),1);

% vector to store forecast curves
curvesforecasts1 = [];
cummulative1 = [];

% generate forecast curves from each bootstrap realization
for count = 1:bootstraps
    r_est = Phatss(count,1);
    p_est = Phatss(count,2);
    k_est = Phatss(count,3);
    [~,x] = ode45(@logisticGrowth,timevect2,Ic,[],r_est,p_est,k_est);
    incidence1=[x(1,1);diff(x(:,1))];
    curvesforecasts1 = [curvesforecasts1 incidence1];
    cummulative1 = [cummulative1 x(:,1)];
end

% forecasted values
meanValue = plims(curvesforecasts1',0.5);
lowerValue = plims(curvesforecasts1',0.025);
upperValue = plims(curvesforecasts1',0.975);

CmeanValue = plims(cummulative1',0.5);
ClowerValue = plims(cummulative1',0.025);
CupperValue = plims(cummulative1',0.975);

% estimate mean and 95% confidence interval from distribution of parameter
param_r = [mean(Phatss(:,1)) plims(Phatss(:,1),0.025) plims(Phatss(:,1),0.975)];
param_p = [mean(Phatss(:,2)) plims(Phatss(:,2),0.025) plims(Phatss(:,2),0.975)];
param_K = [mean(Phatss(:,3)) plims(Phatss(:,3),0.025) plims(Phatss(:,3),0.975)];
par_r(jammingType) = param_r(1);  % mean of r
par_rl(jammingType) = param_r(2); % lower bound of the 95% confidence interval
par_ru(jammingType) = param_r(3); % upper bound of the 95% confidence interval
par_p(jammingType) = param_p(1);  % mean of p
par_pl(jammingType) = param_p(2); % lower bound of the 95% confidence interval
par_pu(jammingType) = param_p(3); % upper bound of the 95% confidence interval
par_k(jammingType) = param_K(1);  % mean of K
par_kl(jammingType) = param_K(2); % lower bound of the 95% confidence interval
par_ku(jammingType) = param_K(3); % upper bound of the 95% confidence interval
cad1 = strcat('{\it r} =',num2str(param_r(1),3),' (95% CI: ',num2str(param_r(2),3),...
    ', ',num2str(param_r(3),3),')');
cad2 = strcat('{\it p} =',num2str(param_p(1),3),' (95% CI: ',num2str(param_p(2),3),...
    ', ',num2str(param_p(3),3),')');
cad3 = strcat('{\it K} =',num2str(param_K(1),2),' (95% CI:',num2str(param_K(2),2),...
    ',',num2str(param_K(3),2),')');

% plot empirical distributions of r, p and K
figure(jammingType+10)
subplot(3,3,1)
histogram(Phatss(:,1),'FaceColor',[0.7 0.7 0.7])
xlabel('{\it r}')
ylabel('Frequency')
title(cad1,'FontWeight','normal','FontSize', 12);
set(gca,'FontSize', 12,'FontName','Times New Roman')
text(-0.25,1.18,'(a)','Units','normalized','fontsize',15,'fontweight','bold',...
    'FontName','Times New Roman');
set(gcf,'color','white')
subplot(3,3,2)
histogram(Phatss(:,2),'FaceColor',[0.7 0.7 0.7])
xlabel('{\it p}')
ylabel('Frequency')
title(cad2,'FontWeight','normal','FontSize', 12);
set(gca,'FontSize', 12,'FontName','Times New Roman')
text(-0.15,1.18,'(b)','Units','normalized','fontsize',15,'fontweight','bold',...
    'FontName','Times New Roman');
set(gcf,'color','white')
subplot(3,3,3)
histogram(Phatss(:,3),'FaceColor',[0.7 0.7 0.7])
xlabel('{\it K}')
ylabel('Frequency')
title(cad3,'FontWeight','normal','FontSize', 12);
set(gca,'FontSize', 12,'FontName','Times New Roman')
text(-0.15,1.18,'(c)','Units','normalized','fontsize',15,'fontweight','bold',...
    'FontName','Times New Roman');
set(gcf,'color','white')

% plot model fit and uncertainty around model fit
subplot(3,3,[4 5 6])
hold on
set(gcf,'color','white')
hold off

% plot first forecasted curves
plot(timevect2,curvesforecasts1,'LineStyle','-','Color',[0.8 0.8 0.8],'Linewidth',1.2)
hold on
axis([0 120 0 max(data(:,2))*2])
plot(timevect2,meanValue,'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5);
plot(timevect2,lowerValue,'LineStyle','-','Color',[0.1 0.1 0.1],'LineWidth',1.2);
plot(timevect2,upperValue,'LineStyle','-','Color',[0.1 0.1 0.1],'LineWidth',1.2);
xlabel('Time (s)');
ylabel('Attack incidence')
set(gca,'FontSize', 12,'FontName','Times New Roman')
text(-0.07,1.1,'(d)','Units','normalized','fontsize',15,'fontweight','bold',...
    'FontName','Times New Roman');

% plot time series data
plot(data(:,1),data(:,2),'o','Color',[0 0.4470 0.7410],'LineWidth',1.2);

% plot vertical line separating calibration versus forecast periods
line2 = [timevect(end) 0;timevect(end) max(data(:,2))*2.5];
plot(line2(:,1),line2(:,2),'LineStyle','--','Color',[0.1 0.1 0.1],'LineWidth',1.2);

% plot cummulative cases and uncertainty around model fit
subplot(3,3,[7 8 9])
hold on
set(gcf,'color','white')
hold off

cummulative2 = data(:,3)';
cummulative2(1) = 1;
timevect3 = data(:,1);
plot(timevect2,cummulative1','LineStyle','-','Color',[0.8 0.8 0.8],'LineWidth',1.2)

hold on
[~,C] = ode45(@logisticGrowth,timevect2,Ic,[],param_r(1),param_p(1),param_K(1));
plot(timevect2,C,'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5);
plot(timevect2,plims(cummulative1',0.025),'LineStyle','-','Color',[0.1 0.1 0.1],'LineWidth',1.2);
plot(timevect2,plims(cummulative1',0.975),'LineStyle','-','Color',[0.1 0.1 0.1],'LineWidth',1.2);
plot(timevect3,cummulative2,'LineStyle','-','Color',[0 0.4470 0.7410],'LineWidth',1.5);

% plot vertical line separating calibration versus forecast periods
line2 = [timevect(end) 0;timevect(end) 50];
plot(line2(:,1),line2(:,2),'LineStyle','--','Color',[0.1 0.1 0.1],'LineWidth',1.2);
axis([0 120 0 50])
xlabel('Time (s)');
ylabel('Cummulative incidence')
set(gca,'FontSize', 12,'FontName','Times New Roman')
text(-0.07,1.25,'(e)','Units','normalized','fontsize',15,'fontweight','bold',...
    'FontName','Times New Roman');
x0 = 200;
y0 = 200;
width = 760;
height = 480;
set(gcf,'position',[x0,y0,width,height])
hold off
fig2name = sprintf('figure2_glgm_%d%d.tiff',escenario,jammingType);
saveas(gcf,fig2name)

% save parameter estimations for r, p, K and Ro
% Ro data
S = data(:,4);% susceptibles
ratio1 = S(1)/S(end);
ratio2 = log(ratio1);
ratio3 = 1 - S(end)/S(1);
RoData = ratio2/ratio3;

% Ro GLGM
ratio1 = S(1)/(S(1) - param_K(1));
ratio2 = log(ratio1);
ratio3 = S(1)/param_K(1);
RoGlgmMean = ratio2*ratio3;

ratio1 = S(1)/(S(1) - param_K(2));
ratio2 = log(ratio1);
ratio3 = S(1)/param_K(2);
RoGlgmLow = ratio2*ratio3;

ratio1 = S(1)/(S(1) - param_K(3));
ratio2 = log(ratio1);
ratio3 = S(1)/param_K(3);
RoGlgmHigh = ratio2*ratio3;

RoGlgm = [RoGlgmMean; RoGlgmLow; RoGlgmHigh];
RoData = [RoData; 0; 0];

% save parameter estimations for r, p, K, Ro
filename1 = sprintf('param_glgm_%d%d.mat',escenario,jammingType);
save(filename1, 'param_p','param_r','param_K','RoGlgm','RoData');

