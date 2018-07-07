clc;
clear all;
format long

%% a)
% Input 30 dow jones stock files
AAPL= 'AAPL.csv';
AXP = 'AXP.csv';
BA = 'BA.csv';
CAT = 'CAT.csv';
CSCO = 'CSCO.csv';
CVX = 'CVX.csv';
DIS = 'DIS.csv';
DWDP = 'DWDP.csv';
GE = 'GE.csv';
GS = 'GS.csv';
HD = 'HD.csv';
IBM = 'IBM.csv';
INTC = 'INTC.csv';
JNJ = 'JNJ.csv';
JPM = 'JPM.csv';
KO = 'KO.csv';
MCD = 'MCD.csv';
MMM = 'MMM.csv';
MRK = 'MRK.csv';
MSFT = 'MSFT.csv';
NKE = 'NKE.csv';
PFE = 'PFE.csv';
PG = 'PG.csv';
TRV = 'TRV.csv';
UNH = 'UNH.csv';
UTX = 'UTX.csv';
V = 'V.csv';
VZ = 'VZ.csv';
WMT = 'WMT.csv';
XOM = 'XOM.csv';


file_name = {AAPL AXP BA CAT CSCO CVX DIS DWDP GE GS HD IBM INTC JNJ JPM KO...
    MCD MMM MRK MSFT NKE PFE PG TRV UNH UTX V VZ WMT XOM};


N = length(file_name);

% Read daily prices
for i = 1:N
% Read daily prices
    if(exist(file_name{1,i},'file'))
      fprintf('\nReading Monthly prices datafile - %s\n', file_name{1,i});
      fid = fopen(file_name{1,i});
      vheader{i} = textscan(fid, '%[^,]%*[^\n]');
      dates{i} = vheader{1}(1:end); % separate d into different cells 
      fclose(fid);
      data_prices{i} = dlmread(file_name{1,i}, ',',[1,5,755,5]);
    else
      error('Daily prices datafile does not exist')
    end
    consolidated_price(:,i) = data_prices{1,i};
end

% 36 months in total, and 35 rates of return
T = length(dates{1,1}{1,1}) - 2;

for i = 1:N
    geo_return = 1;
     for j = 1:T
        returns(j,i) = data_prices{1,i}(j+1) / data_prices{1,i}(j) - 1;
        geo_return = geo_return * (1+returns(j,i));
     end
     arithmetic_average_return(i) = mean(returns(:,i));
     geometic_expected_return(i) = geo_return^(1/T) - 1;
     sta_dev(i) = std(returns(:,i));
end

Q = cov(returns)*(T-1)/T;
mu_i = geometic_expected_return;
R = median(geometic_expected_return);

%% b) MVO -- portfolio 1

weight_wo_1 = zeros(N,1);
weight_1 = zeros(N,1);

c = zeros(1,N); %[0 0 0]
A = - mu_i;
b = - R;
Aeq = ones(1,N); %[1 1 1]
beq = 1;
ub = ones(N,1) * inf; %[inf; inf; inf];
lb_without = zeros(N,1); % [0; 0; 0]% without short selling
lb = -ones(N,1) * inf;

[x_wo, fval_wo, exitflag_wo] = quadprog(2*Q, c, A, b, Aeq, beq, lb_without, ub);
if exitflag_wo == 1
    weight_wo_1(:,1) = x_wo;
    value_wo_1 = fval_wo;
else
    weight_wo_1(:,1) = zeros(N,1);
    value_wo_1 = 0;
end

[x, fval, exitflag] = quadprog(2*Q, c, A, b, Aeq, beq, lb, ub);
if exitflag == 1
    weight_1(:,1) = x;
    value_1= fval;
else
    weight_1(:,1) = zeros(N,1);
    value_1 = 0;
end

mu_wo_1 = mu_i * x_wo;
var_wo_1 = x_wo' * Q * x_wo;
sta_wo_1 = sqrt(var_wo_1);

mu_1 = mu_i * x;
var_1 = x' * Q * x;
sta_1 = sqrt(var_1);
%% c) 

% generate ramdom value
a = 0.95;
b = 1.05;
rdm = (b-a).*rand(30,1) + a;

rdm_mu_i = mu_i .* rdm';

weight_wo_2 = zeros(N,1);

c = zeros(1,N); %[0 0 0]
A = - rdm_mu_i;
b = - R;
Aeq = ones(1,N); %[1 1 1]
beq = 1;
ub = ones(N,1) * inf; %[inf; inf; inf];
lb_without = zeros(N,1); % [0; 0; 0]% without short selling

[x_wo, fval_wo, exitflag_wo] = quadprog(2*Q, c, A, b, Aeq, beq, lb_without, ub);
if exitflag_wo == 1
    weight_wo_2(:,1) = x_wo;
    value_wo_2 = fval_wo;
else
    weight_wo_2(:,1) = zeros(N,1);
    value_wo_2 = 0;
end

[x, fval, exitflag] = quadprog(2*Q, c, A, b, Aeq, beq, lb, ub);
if exitflag == 1
    weight_2(:,1) = x;
    value_2 = fval;
else
    weight_2(:,1) = zeros(N,1);
    value_2 = 0;
end

mu_wo_2 = rdm_mu_i * x_wo;
var_wo_2 = x_wo' * Q * x_wo;
sta_wo_2 = sqrt(var_wo_2);

mu_2 = rdm_mu_i * x;
var_2 = x' * Q * x;
sta_2 = sqrt(var_2);
%% d) 

weight_wo_3 = zeros(N,3);
mu_wo_3 = zeros(3,1);
var_wo_3 = zeros(3,1);
sta_wo_3 = zeros(3,1);

weight_3 = zeros(N,3);
mu_3 = zeros(3,1);
var_3 = zeros(3,1);
sta_3 = zeros(3,1);

for i = 1: 3
    
    a = 0.95;
    b = 1.05;
    rdm = (b-a).*rand(30,1) + a;
    rdm_mu_i = mu_i .* rdm';
    
    c = zeros(1,N); %[0 0 0]
    A = - rdm_mu_i;
    b = - R;
    Aeq = ones(1,N); %[1 1 1]
    beq = 1;
    ub = ones(N,1) * inf; %[inf; inf; inf];
    lb_without = zeros(N,1); % [0; 0; 0]% without short selling
    
    [x_wo, fval_wo, exitflag_wo] = quadprog(2*Q, c, A, b, Aeq, beq, lb_without, ub);
    if exitflag_wo == 1
        weight_wo_3(:,i) = x_wo;
        value_wo_3 (i) = fval_wo;
    else
        weight_wo_3(:,i) = zeros(N,1);
        value_wo_3 (i) = 0;
    end
    
    [x, fval, exitflag] = quadprog(2*Q, c, A, b, Aeq, beq, lb, ub);
    if exitflag == 1
        weight_3(:,i) = x;
        value_3(i)= fval;
    else
        weight_3(:,i) = zeros(N,1);
        value_3(i) = 0;
    end
    
    mu_wo_3(i) = rdm_mu_i * x_wo;
    var_wo_3(i) = x_wo' * Q * x_wo;
    sta_wo_3(i) = sqrt(var_wo_3(i));
    
    mu_3(i) = rdm_mu_i * x;
    var_3(i) = x' * Q * x;
    sta_3(i) = sqrt(var_3(i));
end

% construct the average portfolio
avg_weight_wo = (weight_wo_1 + weight_wo_2 + sum(weight_wo_3,2))/5;
avg_mu_wo = (mu_wo_1 + mu_wo_2 + sum(mu_wo_3))/5;
avg_var_wo = avg_weight_wo' * Q * avg_weight_wo;
avg_sta_wo = sqrt(avg_var_wo);

avg_weight = (weight_1 + weight_2 + sum(weight_3,2))/5;
avg_mu = (mu_1 + mu_2 + sum(mu_3))/5;
avg_var = avg_weight' * Q * avg_weight;
avg_sta = sqrt(avg_var);
%% e) this part has been encoded in the above parts

%% plotting
figure(1);
set(gcf, 'color', 'white');
weight = [weight_wo_1';weight_wo_2';avg_weight_wo'];
c = categorical({'port1','port2','port3'});
bar(c,weight,'stacked')
title('Weight Composition of Three No Shorting Cases');
xlabel('Portfolios');
ylabel('Weight');
ylim([0 1.1])
legend({'AAPL' 'AXP' 'BA' 'CAT' 'CSCO' 'CVX' 'DIS' 'DWDP' 'GE'...
    'GS' 'HD' 'IBM' 'INTC' 'JNJ' 'JPM' 'KO' 'MCD' 'MMM' 'MRK' 'MSFT' ...
    'NKE' 'PFE' 'PG' 'TRV' 'UNH' 'UTX' 'V' 'VZ' 'WMT' 'XOM'});

figure(2);
set(gcf, 'color', 'white');
weight = [weight_1';weight_2';avg_weight'];
c = categorical({'port4','port5','port6'});
bar(c,weight,'stacked')
title('Weight Composition of Three Shorting Cases');
xlabel('Portfolios');
ylabel('Weight');
ylim([0 1.1])
legend({'AAPL' 'AXP' 'BA' 'CAT' 'CSCO' 'CVX' 'DIS' 'DWDP' 'GE'...
    'GS' 'HD' 'IBM' 'INTC' 'JNJ' 'JPM' 'KO' 'MCD' 'MMM' 'MRK' 'MSFT' ...
    'NKE' 'PFE' 'PG' 'TRV' 'UNH' 'UTX' 'V' 'VZ' 'WMT' 'XOM'});