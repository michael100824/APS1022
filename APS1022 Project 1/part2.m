clc;
clear all;
format long
%% Part A
% Data used for each asset class
file_name = {'FSITX.csv' 'TPINX.csv' 'VIG.csv' 'IWD.csv' ...
    'PRNHX.csv' 'VISVX.csv' 'FSIIX.csv' 'EEM.csv'};

N = length(file_name);

% Read monthly prices
for i = 1:N
    if(exist(file_name{1,i},'file'))
      fprintf('\nReading Monthly prices datafile - %s\n', file_name{1,i});
      fid = fopen(file_name{1,i});
      vheader{i} = textscan(fid, '%[^,]%*[^\n]');
      dates{i} = vheader{1}(1:end); % separate d into different cells 
      fclose(fid);
      data_prices{i} = dlmread(file_name{1,i}, ',',[1,5,61,5]);
    else
      error('Daily prices datafile does not exist')
    end
    consolidated_price(:,i) = data_prices{1,i};
end

% US 10-years T bill is considered as risk free rate 
% since our data is from 2013 to 2018, we take the average number 
rf=(1.91+2.58+2.43+2.09+1.88+2.86)/6;

% Calculate monthly return excess for each assets
returns = consolidated_price(2:61,:) ./ consolidated_price(1:61-1,:) - 1;
returns = returns - rf/1200;

% Calculate the expected monthly excess returns for each asset
mu_i = mean(returns);


% Calculate the covariance matrix of the excess return for 8 assets
Q = cov(returns);

% Net asset values for all eight assets. From left to right, they are:
% 'FSITX' 'TPINX' 'VIG' 'IWD' 'PRNHX' 'VISVX' 'FSIIX' 'EEM'
% the values are measured in Billion
net_asset=[36.15 38.28 33.99 35.92 22.85 29.34 22.87 42.17];

% Calculate the market weights for eight assets we chose above
w_mkt=(net_asset/sum(net_asset))';

% Calculate the parameter lambda
lambda=((w_mkt'*mu_i')-(rf/1200))/(w_mkt'*Q*w_mkt);

% Implied Equilibrium excess return vector
pi=lambda*Q*w_mkt;

% Market capitalization weighted portfolio expected excess return
r_mkt=w_mkt'*pi;

% Market capitalization weighted portfolio variance
var_mkt = w_mkt' * Q * w_mkt;

% Beta
beta = [1.04,-0.24,0.85,0.96,1.01,0.98,0.91,1.15];

% the scaling factor for voariability of the new expected return
tau=0.025;

% View 1: FSIIX will have an absolute excess monthly return of 0.5%
% View 2: TPINX will outperform FSITX by 5 basis points
% View 3: VIG and PRNHX will outperform IWD and VISVX by 0.6%

p33 = net_asset(3)/(net_asset(3)+net_asset(5));
p34 = - net_asset(4)/(net_asset(4)+net_asset(6));
p35 = net_asset(5)/(net_asset(3)+net_asset(5));
p36 = - net_asset(6)/(net_asset(4)+net_asset(6));


P=[0,0,0,0,0,0,1,0;-1,1,0,0,0,0,0,0;0,0,p33,p34,p35,p36,0,0];

% Set the view as our expectation
q=[0.5/100;0.05/100;0.6/100];

% o1,o2,o3: our confidence level to view 1,2 and 3
o1=tau*P(1,:)*Q*P(1,:)';
o2=tau*P(2,:)*Q*P(2,:)';
o3=tau*P(3,:)*Q*P(3,:)';

% Construct the Omega matrix using o1,o2,o3 
Ome=diag([o1,o2,o3]);

% Variance of the View Portfolios
v_vp1=P(1,:)*Q*P(1,:)';
v_vp2=P(2,:)*Q*P(2,:)';
v_vp3=P(3,:)*Q*P(3,:)';

% Calculate the expected return of black-litterman
exp_r=inv(inv(tau*Q)+P'*inv(Ome)*P)*(inv(tau*Q)*pi+P'*inv(Ome)*q);

% New covariance of the new combined distribution
cov_BL=inv(inv(tau*Q)+P'*inv(Ome)*P)+Q;


% Black-litterman asset weights
w_bl=inv(lambda * cov_BL)*exp_r;
w_bl=w_bl/sum(w_bl);

% Black_litterman portfolio excess return
r_bl=w_bl'*exp_r;

% Black_litterman portfolio excess return variance
var_bl=w_bl'*cov_BL*w_bl;

% Difference from black_litterman mu and pi
mu_diff=exp_r-pi;

% Difference from black_litterman portfolio and market portfolio
w_diff=w_bl-w_mkt;

% Market beta normalizer
beta_norm=beta*w_mkt;

% Black-litterman beta
beta_bl=beta*w_bl/beta_norm;

% Black-litterman Residual Return
rret_bl=r_bl-beta_bl*w_mkt'*exp_r;

% Black-litterman Residual Risk
rrsk_bl=sqrt(var_bl-(beta_bl^2)*var_mkt);

% Active Portfolio Beta
beta_pa=beta_bl-1;

% Black-litterman Active Return
aret_bl=r_bl-w_mkt'*exp_r;

% Black-litterman Active Risk
arsk_bl=sqrt((rrsk_bl^2)+(beta_pa^2)*var_mkt);

% Market portfolio sharpe ratio
sr_mkt=r_mkt/sqrt(var_mkt);

% Black-litterman sharpe ratio
sr_bl=r_bl/sqrt(var_bl);

% Black-litterman information ratio:
ir_bl=aret_bl/arsk_bl;

%%%%%%%%%%%%%%%%%%%% ERC portfolio
n=8;
X=repmat(1.0/n, n, 1);
Sigma=Q;

f = @(X) var(X.*(Sigma*X))*10^13; 

w_erc = fmincon(f,X,[],[],ones(1,length(X)),1);
for i=1:8
    l=Q*w_erc;
    var{i}=w_erc(i)*l(i);
end

% ERC profolio expected excess return
r_erc=mu_i*w_erc;

% Return difference between ERC and black-litterman portfolio
r_ercbl=r_erc-r_bl;

% ERC portfolio variance
var_erc=w_erc'*Q*w_erc;

% Variance difference between ERC and black-litterman portfolio
var_ercbl=var_erc-var_bl;

% Weights difference between ERC and black-litterman portfolio
w_ercbl=w_erc-w_bl;