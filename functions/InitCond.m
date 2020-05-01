%--------------------------------------------------------------------------

%InitCond()      Calculates initial conditions for parameter estimation
%
%  Description:
%    Given standardized data and model information, InitCond() creates
%    initial parameter estimates. These are intial inputs in the EM
%    algorithm, which re-estimates these parameters using Kalman filtering
%    techniques.
%
%Inputs:
%  - x:      Standardized data
%  - r:      Number of common factors for each block
%  - p:      Number of lags in transition equation
%  - frq: frequency mix
%
%Output:
%  - A:   Transition matrix
%  - C:   Observation matrix
%  - Q:   Covariance for transition equation residuals
%  - R:   Covariance for observation equation residuals
%  - Z_0: Initial value of state
%  - V_0: Initial value of covariance matrix


function [ A, H, Q, R, V_0] = InitCond(X,m,p,frq,isdiff)

[T,k] = size(X);
xBal = zeros(T,k);
for j=1:k
    xBal(:,j) = spline_fill_centered(X(:,j));
end
[H, ~] = eigs(cov(xBal), m, 'lm'); % Initial guess for loadings
z = xBal*H; % Initial guess for factors
Z = stack_obs(z,p,true);
sV = size(Z,2);
B = (z(p+1:T,:)'*Z)/(Z'*Z + eye(sV)); %transition matrix
E = z(p+1:T,:)-Z*B';
%Adjusting for differenced low frequency data
lags = frq;
lags(isdiff,:) = arrayfun(@(x)(2*x-1),frq(isdiff,:));
pp = max([lags;p]);
sA = m*pp; %size of A matrix
Q = zeros(sA,sA);
Q(1:m,1:m) = E'*E/(T-p);
%Shocks to observations
E = xBal - z*H';
E = E(2:T,:);
R = mean(E.^2) + 1;

if pp>p
    B = [B, zeros(m,m*(pp-p))];
end
A = zeros(sA);
A(1:m*pp, 1:m*pp) = comp_form(B);

%Initial variance following Hamilton 1994
V_0 = long_run_var(A,Q);
return
end