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

function [ A, C, Q, R, V_0] = InitCond(x,r,p,optNaN,frq,isdiff)

% OPTS.disp=0;  % Turns off diagnostic information for eigenvalue computation
[xBal,~] = remNaNs_spline(x,optNaN);  % Spline without NaNs

[T,n] = size(xBal);  % Time T series number n
[C, ~] = eigs(cov(xBal), r, 'lm'); % Initial guess for loadings
z = xBal*C; % Initial guess for factors
Z = stack_obs(z,p,true);
sV = size(Z,2);
B = (z(p+1:T,:)'*Z)/(Z'*Z + T*eye(sV)); %transition matrix
E = z(p+1:T,:)-Z*B';
%Adjusting for differenced low frequency data
lags = frq;
lags(isdiff,:) = arrayfun(@(x)(2*x-1),frq(isdiff,:));
pp = max([lags;p]);
sA = r*pp+n; %size of A matrix
Q = zeros(sA,sA);
Q(1:r,1:r) = E'*E/(T-p);
%Shocks to observations
E = xBal - z*C';
a = zeros(1,n);
for j = 1:n
    a(j) = E(2:T,j)'*E(1:(T-1),j)/(E(1:T-1,j)'*E(1:(T-1),j));
end
E = E(2:T,:) - repmat(a,T-1,1).*E(1:T-1,:);
R = mean(E.^2);
Q(r*pp+1:sA, r*pp+1:sA) = diag(R);
if pp>p
    B = [B, zeros(r,r*(pp-p))];
end
A = zeros(sA);
A(1:r*pp, 1:r*pp) = comp_form(B);
A(r*pp+1:sA, r*pp+1:sA) = diag(a);
%Initial variance following Hamilton 1994
xx = eye(sA^2) - kron(A,A);
vQ = reshape(Q, (sA)^2, 1);
V_0 = xx\vQ;
V_0 = reshape(V_0,sA,sA); 
return
end