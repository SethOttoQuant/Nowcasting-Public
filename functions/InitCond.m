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

function [ A, C, Q, R, V_0] = InitCond(x,r,p,optNaN,frq,is_diff)

OPTS.disp=0;  % Turns off diagnostic information for eigenvalue computation
[xBal,~] = remNaNs_spline(x,optNaN);  % Spline without NaNs

[T,~] = size(xBal);  % Time T series number N
[C, ~] = eigs(cov(xBal), r, 'lm'); % Initial guess for loadings
z = xBal*C; % Initial guess for factors
Z = stack_obs(z,p,true);
B = (z(p+1:T,:)'*Z)/(Z'*Z); %transition matrix
E = z(p+1:T,:)-Z*B';
%Adjusting for differenced low frequency data
lags = frq;
lags(is_diff,:) = arrayfun(@(x)(2*x-1),frq(is_diff,:));
pp = max([lags;p]);
Q = zeros(r*pp,r*pp);
Q(1:r,1:r) = E'*E/(T-p);
E = xBal - z*C';
R = mean(E.^2);
if pp>p
    B = [B, zeros(r,r*(pp-p))];
end
A = comp_form(B);
%Initial variance following Hamilton 1994
xx = eye((r*pp)^2) - kron(A,A);
vQ = reshape(Q, (r*pp)^2, 1);
V_0 = xx\vQ;
V_0 = reshape(V_0,r*pp,r*pp); 
return
end