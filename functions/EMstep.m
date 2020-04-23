
%% PROCEDURES -------------------------------------------------------------

function  [H_new, R_new, A_new, Q_new, loglik] = EMstep(Y, A, H, Q, R, p, frq, isdiff)
%EMstep    Applies EM algorithm for parameter reestimation
%
%
%  Description:
%    EMstep reestimates parameters based on the Estimation Maximization (EM)
%    algorithm. This is a two-step procedure:
%    (1) E-step: the expectation of the log-likelihood is calculated using
%        previous parameter estimates.
%    (2) M-step: Parameters are re-estimated through the maximisation of
%        the log-likelihood (maximize result from (1)).
%
%    See "Maximum likelihood estimation of factor models on data sets with
%    arbitrary pattern of missing data" for details about parameter
%    derivation (Banbura & Modugno, 2010). This procedure is in much the
%    same spirit.
%
%  Input:
%    y:      Series data
%    A:      Transition matrix
%    H:      Observation matrix
%    Q:      Covariance for transition equation residuals
%    R:      Covariance for observation matrix residuals
%    p:      Number of lags in transition equation
%    frq:    frequency of each series in y
%    isdiff: logical (T/F) indicating if series in y is differenced
%    
%
%  Output:
%    H_new: Updated observation matrix
%    R_new: Updated covariance matrix for residuals of observation matrix
%    A_new: Updated transition matrix
%    Q_new: Updated covariance matrix for residuals for transition matrix
%    loglik: Log likelihood
%
% References:
%   "Maximum likelihood estimation of factor models on data sets with
%   arbitrary pattern of missing data" by Banbura & Modugno (2010).
%   Abbreviated as BM2010
%
%

%% Initialize preliminary values

% Store series/model values
[nobs, T] = size(Y);
[~,m] = size(H);
sA = size(A,1); 
sV = sA-nobs;
pp = sV/m;
ar_start = sV+1;

%% Normalize
chlky = chol(Q(1:m,1:m),'lower')/sqrt(10);
scl = kron(eye(pp),eye(m)/chlky); 
Iscl = kron(eye(pp),chlky);
Q(1:m,1:m) = 10*eye(m); %due to normalization
aa = scl*A(1:sV,1:sV)*Iscl;
A(1:m,1:sV) = aa(1:m,1:sV);
H = H*chlky;

%Initial variance following Hamilton 1994
xx = eye(sA^2) - kron(A,A);
vQ = reshape(Q, (sA)^2, 1);
V_0 = xx\vQ;
V_0 = reshape(V_0,sA,sA); 
V_0 = (V_0 + V_0')/2; %reduce rounding error

Z_0 = zeros(sA,1);

%% ESTIMATION STEP: Compute the (expected) sufficient statistics for a single Kalman filter sequence

% Running the Kalman filter and smoother with current parameters
% Note that log-liklihood is NOT re-estimated after the runKF step: This
% effectively gives the previous iteration's log-likelihood
% For more information on output, see runKF
HJ =get_CJ(H, frq, isdiff, p); 
HJ = [HJ, eye(nobs)];

% A = Spec.A;
% HJ = Spec.HJ;
% Q = Spec.Q;
% R = Spec.R;

[Zsmooth, Vsmooth, VVsmooth, loglik] = runKF(Y, A, HJ, Q, diag(R), Z_0, V_0);
% Vsmooth gives the variance of contemporaneous factors
% VVsmooth gives the covariance of factors at one lag for Watson Engle
% adjustments

% Zsmooth(1:sA,:)*Zsmooth(1:sA,:)'/T

%% MAXIMIZATION STEP (TRANSITION EQUATION)
% See (Banbura & Modugno, 2010) for details.

%%% 2A. UPDATE FACTOR PARAMETERS  ----------------------------
% Initialize output
A_new = A;
Q_new = Q;
mp = m*p;
% ESTIMATE FACTOR PORTION OF Q, A
% Note: EZZ, EZZ_BB, EZZ_FB are parts of equations 6 and 8 in BM 2010

% E[f_t*f_t' | Omega_T]
EZZ = Zsmooth(1:mp, 2:end) * Zsmooth(1:mp, 2:end)'...
    + sum(Vsmooth(1:mp, 1:mp, 2:end) ,3); % WE adjustment

EZZ = (EZZ + EZZ')/2;

% E[f_{t-1}*f_{t-1}' | Omega_T]
EZZ_BB = Zsmooth(1:mp, 1:end-1)*Zsmooth(1:mp, 1:end-1)'...
        + sum(Vsmooth(1:mp, 1:mp, 1:end-1), 3); % WE adjustment
    
EZZ_BB = (EZZ_BB + EZZ_BB')/2;

% E[f_t*f_{t-1}' | Omega_T]
EZZ_FB = Zsmooth(1:m, 2:end)*Zsmooth(1:mp, 1:end-1)'...
    + sum(VVsmooth(1:m, 1:mp, :), 3); % WE adjustment

% Equation 6: Estimate VAR(p) for factor
A_new(1:m,1:mp) = EZZ_FB/EZZ_BB; % VAR coeficients

% Equation 8: Covariance matrix of residuals of VAR
q = (EZZ(1:m,1:m) - A_new(1:m,1:mp)* EZZ_FB') / T; %shortcut --- very clever (again)
Q_new(1:m,1:m) =  (q + q')/2;

%%% Update AR(1) components
% Below 3 estimate the idiosyncratic component (for eqns 6, 8 BM 2010)
% Just a bunch of AR(1) regressions

% E[f_t*f_t' | \Omega_T]
EZZ = diag(diag(Zsmooth(ar_start:end, 2:end) * Zsmooth(ar_start:end, 2:end)'))...
    + diag(diag(sum(Vsmooth(ar_start:end, ar_start:end, 2:end), 3)));

% E[f_{t-1}*f_{t-1}' | \Omega_T]
EZZ_BB = diag(diag(Zsmooth(ar_start:end, 1:end-1)* Zsmooth(ar_start:end, 1:end-1)'))...
       + diag(diag(sum(Vsmooth(ar_start:end, ar_start:end, 1:end-1), 3)));

% E[f_t*f_{t-1}' | \Omega_T]
EZZ_FB = diag(diag(Zsmooth(ar_start:end, 2:end)*Zsmooth(ar_start:end, 1:end-1)'))...
       + diag(diag(sum(VVsmooth(ar_start:end, ar_start:end, :), 3)));

A_new(ar_start:end, ar_start:end) = EZZ_FB * diag(1./diag((EZZ_BB)));  % Equation 6
Q_new(ar_start:end, ar_start:end) = (EZZ - A_new(ar_start:end, ar_start:end)*EZZ_FB') / T;           % Equation 8

% Z_0 = Zsmooth(:,1); %zeros(size(Zsmooth,1),1); %update initial condition
% V_0 = squeeze(Vsmooth(:,:,1));


%% 3 MAXIMIZATION STEP (observation equation)

%%% INITIALIZATION AND SETUP ----------------------------------------------

% LOADINGS
H_new = H;

for j = 1:nobs % Loop through observables
    fq = frq(j);
    y = Y(j,:) - Zsmooth(sV+j,2:end);
    y_idx = ~isnan(y);
    y_obs = y(y_idx);
    if fq==1
        Z_obs = Zsmooth(1:m,2:end); %drop pre sample value Z_0
        Z_obs = Z_obs(:,y_idx); %Z_obs where y observed
        V_obs = sum(Vsmooth(1:m,1:m,logical([0,y_idx])),3); %Vsmooth where y observed
    else
        J = helper_mat(fq,isdiff(j),m,sV);
        Z_obs = J*Zsmooth(1:sV,2:end);
        Z_obs = Z_obs(:,y_idx);
        V_obs = J*sum(Vsmooth(1:sV,1:sV,logical([0,y_idx])),3)*J';
    end
    V_ar = sum(Vsmooth(sV+j,sV+j,logical([0,y_idx])),3); %WE adjustment for AR(1) error term
    EZZ = Z_obs*Z_obs' + V_obs;
    EZZ = (EZZ + EZZ')/2;
    H_new(j,:) = (y_obs*Z_obs')/ EZZ;
    R(j) = (sum((y_obs-H_new(j,:)*Z_obs).^2)+ H_new(j,:)*V_obs*H_new(j,:)' + V_ar)/size(y_obs,2);
end
R_new = R;

return
end