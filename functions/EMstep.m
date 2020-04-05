
%% PROCEDURES -------------------------------------------------------------
% Note: Kalman filter (runKF()) is in the 'functions' folder

function  [C_new, R_new, A_new, Q_new, Z_0, V_0, loglik] = EMstep(Y, A, C, Q, R, Z_0, V_0, p, frq, is_diff)
%EMstep    Applies EM algorithm for parameter reestimation
%
%  Syntax:
%    [C_new, R_new, A_new, Q_new, Z_0, V_0, loglik]
%    = EMstep(y, A, C, Q, R, Z_0, V_0, r, p, R_mat, q, nQ, i_idio, blocks)
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
%    C:      Observation matrix
%    Q:      Covariance for transition equation residuals
%    R:      Covariance for observation matrix residuals
%    Z_0:    Initial values of factors
%    V_0:    Initial value of factor covariance matrix
%    r:      Number of common factors for each block (e.g. vector [1 1 1 1])
%    p:      Number of lags in transition equation
%    R_mat:  Estimation structure for quarterly variables (i.e. "tent")
%    q:      Constraints on loadings
%    nQ:     Number of quarterly series
%    i_idio: Indices for monthly variables
%    blocks: Block structure for each series (i.e. for a series, the structure
%            [1 0 0 1] indicates loadings on the first and fourth factors)
%
%  Output:
%    C_new: Updated observation matrix
%    R_new: Updated covariance matrix for residuals of observation matrix
%    A_new: Updated transition matrix
%    Q_new: Updated covariance matrix for residuals for transition matrix
%    Z_0:   Initial value of state
%    V_0:   Initial value of covariance matrix
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
[~, T] = size(Y);
[k,r] = size(C);
m = size(Z_0,1); 
pp = m/r;

chlky = chol(Q(1:r,1:r),'lower');
scl = kron(eye(pp),eye(r)/chlky); 
Iscl = kron(eye(pp),chlky);
Q(1:r,1:r) = eye(r); %due to normalization
A = scl*A*Iscl;
C = C*chlky;

%% ESTIMATION STEP: Compute the (expected) sufficient statistics for a single Kalman filter sequence

% Running the Kalman filter and smoother with current parameters
% Note that log-liklihood is NOT re-estimated after the runKF step: This
% effectively gives the previous iteration's log-likelihood
% For more information on output, see runKF
CJ = get_CJ(C, frq, is_diff, p);
[Zsmooth, Vsmooth, VVsmooth, loglik] = runKF(Y, A, CJ, Q, diag(R), Z_0, V_0, r);
% Vsmooth gives the variance of contemporaneous factors
% VVsmooth gives the covariance of factors at one lag for Watson Engle
% adjustments

%% MAXIMIZATION STEP (TRANSITION EQUATION)
% See (Banbura & Modugno, 2010) for details.

%%% 2A. UPDATE FACTOR PARAMETERS  ----------------------------

rp = r*p;
% ESTIMATE FACTOR PORTION OF Q, A
% Note: EZZ, EZZ_BB, EZZ_FB are parts of equations 6 and 8 in BM 2010

% E[f_t*f_t' | Omega_T]
EZZ = Zsmooth(1:rp, 2:end) * Zsmooth(1:rp, 2:end)'...
    +sum(Vsmooth(1:rp, 1:rp, 2:end) ,3); % WE adjustment

% E[f_{t-1}*f_{t-1}' | Omega_T]
EZZ_BB = Zsmooth(1:rp, 1:end-1)*Zsmooth(1:rp, 1:end-1)'...
        +sum(Vsmooth(1:rp, 1:rp, 1:end-1), 3); % WE adjustment

% E[f_t*f_{t-1}' | Omega_T]
EZZ_FB = Zsmooth(1:r, 2:end)*Zsmooth(1:rp, 1:end-1)'...
    +sum(VVsmooth(1:r, 1:rp, :), 3); % WE adjustment

% Equation 6: Estimate VAR(p) for factor
A(1:r,1:rp) = EZZ_FB/EZZ_BB; % AR coeficients

% Equation 8: Covariance matrix of residuals of VAR
Q(1:r,1:r) = (EZZ(1:r,1:r) - A(1:r,1:rp)* EZZ_FB') / T; %shortcut --- very clever (again)

Z_0 = Zsmooth(:,1); %zeros(size(Zsmooth,1),1); %update initial condition
V_0 = squeeze(Vsmooth(:,:,1));
% Initialize output
A_new = A;
Q_new = Q;

%% 3 MAXIMIZATION STEP (observation equation)

%%% INITIALIZATION AND SETUP ----------------------------------------------

% LOADINGS
C_new = C;

for j = 1:k % Loop through observables
    fq = frq(j);
    y = Y(j,:);
    y_idx = ~isnan(y);
    y_obs = y(y_idx);
    if fq==1
        Z_obs = Zsmooth(1:r,2:end); %drop pre sample value Z_0
        Z_obs = Z_obs(:,y_idx);
        V_obs = sum(Vsmooth(1:r,1:r,logical([0,y_idx])),3);
    else
        J = helper_mat(fq,is_diff(j),r,m);
        Z_obs = J*Zsmooth(:,2:end);
        Z_obs = Z_obs(:,y_idx);
        V_obs = J*sum(Vsmooth(:,:,logical([0,y_idx])),3)*J';
    end
    EZZ = Z_obs*Z_obs' + V_obs;
    C_new(j,:) = (y_obs*Z_obs')/ EZZ;
    R(j) = ((y_obs-C_new(j,:)*Z_obs)*(y_obs-C_new(j,:)*Z_obs)' + C_new(j,:)*V_obs*C_new(j,:)')/size(y_obs,2);
end
R_new = R;

return
end