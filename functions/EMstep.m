
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
[k, T] = size(Y);
[~,m] = size(H);
sA = size(A,1); 
pp = sA/m;

%Initial variance following Hamilton 1994
V_0 = long_run_var(A,Q);

Z_0 = zeros(sA,1);

%% ESTIMATION STEP: Compute the (expected) sufficient statistics for a single Kalman filter sequence

% Running the Kalman filter and smoother with current parameters
% Note that log-liklihood is NOT re-estimated after the runKF step: This
% effectively gives the previous iteration's log-likelihood
% For more information on output, see runKF
HJ =get_CJ(H, frq, isdiff, p); 

% A = Spec.A;
% HJ = Spec.HJ;
% Q = Spec.Q;
% R = Spec.R;

[Zsmooth, Vsmooth, VVsmooth, loglik] = runKF(Y, A, HJ, Q, diag(R), Z_0, V_0);
% Vsmooth gives the variance of contemporaneous factors
% VVsmooth gives the covariance of factors at one lag for Watson Engle
% adjustments

%% Normalize
% This normalization is particularly effective at preventing rounding
% errors, but the specific normalization we use should not be that crucial
if pp>1
    head_z = fliplr(reshape(Zsmooth(:, 1), m, pp));
    head_z = head_z(:,1:pp-1);
    ZZsmooth = [head_z, Zsmooth(1:m,:)];
end
  iscl = ZZsmooth*ZZsmooth'/(T+pp);
  iscl = (iscl+iscl')/2;
  iscl = chol(iscl,  "lower");
  scl = pinv(iscl); 
  Scl = kron(eye(p), scl);
  sscl = kron(eye(pp), scl);
  
  Zsmooth = sscl*Zsmooth;

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
    + Scl*sum(Vsmooth(1:mp, 1:mp, 2:end) ,3)*Scl'; % WE adjustment

EZZ = (EZZ + EZZ')/2;

% E[f_{t-1}*f_{t-1}' | Omega_T]
EZZ_BB = Zsmooth(1:mp, 1:end-1)*Zsmooth(1:mp, 1:end-1)'...
        + Scl*sum(Vsmooth(1:mp, 1:mp, 1:end-1), 3)*Scl'; % WE adjustment
    
EZZ_BB = (EZZ_BB + EZZ_BB')/2;

% E[f_t*f_{t-1}' | Omega_T]
EZZ_FB = Zsmooth(1:m, 2:end)*Zsmooth(1:mp, 1:end-1)'...
    + scl*sum(VVsmooth(1:m, 1:mp, :), 3)*Scl'; % WE adjustment

% Equation 6: Estimate VAR(p) for factor
A_new(1:m,1:mp) = EZZ_FB/EZZ_BB; % VAR coeficients

% Equation 8: Covariance matrix of residuals of VAR
q = (EZZ(1:m,1:m) - A_new(1:m,1:mp)* EZZ_FB') / T; %shortcut --- very clever (again)
Q_new(1:m,1:m) =  (q + q')/2;

%% 3 MAXIMIZATION STEP (observation equation)

%%% INITIALIZATION AND SETUP ----------------------------------------------

% LOADINGS
H_new = H;

for j = 1:k % Loop through observables
    fq = frq(j);
    y = Y(j,:);
    y_idx = ~isnan(y);
    y_obs = y(y_idx);
    if fq==1
        Z_obs = Zsmooth(1:m,2:end); %drop pre sample value Z_0
        Z_obs = Z_obs(:,y_idx); %Z_obs where y observed
        V_obs = sum(Vsmooth(1:m,1:m,logical([0,y_idx])),3); %Vsmooth where y observed
    else
        J = helper_mat(fq,isdiff(j),m,sA);
        Z_obs = J*Zsmooth(1:sA,2:end);
        Z_obs = Z_obs(:,y_idx);
        V_obs = J*sum(Vsmooth(1:sA,1:sA,logical([0,y_idx])),3)*J';
    end
    V_obs = scl*V_obs*scl';
    EZZ = Z_obs*Z_obs' + V_obs;
    EZZ = (EZZ + EZZ')/2;
    H_new(j,:) = (y_obs*Z_obs')/ EZZ;
    R(j) = (sum((y_obs-H_new(j,:)*Z_obs).^2)+ H_new(j,:)*V_obs*H_new(j,:)')/size(y_obs,2);
end
R_new = R;

return
end