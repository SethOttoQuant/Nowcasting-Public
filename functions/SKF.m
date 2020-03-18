%______________________________________________________________________
function S = SKF(Y, A, C, Q, R, Z_0, V_0, r)
%SKF    Applies Kalman filter
%
%  Syntax:
%    S = SKF(Y, A, C, Q, R, Z_0, V_0)
%
%  Description:
%    SKF() applies the Kalman filter

%  Input parameters:
%    Y: k-by-nobs matrix of input data
%    A: m-by-m transition matrix 
%    C: k-by-m observation matrix
%    Q: m-by-m covariance matrix for transition equation residuals (mu_t)
%    R: k-by-k covariance for observation matrix residuals (e_t)
%    Z_0: 1-by-m vector, initial value of state
%    V_0: m-by-m matrix, initial value of state covariance matrix
%    r: number of factors
%
%  Output parameters:
%    S.Zm: m-by-nobs matrix, prior/predicted factor state vector
%          (S.Zm(:,t) = Z_t|t-1)
%    S.ZmU: m-by-(nobs+1) matrix, posterior/updated state vector
%           (S.Zm(t+1) = Z_t|t)
%    S.Vm: m-by-m-by-nobs array, prior/predicted covariance of factor
%          state vector (S.Vm(:,:,t) = V_t|t-1)  
%    S.VmU: m-by-m-by-(nobs+1) array, posterior/updated covariance of
%           factor state vector (S.VmU(:,:,t+1) = V_t|t)
%    S.loglik: scalar, value of likelihood function
%    S.k_t: k-by-m Kalman gain
  
%% INITIALIZE OUTPUT VALUES ---------------------------------------------
  % Output structure & dimensions of state space matrix
  m = size(C,2);
  
  % Outputs time for data matrix. "number of observations"
  [k,nobs]  = size(Y);
  
  % Instantiate output
  S.Zm  = nan(m, nobs);       % Z_t | t-1 (prior)
  S.Vm  = nan(m, m, nobs);    % V_t | t-1 (prior)
  S.ZmU = nan(m, nobs+1);     % Z_t | t (posterior/updated)
  S.VmU = nan(m, m, nobs+1);  % V_t | t (posterior/updated)
  S.loglik = 0;
  S.UD = zeros(m,k,nobs);

%% SET INITIAL VALUES ----------------------------------------------------
  Zu = Z_0;  % Z_0|0 (In below loop, Zu gives Z_t | t)
  Vu = V_0;  % V_0|0 (In below loop, Vu guvse V_t | t)
  
  % Store initial values
  S.ZmU(:,1)    = Zu;
  S.VmU(:,:,1)  = Vu;

%% KALMAN FILTER PROCEDURE ----------------------------------------------
  for t = 1:nobs
      %%% CALCULATING PRIOR DISTIBUTION----------------------------------
      
      % Use transition eqn to create prior estimate for factor
      % i.e. Z = Z_t|t-1
      Z   = A * Zu;
      
      % Prior covariance matrix of Z (i.e. V = V_t|t-1)
      %   Var(Z) = Var(A*Z + u_t) = Var(A*Z) + Var(\epsilon) = 
      %   A*Vu*A' + Q
      V   = A * Vu* A' + Q; 
      V   =  0.5 * (V+V');  % Trick to make symmetric
      
      %%% CALCULATING POSTERIOR DISTRIBUTION ----------------------------
       
      % Removes missing series: These are removed from Y, C, and R
      [Y_t, C_t, R_t, ~] = MissData(Y(:,t), C, R); 

      % Check if y_t contains no data. If so, replace Zu and Vu with prior.
      if isempty(Y_t)
          Zu = Z;
          Vu = V;
      else  
          % Steps for variance and population regression coefficients:
          % Var(c_t*Z_t + e_t) = c_t Var(A) c_t' + Var(u) = c_t*V *c_t' + R
          VC  = V * C_t';  
          iF  = pinv(C_t * VC + R_t); %variance of observables
          
          % Matrix of population regression coefficients (QuantEcon eqn #4)
          VCF = VC*iF;  

          % Gives difference between actual and predicted observation
          % matrix values
          innov  = Y_t - C_t*Z;
          
          % Update estimate of factor values (posterior)
          Zu  = Z  + VCF * innov;
          
          % Update covariance matrix (posterior) for time t
          Vu  = V  - VCF * VC';
          Vu   =  0.5 * (Vu+Vu'); % Approximation trick to make symmetric
          
          % Update log likelihood 
          S.loglik = S.loglik + 0.5*(log(det(iF))  - innov'*iF*innov);
      end
      
      %%% STORE OUTPUT----------------------------------------------------
      
      % Store covariance and observation values for t-1 (priors)
      S.Zm(:,t)   = Z;
      S.Vm(:,:,t) = V;

      % Store covariance and state values for t (posteriors)
      % i.e. Zu = Z_t|t   & Vu = V_t|t
      S.ZmU(:,t+1)    = Zu;
      S.VmU(:,:,t+1)  = Vu;
      % Store Update Contribution for contemporanious factors
      S.UD(:,isfinite(Y(:,t)),t) = VCF .* kron(ones(m,1), innov');
  end 
  % Store Kalman gain k_t
  if isempty(Y_t)
      S.k_t = zeros(m,m);
  else
      S.k_t = VCF * C_t;
  end
  
end