function Res = dfm(X,X_pred,m,p,frq,isdiff,threshold)
%DFM()    Runs the dynamic factor model
%
%  Input arguments:
%    X: Transformed input data (i.e. log differenced where needed).
%    m: number of factors
%    p: number of lags
%    isdiff: logical (T/F) indicating if series in corresponding column is
%    differenced (for mixed frequency aggregation)
%
% Output Arguments:
%
%   Res - structure of model results with the following fields
%       .X_sm  Kalman-smoothed data where missing values are replaced by their expectation
%       .Z  Smoothed states. 
%       .H  Observation matrix. For low frequency data, C gives loadings for aggregated factors
%       .HJ Observation matrix times appropriate helper matrix J for each
%       series. This is what is actually used to get fitted values of
%       observables. 
%       .R: Covariance for observation matrix residuals
%       .A: Transition matrix
%       .Q: Covariance for transition equation residuals.
%       .Mx: Series mean
%       .Wx: Series standard deviation
%       .Z_0: Initial value of state
%       .V_0: Initial value of covariance matrix
%       .r: Number of common factors for each block
%       .p: Number of lags in transition equation
%
% References:
%
%   Marta Banbura, Domenico Giannone and Lucrezia Reichlin
%   Nowcasting (2010)
%   Michael P. Clements and David F. Hendry, editors,
%   Oxford Handbook on Economic Forecasting.

%% Store model parameters ------------------------------------------------

fprintf('Estimating the dynamic factor model (DFM) ... \n\n');

[T,N] = size(X);

if(nargin < 5)
    threshold = 1e-4;  % EM loop threshold (default value)
end

%% Prepare data -----------------------------------------------------------
Mx = mean(X,'omitnan');
Wx = std(X,'omitnan');
xNaN = 10*(X-repmat(Mx,T,1))./repmat(Wx,T,1);  % Standardize series, ie scale()

frq = set_frequencies(frq);
%% Initial Conditions -----------------------------------------------------

optNaN.method = 2; % Remove leading and closing zeros
optNaN.k = 3;      % Setting for filter(): See remNaN_spline

[A, H, Q, R] = InitCond(xNaN,m,p,optNaN,frq,isdiff);

% Initialize EM loop values
previous_loglik = -inf;
num_iter = 0;
converged = 0;
max_iter = 5000;

% Y for the estimation is WITH missing data
Y = xNaN'; %transpose for faster column-wise access

%% EM LOOP ----------------------------------------------------------------

%The model can be written as
%y = C*Z + e;
%Z = A*Z(-1) + v
%where y is NxT, Z is (pr)xT, etc

% Remove the leading and ending nans
%optNaN.method = 3;
%y_est = remNaNs_spline(xNaN,optNaN)';

while (num_iter < max_iter) && ~converged % Loop until converges or max iter.

    [H_new, R_new, A_new, Q_new, loglik] = ...  % Applying EM algorithm
        EMstep(Y, A, H, Q, R, p, frq, isdiff);

    H = H_new;
    R = R_new;
    A = A_new;
    Q = Q_new;

    if num_iter > 2  % Checking convergence
        [converged, ~] = ...
            em_converged(loglik, previous_loglik, threshold, 1);
    end

    if (mod(num_iter,10) == 0) && (num_iter > 0)  % Print updates to command window
        disp(['Now running the ',num2str(num_iter),...
              'th iteration of max ', num2str(max_iter)]);
        disp(['  Loglik','   (% Change)'])
        disp([num2str(loglik),'   (', sprintf('%6.2f',100*((loglik-previous_loglik)/previous_loglik)) '%)'])
    end
    %LL = [LL loglik];
    previous_loglik = loglik;
    num_iter =  num_iter + 1;

end

if(num_iter < max_iter)
    disp(['Successful: Convergence at ', num2str(num_iter), ' iterations'])
else
   disp('Stopped because maximum iterations reached')
end

% Final run of the Kalman filter
% Normalize
sA = size(A,1);
sV = (sA - N);
pp = sV/m;
chlky = chol(Q(1:m,1:m),'lower');
scl = kron(eye(pp),eye(m)/chlky); 
Iscl = kron(eye(pp),chlky);
Q(1:m,1:m) = eye(m); %due to normalization
A(1:sV,1:sV) = scl*A(1:sV,1:sV)*Iscl;
H = H*chlky;

%Initial variance following Hamilton 1994
xx = eye(sA^2) - kron(A,A);
vQ = reshape(Q, (sA)^2, 1);
V_0 = xx\vQ;
V_0 = reshape(V_0,sA,sA); 

Z_0 = zeros(sA,1);

% Run filter/smoother
HJ = [get_CJ(H,frq,isdiff,p), eye(N)];
[Zsmooth, ~, ~, LogLik, Update] = runKF(X', A, HJ, Q, diag(R), Z_0, V_0);
Zsmooth = Zsmooth(:, 2:end)'; % Drop pre-sample values 

plot(dates,Zsmooth(:,1:2))


x_sm = Zsmooth * HJ';  % Get smoothed X

%%  Loading the structure with the results --------------------------------
Res.x_sm = x_sm;
Res.X_sm = repmat(Wx,T,1).*x_sm/10 + repmat(Mx,T,1);  % Unstandardized, smoothed
Res.Z = Zsmooth;
Res.H = H;
Res.HJ = HJ;
Res.R = R;
Res.A = A;
Res.Q = Q;
Res.Mx = Mx;
Res.Wx = Wx;
Res.Z_0 = Z_0;
Res.V_0 = V_0;
Res.m = m;
Res.p = p;
Res.ll = real(LogLik);
Res.Update = Update;

%% Display output
% Table with names and factor loadings

fprintf('\n\n\n');

try
disp('Table 4: Factor Loadings');
disp(array2table(Res.H,...  % Only select lag(0) terms
     'RowNames', strrep(Spec.colnames(1:N),' ','_')));
fprintf('\n\n\n');
catch
end

% Table with AR model on factors (factors with AR parameter and variance of residuals)

try
disp('Table 5: Autoregressive Coefficients')
disp(table(A(1:m,1:m*p), ...  
           Q(1:m,1:m), ...
           'VariableNames', {'AR_Coefficient', 'Variance_Residual'}));
fprintf('\n\n\n');
catch
end

end
