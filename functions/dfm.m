function Res = dfm(X,Spec,threshold)
%DFM()    Runs the dynamic factor model
%
%  Syntax:
%    Res = DFM(X,Par)
%
%  Description:
%   DFM() inputs the organized and transformed data X and parameter structure Par.
%   Then, the function outputs dynamic factor model structure Res and data
%   summary statistics (mean and standard deviation).
%
%  Input arguments:
%    X: Transformed input data (i.e. log differenced where needed).
%      Par.Transformation: for mixed frequency aggregation 
%      Par.Frequency: frequency of input data (flexible here)
%      Par.p: Number of lags in transition matrix
%      Par.r: Number of factors (all are treated as global)
%
% Output Arguments:
%
%   Res - structure of model results with the following fields
%       .X_sm  Kalman-smoothed data where missing values are replaced by their expectation
%       .Z  Smoothed states. 
%       .C  Observation matrix. For low frequency data, C gives loadings for aggregated factors
%       .CJ Observation matrix times appropriate helper matrix J for each
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


% DFM input specifications: See documentation for details
if(isfield(Spec, 'p'))
    Par.p = Spec.p;
else
    Par.p = 1;   % Number of lags in autoregressive of factor (same for all factors)
end                        
if(isfield(Spec, 'r'))
    Par.r = Spec.r;
else
    Par.r = 1;   % Number of lags in autoregressive of factor (same for all factors)
end  

fprintf('Estimating the dynamic factor model (DFM) ... \n\n');

[T,N] = size(X);
r = Par.r; %number of factors
p = Par.p;
frq = set_frequencies(Spec.Frequency);
isdiff = is_diff(Spec.Transformation);

if(nargin < 3)
    threshold = 1e-5;  % EM loop threshold (default value)
end

%% Prepare data -----------------------------------------------------------
Mx = mean(X,'omitnan');
Wx = std(X,'omitnan');
xNaN = (X-repmat(Mx,T,1))./repmat(Wx,T,1);  % Standardize series, ie scale()

%% Initial Conditions -----------------------------------------------------

optNaN.method = 2; % Remove leading and closing zeros
optNaN.k = 3;      % Setting for filter(): See remNaN_spline

[A, C, Q, R, V_0] = InitCond(xNaN,r,p,optNaN,frq,isdiff);
Z_0 = zeros(size(A,1),1);

% Initialize EM loop values
previous_loglik = -inf;
num_iter = 0;
LL = -inf;
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

    [C_new, R_new, A_new, Q_new, Z_0, V_0, loglik] = ...  % Applying EM algorithm
        EMstep(Y, A, C, Q, R, Z_0, V_0, p, frq, isdiff);

    C = C_new;
    R = R_new;
    A = A_new;
    Q = Q_new;

    if num_iter > 2  % Checking convergence
        [converged, decrease(num_iter + 1)] = ...
            em_converged(loglik, previous_loglik, threshold, 1);
    end

    if (mod(num_iter,10) == 0) && (num_iter > 0)  % Print updates to command window
        disp(['Now running the ',num2str(num_iter),...
              'th iteration of max ', num2str(max_iter)]);
        disp(['  Loglik','   (% Change)'])
        disp([num2str(loglik),'   (', sprintf('%6.2f',100*((loglik-previous_loglik)/previous_loglik)) '%)'])
    end
    LL = [LL loglik];
    previous_loglik = loglik;
    num_iter =  num_iter + 1;

end

if(num_iter < max_iter)
    disp(['Successful: Convergence at ', num2str(num_iter), ' iterations'])
else
   disp('Stopped because maximum iterations reached')
end

% Final run of the Kalman filter
CJ = [get_CJ(C,frq,isdiff,p), eye(N)];
Zsmooth = runKF(Y, A, CJ, Q, diag(R), Z_0, V_0, r)';

x_sm = Zsmooth(2:end,:) * CJ';  % Get smoothed X

%%  Loading the structure with the results --------------------------------
Res.x_sm = x_sm;
Res.X_sm = 10*repmat(Wx,T,1).*x_sm + repmat(Mx,T,1);  % Unstandardized, smoothed
Res.Z = Zsmooth(2:end,:);
Res.C = C;
Res.CJ = CJ;
Res.R = R;
Res.A = A;
Res.Q = Q;
Res.Mx = Mx;
Res.Wx = Wx;
Res.Z_0 = Z_0;
Res.V_0 = V_0;
Res.r = r;
Res.p = p;

%% Display output
% Table with names and factor loadings

fprintf('\n\n\n');

try
disp('Table 4: Factor Loadings');
disp(array2table(Res.C,...  % Only select lag(0) terms
     'RowNames', strrep(Spec.SeriesName(1:k),' ','_')));
fprintf('\n\n\n');
catch
end

% Table with AR model on factors (factors with AR parameter and variance of residuals)

try
disp('Table 5: Autoregressive Coefficients')
disp(table(A(1:r,1:r*p), ...  
           Q(1:r,1:r), ...
           'VariableNames', {'AR_Coefficient', 'Variance_Residual'}));
          % 'RowNames',      strrep(Spec.BlockNames,' ','_')));
fprintf('\n\n\n');
catch
end

end
