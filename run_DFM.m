%%% Dynamic factor model (DFM) and Nowcast  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is a simplified version of the New York Fed Nowcasting Code on
% [Github](https://github.com/FRBNY-TimeSeriesAnalysis/Nowcasting).
% This version does not include blocks for data types (i.e. global, soft,
% real, labor) but does include AR(1) components for series.
% It also accepts a wider range of frequency mixes,
% accepting inputs y, q, m, bw, w, d. Key differences are highlighted
% below.
%
% Note: data must already be correctly formatted and dated. For example,
% for monthly quarterly data, quarterly observations must be refference
% dated by the last day in the quarter. For monthly-weekly data, monthly
% observations must be refference dated by the last week in the month, and
% so on. 

%% Clear workspace and set paths.
close all; clear; clc;
addpath('functions');
addpath('data')

%% Read data and structure
% Required parameters are:
% Spec.p         number of lags
% Spec.m         number of factors
% Spec.isdiff    logical (T/F) indicating whether input data was differenced
% Optional parameters are:
% Spec.colnames  names of columns in data
% Spec.A         Previously estimated transition matrix
% Spec.Q         Previously estimated covar of shocks to factors
% Spec.H         Previously estimated loadings (observation equation)
% Spec.HJ        Previously estimated loadings in mixed frequency format 
% Spec.R         Previously estimated shocks to observations

fname = 'model.json';
Spec = jsondecode(fileread(fname));

% Read data

% Data to estimate the model excludes outliers. However, predictions are
% made with data that includes outliers; it is then up to the data
% scientist how to treat extreme values

% If the first row does not have any observations Matlab drops it. We need 
% prevent this behavior so that rows correspond to dates.
opts = delimitedTextImportOptions;
opts.DataLines = 1;
opts.VariableTypes = 'double';
X = readmatrix('estimation_data_no_outliers.csv', opts);
X_pred = readmatrix('prediction_data_with_outliers.csv', opts);
Trend = readmatrix('low_frequency_trends.csv', opts);
dates = table2array(readtable('dates.csv'));
datesM = datenum(dates); %dates in Matlab format

%% Run dynamic factor model (DFM) and save estimation output as 'ResDFM'.
threshold = 1e-6; % Set to 1e-6 for more robust estimates

% A key difference from the NY Fed code is the use of the helper matrix J:
%
% J is a helper matrix that performs the appropriate aggregation following
% Mariano Murasawa (2003). For example, for a differenced quarterly series
% and a two factor model where data is either monthly or quarterly (so that
% the transition equation is quarterly) we have
% disp(helper_mat(3,true,2,10))
% This matrix is generated automatically using 
%     - frq: the number of high frequency periods in each series. For
%     monthly-quarterly data, a monthly series will be labled 1 and a
%     quarterly series 3.
%     - isdiff: true/false is the data differenced. 
%     - p: the number of lags
%     - m: the total number of factors in the models. For example, with
%     three (contemporaneous) factors and monthly-quarterly data, where some quarterly series
%     are differenced, the companion form of the model must have at least 15
%     factors.  

m = Spec.m;
p = Spec.p;
isdiff = logical(Spec.isdiff);
frq = Spec.frq;

% Arguments are entered explicity so users can change them easily
Res = dfm(X, X_pred, m, p, frq, isdiff, threshold); %Estimate the model
save('ResDFM','Res','Spec');

%% Compare with results estimated by OttoQuant

% Fitted log likelihoods should be very similar
disp(['C++ fitted log likelihood is ', num2str(Spec.fitted_loglikelihood)])
disp(['Matlab fitted log likelihood is ', num2str(Res.fitted_loglikelihood)])

% However, models with similar likelihoods can lead to different results
% when fitted with different data due to the latency of factors
disp(['C++ forecast log likelihood is ', num2str(Spec.forecast_loglikelihood)])
disp(['Matlab forecast log likelihood is ', num2str(Res.forecast_loglikelihood)])

% The best way to address this issue of model identification is through
% Bayesian estimation by simulation. When models with similar likelihoods
% yeild different results, this uncertainty will be captured in the
% increased vairance of draws for predicted values. 

% We can compare the smoothed, filtered C++ estimates with those from Matlab
T = size(X_pred,1);
xpNaN = 10*(X_pred-repmat(Res.Mx,T,1))./repmat(Res.Wx,T,1); 
V_0 = long_run_var(Spec.A,Spec.Q);
Z_0 = zeros(size(Spec.A,1),1);
[Zsmooth_cpp, ~, ~, LogLik_cpp, ~] = runKF(xpNaN', Spec.A, Spec.HJ, Spec.Q, diag(Spec.R), Z_0, V_0);
Zsmooth_cpp = Zsmooth_cpp(:,2:end)'; % drop pre-sample values

figure('Name','Common Factors');
hold on
plot(datesM,Res.Z(:,1));
plot(datesM,Zsmooth_cpp(:,1));
xlim(datesM([1 end])); datetick('x','yyyy','keeplimits');
hold off
pause(5); % to display plot

%% Plot projection of common factors onto the first series

% Data with which we fitted the model have long run trends removed. We'll
% add these trends back in to plot our fitted estimates

%Using a scatter plot for true observations is confusing, so well use a
%line plot with any missing observations filled via cubic spline.
plot_idx = 1; %Change this to plot a different series
x_plot = spline_fill_centered(X_pred(:,plot_idx)) + Trend(:,plot_idx);

figure('Name', Spec.colnames{plot_idx});
hold on
plot(datesM,Res.X_sm(:,plot_idx) + Trend(:,plot_idx), 'r');
plot(datesM,Res.X_upper(:,plot_idx) + Trend(:,plot_idx), '--r');
plot(datesM,Res.X_lower(:,plot_idx) + Trend(:,plot_idx), '--r');
plot(datesM,x_plot, 'b');
xlim(datesM([1 end])); datetick('x','yyyy','keeplimits');
hold off
pause(5); % to display plot




