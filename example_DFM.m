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
dates = table2array(readtable('dates.csv'));
datesM = datenum(dates); %dates in Matlab format

%% Run dynamic factor model (DFM) and save estimation output as 'ResDFM'.
threshold = 1e-5; % Set to 1e-5 for more robust estimates

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

%% Plot Factors 
figure('Name','Common Factors');
h = plot(datesM,Res.Z(:,1:m));
xlim(datesM([1 end])); datetick('x','yyyy','keeplimits');
pause(5); % to display plot

Spec.colnames

%plot(datesM,Res.Z(:,1)*Res.H(idxSeries,1),'k','LineWidth',1.5)


%% Plot projection of common factor onto Payroll Employment and GDP.
idxSeries = strcmp('PAYEMS',SeriesID); t_obs = ~isnan(X(:,idxSeries));
figure('Name','Common Factor Projection');

subplot(2,1,1); % projection of common factor onto PAYEMS
CommonFactor = Res.X_sm(:,idxSeries);
plot(Time,CommonFactor,'k'); hold on;
plot(Time(t_obs),X(t_obs,idxSeries),'b'); box on;
title(SeriesName{idxSeries}); xlim(Time([1 end])); datetick('x','yyyy','keeplimits');
ylabel({Units{idxSeries}, UnitsTransformed{idxSeries}});
legend('common component','data'); legend boxoff;

subplot(2,1,2); % projection of common factor onto GDPC1
idxSeries = strcmp('GDPC1',SeriesID); t_obs = ~isnan(X(:,idxSeries));
CommonFactor = Res.X_sm(:,idxSeries);
plot(Time,CommonFactor,'k'); hold on;
plot(Time(t_obs), X(t_obs,idxSeries),'b'); box on;
title(SeriesName{idxSeries}); xlim(Time([1 end])); datetick('x','yyyy','keeplimits');
ylabel({Units{idxSeries},UnitsTransformed{idxSeries}});
legend('common component','data'); legend boxoff;

%% Get Nowcast Update for series 1 (example)



new_vintage = '2016-12-23'; % vintage dataset to use for estimation
new_datafile = fullfile('data',country,[new_vintage '.xls']);
[X_new,Time_new,Z_new] = load_data(new_datafile,Spec,sample_start); 
%X - transformed data
%Time - time in Matlab format
%Z - raw data

%Make sure old and new data is the same size
[N,k] = size(X_new);
T_old = size(X,1); T = size(X_new,1);
if T > T_old
    X = [X; NaN(T-T_old,N)];
end

%Identify periods where old and new data are different due either to new
%releases or revisions
% is_different = @(x,y) any(x(isfinite(x))~=y(isfinite(x))) || any(~isfinite(x) & isfinite(y));
% update_dates = false(T_new,1);

Xs_old = (X-repmat(Res.Mx,T,1))./repmat(Res.Wx,T,1); %scaled version
Xs_new = (X_new-repmat(Res.Mx,T,1))./repmat(Res.Wx,T,1); %scaled version

% Updates in this version come directly from the Kalman Filter. I've 
% updated SKF.m to include the gain times the innovation, using this
% (element wise) multiplication to give the update each series contributes
% to factors. At the end of the sample predictions plus updates will equal
% factor values exactly, earlier dates may be slightly different due to
% smoothing. 

%Filter old data to get Kalman update
UD_old = SKF(Xs_old', Res.A, Res.CJ, Res.Q, diag(Res.R), Res.Z_0, Res.V_0, Spec.r);  
%Filter new data to get Kalman update
UD_new = SKF(Xs_new', Res.A, Res.CJ, Res.Q, diag(Res.R), Res.Z_0, Res.V_0, Spec.r);

%The news content is the difference in the updates at each period
News = UD_new.UD - UD_old.UD; 

%For GDP in particular
series = 'GDPC1';
i_series = find(strcmp(series,Spec.SeriesID)); %series index
GDP_news = zeros(T,k);
for t=1:T
    GDP_news(t,:) = Res.Wx(i_series)*Res.CJ(i_series,:)*squeeze(News(:,:,t));
end

%Print Results for variables refference dated November 2016
idx_date = '01-Nov-2016';
idx = find(strcmp(cellstr(datestr(Time_new)), idx_date));
is_news = GDP_news(idx,:)';
news_idx = find(is_news ~= 0);

fprintf('\n\n\n');

try
disp('Table 6: News Contribution to GDP');
disp(array2table(is_news(news_idx),'RowNames', strrep(Spec.SeriesName(news_idx),' ','_')));
fprintf('\n\n\n');
catch
end



















