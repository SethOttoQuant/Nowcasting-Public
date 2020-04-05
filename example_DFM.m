%%% Dynamic factor model (DFM) and Nowcast  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is a simplified version of the New York Fed Nowcasting Code on
% [Github](https://github.com/FRBNY-TimeSeriesAnalysis/Nowcasting).
% This version does not include blocks for data types (i.e. global, soft,
% real, labor) but does include AR(1) components for series.
% It also accepts a wider range of frequency mixes,
% accepting inputs y, q, m, bw, w, d. Key differences are highlighted
% below. This file both estimates the model and calculates the
% "News" component of revised and new data between 2016-12-16 and
% 2016-12-23.
%
% Note: data must already be correctly formatted and dated. For example,
% for monthly quarterly data, quarterly observations must be refference
% dated by the last day in the quarter. For monthly-weekly data, monthly
% observations must be refference dated by the last week in the month, and
% so on. 

%% Clear workspace and set paths.
close all; clear; clc;
addpath('functions');


%% User inputs.
vintage = '2016-12-16'; % vintage dataset to use for estimation
country = 'US';         % United States macroeconomic data
sample_start  = datenum('2000-01-01','yyyy-mm-dd'); % estimation sample


%% Load model specification and dataset.
% Load model specification structure `Spec`. Note the Spec file has not
% changed from the original NY Fed example, but some inputs are not used. 
Spec = load_spec('Spec_US_example.xls');
% Parse `Spec`
SeriesID = Spec.SeriesID; SeriesName = Spec.SeriesName; Units = Spec.Units; UnitsTransformed = Spec.UnitsTransformed;
% Load data
datafile = fullfile('data',country,[vintage '.xls']);
[X,Time,Z] = load_data(datafile,Spec,sample_start); 
%X - transformed data
%Time - time in Matlab format
%datestr(Time)
%Z - raw data
summarize(X,Time,Spec,vintage); % summarize data


%% Plot raw and transformed data.
% Industrial Production (INDPRO) <fred.stlouisfed.org/series/INDPRO>
idxSeries = strcmp('INDPRO',SeriesID); t_obs = ~isnan(X(:,idxSeries));
figure('Name',['Data - ' SeriesName{idxSeries}]);

subplot(2,1,1); box on;
plot(Time(t_obs),Z(t_obs,idxSeries)); title('raw observed data');
ylabel(Units{idxSeries}); xlim(Time([1 end])); datetick('x','yyyy','keeplimits');

subplot(2,1,2); box on;
plot(Time(t_obs),X(t_obs,idxSeries)); title('transformed data');
ylabel(UnitsTransformed{idxSeries}); xlim(Time([1 end])); datetick('x','yyyy','keeplimits');
pause(1); % to display plot


%% Run dynamic factor model (DFM) and save estimation output as 'ResDFM'.
threshold = 1e-4; % Set to 1e-5 for more robust estimates

% Orignially model parameters were hard coded. In this version they are
% entered below.
% In this version they can be included in Spec, ie
Spec.p = 3; %number of lags
Spec.r = 3; %number of factors (no blocks here)
% Other importnat elements of Spec:
% Spec.Frequency

% A key difference here is the use of the helper matrix J:
%
% J is a helper matrix that performs the appropriate aggregation following
% Mariano Murasawa (2003). For example, for a differenced quarterly series
% and a two factor model where data is either monthly or quarterly (so that
% the transition equation is quarterly) we have
disp(helper_mat(3,true,2,10))
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

Res = dfm(X,Spec,threshold); %Estimate the model
save('ResDFM','Res','Spec');


plot(Time,Res.Z(:,1:Spec.r));

%% Plot common factor and standardized data.
idxSeries = strcmp('INDPRO',SeriesID);
figure('Name','Common Factor and Standardized Data');
plot(Time,Res.x_sm,':'); hold on;
h = plot(Time,Res.Z(:,1)*Res.C(idxSeries,1),'k','LineWidth',1.5); box on;
xlim(Time([1 end])); datetick('x','yyyy','keeplimits');
legend(h,'common factor'); legend boxoff;
pause(5); % to display plot


%% Plot projection of common factor onto Payroll Employment and GDP.
idxSeries = strcmp('PAYEMS',SeriesID); t_obs = ~isnan(X(:,idxSeries));
figure('Name','Common Factor Projection');

subplot(2,1,1); % projection of common factor onto PAYEMS
CommonFactor = Res.CJ(idxSeries,:)*Res.Z'*Res.Wx(idxSeries)+Res.Mx(idxSeries);
plot(Time,CommonFactor,'k'); hold on;
plot(Time(t_obs),X(t_obs,idxSeries),'b'); box on;
title(SeriesName{idxSeries}); xlim(Time([1 end])); datetick('x','yyyy','keeplimits');
ylabel({Units{idxSeries}, UnitsTransformed{idxSeries}});
legend('common component','data'); legend boxoff;

subplot(2,1,2); % projection of common factor onto GDPC1
idxSeries = strcmp('GDPC1',SeriesID); t_obs = ~isnan(X(:,idxSeries));
CommonFactor = Res.CJ(idxSeries,:)*Res.Z'*Res.Wx(idxSeries)+Res.Mx(idxSeries);
plot(Time,CommonFactor,'k'); hold on;
plot(Time(t_obs), X(t_obs,idxSeries),'b'); box on;
title(SeriesName{idxSeries}); xlim(Time([1 end])); datetick('x','yyyy','keeplimits');
ylabel({Units{idxSeries},UnitsTransformed{idxSeries}});
legend('common component','data'); legend boxoff;

%% Get Nowcast Update (simple)

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



















