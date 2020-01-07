function bloodflow_analysis
%BLOODFLOW_ANALYSIS   Statistical analysis of cerebral blood flow.
%   BLOODFLOW_ANALYSIS performs time series analysis on cerebral blood
%   flow (CBF) traces. In the example analysis, CBF after common carotid
%   artery occlusion is analyzed in two groups: microglia-depleted and
%   control.
%
%   BLOODFLOW_ANALYSIS performs a Principal Component Analysis on the CBF
%   traces and compares the first 3 PC scores between control and effect
%   group (in the example: microglia-depleted). The data is plotted in the
%   3D space of the first 3 PC scores.
%
%   BLOODFLOW_ANALYSIS performs curve fitting on the traces. The initial
%   fast deflections and the late slow drift components are fitted
%   separately with a biexponential and a 3rd order polynomial function,
%   respectively.
%
%   BLOODFLOW_ANALYSIS was optimized for the example data. 'Data
%   preprocessing' code should be customized for other data formats.
%   If automated curve fitting failed, the CFTOOL function was used to
%   derive the optimal fits (with lowest error).
%
%   Outputs:
%       Figures:
%              * avg data by groups with SE
%              * PCA scores, 3D plot
%              * exp curve fitting for fast changes (first part)
%              * polynomial curve fitting for slow changes (second part)
%
%   BLOODFLOW_ANALYSIS requires Curve Fitting Toolbox.
%
%   See also PCA, FIT, CFTOOL, INPAINT_NANS and SMOOTH.

% Diana Balazsfi and Balazs Hangya
% Institute of Experimental Medicine, Budapest
% balazsfi.diana@koki.mta.hu
% last modified 06.11.2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data
% Example: normalized perfusion unit changes in time:
% Rows contain: animals (merged brain areas-ipsilateral hemisphere)
% Columns: time (1000 data points)
load CBF_example_data.mat  % animals excluded if there is no data from one occl

% Area indeces
occl1_ctrl_inx = [1:9 28:36 55:63];
occl2_ctrl_inx = [10:18 37:45 64:72];
occl3_ctrl_inx = [19:27 46:54 73:81];
occl1_plx_inx = [1:10 31:40 61:70];
occl2_plx_inx = [11:20 41:50 71:80];
occl3_plx_inx = [21:30 51:60 81:90];

% Interpolate nans
ctrl = inpaint_nans(raw_control);
plx = inpaint_nans(raw_plx);

% Number of animals and lenght of time
[numctrl,timectrl] = size(ctrl);
[numplx,timeplx] = size(plx);
NumCtrl = length(occl1_ctrl_inx);
NumPlx = length(occl1_plx_inx);

% Normalize data (centralization)
ctrlmean = mean(ctrl,2);
plxmean = mean(plx,2);
corrctrl = ctrl - repmat(ctrlmean,1,timectrl); % ctrl
corrplx = plx - repmat(plxmean,1,timeplx); % plx

% Get smoothed data (15 points window moving average)
corrctrls = nan(numctrl, timectrl);
for n = 1:numctrl
    corrctrls(n,:) = smooth(corrctrl(n,:),'linear',31);  % ctrl, smoothed
end

corrplxs = nan(numplx, timeplx);
for n = 1:numplx
    corrplxs(n,:) = smooth(corrplx(n,:),'linear',31); % plx, smoothed
end

% Calculate mean values per occlusions
% Not smoothed
poccl1ctrl = corrctrl(occl1_ctrl_inx,:)';  % definition of groups per occlusions
poccl2ctrl = corrctrl(occl2_ctrl_inx,:)';
poccl3ctrl = corrctrl(occl3_ctrl_inx,:)';
poccl1plx = corrplx(occl1_plx_inx,:)';
poccl2plx = corrplx(occl2_plx_inx,:)';
poccl3plx = corrplx(occl3_plx_inx,:)';

mpoccl1ctrl = mean(poccl1ctrl,2);  % calculating mean
mpoccl2ctrl = mean(poccl2ctrl,2);  % ctrl
mpoccl3ctrl = mean(poccl3ctrl,2);

mpoccl1plx = mean(poccl1plx,2);   % plx
mpoccl2plx = mean(poccl2plx,2);
mpoccl3plx = mean(poccl3plx,2);

% Smoothed data
poccl1ctrls = corrctrls(occl1_ctrl_inx,:)';  % definition of groups per occlusions
poccl2ctrls = corrctrls(occl2_ctrl_inx,:)';
poccl3ctrls = corrctrls(occl3_ctrl_inx,:)';
poccl1plxs = corrplxs(occl1_plx_inx,:)';
poccl2plxs = corrplxs(occl2_plx_inx,:)';
poccl3plxs = corrplxs(occl3_plx_inx,:)';

mpoccl1ctrls = mean(poccl1ctrls,2);  % calculating mean
mpoccl2ctrls = mean(poccl2ctrls,2);  % ctrl
mpoccl3ctrls = mean(poccl3ctrls,2);

mpoccl1plxs = mean(poccl1plxs,2);  % plx
mpoccl2plxs = mean(poccl2plxs,2);
mpoccl3plxs = mean(poccl3plxs,2);

% Calculating standard error for non-smoothed data (errorshade in figs)
time = 1:1:size(poccl1ctrl,1);

se_corrctrl1 = std(poccl1ctrl,[],2) / sqrt(size(poccl1ctrl,2)); % occl1, standard error
se_corrplx1 = std(poccl1plx,[],2) / sqrt(size(poccl1plx,2));

se_corrctrl2 = std(poccl2ctrl,[],2) / sqrt(size(poccl2ctrl,2)); % occl2
se_corrplx2 = std(poccl2plx,[],2) / sqrt(size(poccl2plx,2));

se_corrctrl3 = std(poccl3ctrl,[],2) / sqrt(size(poccl3ctrl,2)); % occl3
se_corrplx3 = std(poccl3plx,[],2) / sqrt(size(poccl3plx,2));

% Figure: plot raw data (non-smoothed) with errorshade
figure % occl1
errorshade(time, mpoccl1ctrl, se_corrctrl1, 'LineColor', 'k', 'ShadeColor', 'k', 'LineWidth', 1)
hold on
errorshade(time, mpoccl1plx, se_corrplx1, 'LineColor', 'r', 'ShadeColor', 'r', 'LineWidth', 1);

figure % occl2
errorshade(time, mpoccl2ctrl, se_corrctrl2, 'LineColor', 'k', 'ShadeColor', 'k', 'LineWidth', 1)
hold on
errorshade(time, mpoccl2plx, se_corrplx2, 'LineColor', 'r', 'ShadeColor', 'r', 'LineWidth', 1);

figure % occl3
errorshade(time, mpoccl3ctrl, se_corrctrl3, 'LineColor', 'k', 'ShadeColor', 'k', 'LineWidth', 1)
hold on
errorshade(time, mpoccl3plx, se_corrplx3, 'LineColor', 'r', 'ShadeColor', 'r', 'LineWidth', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRINCIPAL COMPONENT ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PCA on raw data (not smoothed)
[coeff,score,latent,tsquared,explained,mu] = pca([corrctrl; corrplx]);
score_ctrl = score(1:size(corrctrl,1),:);   % PCA scores corresponding to ctrl measurements
score_plx = score(size(corrctrl,1)+1:end,:);   % PCA scores corresponding to plx measurements

% Figure: plot first three PC scores in 3D
figure % occl1
scatter3(score_ctrl(occl1_ctrl_inx,1),score_ctrl(occl1_ctrl_inx,2),score_ctrl(occl1_ctrl_inx,3),...
    60,'o','MarkerEdgeColor','k','MarkerFaceColor',[0 0 0]);
hold on
scatter3(score_plx(occl1_plx_inx,1),score_plx(occl1_plx_inx,2),score_plx(occl1_plx_inx,3),...
    60,'o','MarkerEdgeColor','r','MarkerFaceColor',[1 0 0]);

figure % occl2  
scatter3(score_ctrl(occl2_ctrl_inx,1),score_ctrl(occl2_ctrl_inx,2),score_ctrl(occl2_ctrl_inx,3),...
    60,'o','MarkerEdgeColor','k','MarkerFaceColor',[0 0 0]);
hold on
scatter3(score_plx(occl2_plx_inx,1),score_plx(occl2_plx_inx,2),score_plx(occl2_plx_inx,3),...
    60,'o','MarkerEdgeColor','r','MarkerFaceColor',[1 0 0]);

figure; % occl3    
scatter3(score_ctrl(occl3_ctrl_inx,1),score_ctrl(occl3_ctrl_inx,2),score_ctrl(occl3_ctrl_inx,3),...
    60,'o','MarkerEdgeColor','k','MarkerFaceColor',[0 0 0]);
hold on
scatter3(score_plx(occl3_plx_inx,1),score_plx(occl3_plx_inx,2),score_plx(occl3_plx_inx,3),...
    60,'o','MarkerEdgeColor','r','MarkerFaceColor',[1 0 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYZE LOCAL EXTREMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate min and max on smoothed data (separate exponential and polynom curve))
% Get mean of occlusions (ctrl and plx merged)
peroccl1 = [poccl1ctrls poccl1plxs];
peroccl2 = [poccl2ctrls poccl2plxs];
peroccl3 = [poccl3ctrls poccl3plxs];

mperoccl1 = mean(peroccl1,2);
mperoccl2 = mean(peroccl2,2);
mperoccl3 = mean(peroccl3,2);

% Get reference value (based on average of all data)
refallctrl = [poccl1ctrls poccl2ctrls poccl3ctrls]; % all ctrl
refallplx = [poccl1plxs poccl2plxs poccl3plxs]; % all plx
refall = [refallctrl refallplx]; % all animals
mrefall = mean(refall,2);

% Get min location of grand average
[~, Iallmin] = min(mrefall); 

% Get max location of grand average
[~, Iallmax] = max(mrefall(Iallmin:end,:)); % search maximum value after the minimum value; 
Iallmax = Iallmax + (Iallmin - 1); % get real position/index of value; 

% Get max location by animals in plx group
[~, I] = max(refallplx(Iallmin:end,:));

% Get sd of max values by animals
se_c = std(I,1,2) / sqrt(length(I)); % standard error
mse = Iallmax + round(se_c); % add se to max location to define search window for individual local maxima

% Find min position in groups
[~,I1c] = min(mpoccl1ctrls(1:Iallmax)); % in ctrls by occlusions
[~,I2c] = min(mpoccl2ctrls(1:Iallmax));
[~,I3c] = min(mpoccl3ctrls(1:Iallmax));

[~,I1p] = min(mpoccl1plxs(1:Iallmax));  % in plx by occlusions
[~,I2p] = min(mpoccl2plxs(1:Iallmax));
[~,I3p] = min(mpoccl3plxs(1:Iallmax));

[~,I1o] = min(mperoccl1(1:Iallmax));  % in all animals (ctrl+plx) by occlusions
[~,I2o] = min(mperoccl2(1:Iallmax));
[~,I3o] = min(mperoccl3(1:Iallmax));

% Find max position in groups (from min position to occl avg max pos.+its se)
[~,I1cx] = max(mpoccl1ctrls(I1c:mse)); % in ctrls by occlusions
[~,I2cx] = max(mpoccl2ctrls(I2c:mse)); 
[~,I3cx] = max(mpoccl3ctrls(I3c:mse));

[~,I1px] = max(mpoccl1plxs(I1p:mse)); % in plx by occlusions
[~,I2px] = max(mpoccl2plxs(I2p:mse));
[~,I3px] = max(mpoccl3plxs(I3p:mse));

[~,I1ox] = max(mperoccl1(I1o:mse)); % in all animals (ctrl+plx) by occlusions
[~,I2ox] = max(mperoccl2(I2o:mse));
[~,I3ox] = max(mperoccl3(I3o:mse));

% Calculate real position of the value(added corresponding (per occl) min index)
I1cx = I1cx + (I1c - 1); % in ctrls by occlusions
I2cx = I2cx + (I2c - 1);
I3cx = I3cx + (I3c - 1);

I1px = I1px + (I1p - 1); % in plx by occlusions
I2px = I2px + (I2p - 1);
I3px = I3px + (I3p - 1);

I1ox = I1ox + (I1o - 1); % in all animals (ctrl+plx) by occlusions
I2ox = I2ox + (I2o - 1);
I3ox = I3ox + (I3o - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT CBF CURVES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fitting of exponential curve for the fast changes
% Define time by occlusions
timeexp1 = (1:1:I1ox)'; % occl1
timeexp2 = (1:1:I2ox)'; % occl2
timeexp3 = (1:1:I3ox)'; % occl3

% Define groups per occlusions
eoccl1ctrl = corrctrl(occl1_ctrl_inx,1:I1ox)';
eoccl2ctrl = corrctrl(occl2_ctrl_inx,1:I2ox)';
eoccl3ctrl = corrctrl(occl3_ctrl_inx,1:I3ox)';
eoccl1plx = corrplx(occl1_plx_inx,1:I1ox)';
eoccl2plx = corrplx(occl2_plx_inx,1:I2ox)';
eoccl3plx = corrplx(occl3_plx_inx,1:I3ox)';

% Invert and transform curve
peoccl1ctrl = -eoccl1ctrl;
peoccl2ctrl = -eoccl2ctrl;
peoccl3ctrl = -eoccl3ctrl;

peoccl1plx = -eoccl1plx;
peoccl2plx = -eoccl2plx;
peoccl3plx = -eoccl3plx;

minppeoccl1ctrl = min(peoccl1ctrl,[],1);
minppeoccl2ctrl = min(peoccl2ctrl,[],1);
minppeoccl3ctrl = min(peoccl3ctrl,[],1);

minppeoccl1plx = min(peoccl1plx,[],1);
minppeoccl2plx = min(peoccl2plx,[],1);
minppeoccl3plx = min(peoccl3plx,[],1);

ppeoccl1ctrl = peoccl1ctrl - repmat(minppeoccl1ctrl,size(peoccl1ctrl,1),1);
ppeoccl2ctrl = peoccl2ctrl - repmat(minppeoccl2ctrl,size(peoccl2ctrl,1),1);
ppeoccl3ctrl = peoccl3ctrl - repmat(minppeoccl3ctrl,size(peoccl3ctrl,1),1);

ppeoccl1plx = peoccl1plx - repmat(minppeoccl1plx,size(peoccl1plx,1),1);
ppeoccl2plx = peoccl2plx - repmat(minppeoccl2plx,size(peoccl2plx,1),1);
ppeoccl3plx = peoccl3plx - repmat(minppeoccl3plx,size(peoccl3plx,1),1);

% Fit exponential curve
% Ctrl
[eoccl1c, eoccl2c, eoccl3c] = deal(nan(NumPlx,4));
zeoccl1c = nan(I1ox,NumPlx);
zeoccl2c = nan(I2ox,NumPlx);
zeoccl3c = nan(I3ox,NumPlx);
[eoccl1cg, eoccl2cg, eoccl3cg] = deal(nan(NumPlx,3));
for n = 1:NumCtrl
    [f1c,g] = fit(timeexp1,ppeoccl1ctrl(1:I1ox,n),'exp2'); % occl1
    eoccl1c(n,:) = coeffvalues(f1c);
    zeoccl1c(:,n) = feval(f1c,1:I1ox);
    eoccl1cg(n,1) = g.rsquare;
    eoccl1cg(n,2) = g.adjrsquare;
    eoccl1cg(n,3) = g.rmse;

    [f2c,g] = fit(timeexp2,ppeoccl2ctrl(1:I2ox,n),'exp2'); % occl2
    eoccl2c(n,:) = coeffvalues(f2c);
    zeoccl2c(:,n) = feval(f2c,1:I2ox);
    eoccl2cg(n,1) = g.rsquare;
    eoccl2cg(n,2) = g.adjrsquare;
    eoccl2cg(n,3) = g.rmse;

    [f3c,g] = fit(timeexp3,ppeoccl3ctrl(1:I3ox,n),'exp2'); % occl3
    eoccl3c(n,:) = coeffvalues(f3c);
    zeoccl3c(:,n) = feval(f3c,1:I3ox);
    eoccl3cg(n,1) = g.rsquare;
    eoccl3cg(n,2) = g.adjrsquare;
    eoccl3cg(n,3) = g.rmse;
end

% Plx
[eoccl1p, eoccl2p, eoccl3p] = deal(nan(NumPlx,4));
zeoccl1p = nan(I1ox,NumPlx);
zeoccl2p = nan(I2ox,NumPlx);
zeoccl3p = nan(I3ox,NumPlx);
[eoccl1pg, eoccl2pg, eoccl3pg] = deal(nan(NumPlx,3));
for n = 1:NumPlx
    [f1p,g] = fit(timeexp1,ppeoccl1plx(1:I1ox,n),'exp2'); % occl1
    eoccl1p(n,:) = coeffvalues(f1p);
    zeoccl1p(:,n) = feval(f1p,1:I1ox);
    eoccl1pg(n,1) = g.rsquare;
    eoccl1pg(n,2) = g.adjrsquare;
    eoccl1pg(n,3) = g.rmse;

    [f2p,g] = fit(timeexp2,ppeoccl2plx(1:I2ox,n),'exp2'); % occl2
    eoccl2p(n,:) = coeffvalues(f2p);
    zeoccl2p(:,n) = feval(f2p,1:I2ox);
    eoccl2pg(n,1) = g.rsquare;
    eoccl2pg(n,2) = g.adjrsquare;
    eoccl2pg(n,3) = g.rmse;

    [f3p,g] = fit(timeexp3,ppeoccl3plx(1:I3ox,n),'exp2'); % occl3
    eoccl3p(n,:) = coeffvalues(f3p);
    zeoccl3p(:,n) = feval(f3p,1:I3ox);
    eoccl3pg(n,1) = g.rsquare;
    eoccl3pg(n,2) = g.adjrsquare;
    eoccl3pg(n,3) = g.rmse;
end

% Figure: Calculate avg of the first/fast changed section (per occlusions)
meoccl1ctrl = mean(eoccl1ctrl,2);
meoccl2ctrl = mean(eoccl2ctrl,2);
meoccl3ctrl = mean(eoccl3ctrl,2);

meoccl1plx = mean(eoccl1plx,2);
meoccl2plx = mean(eoccl2plx,2);
meoccl3plx = mean(eoccl3plx,2);

% Shift data
meoccl3ctrl = -meoccl3ctrl;
meoccl3plx = -meoccl3plx;
pmeoccl3ctrl = meoccl3ctrl - min(meoccl3ctrl);
pmeoccl3plx = meoccl3plx - min(meoccl3plx);

% Fit exponential curve
f3 = fit(timeexp3,pmeoccl3ctrl(1:I3ox,1),'exp2');
zpmeoccl3ctrl(:,1) = feval(f3,1:I3ox);
ff3 = fit(timeexp3,pmeoccl3plx(1:I3ox,1),'exp2');
zpmeoccl3plx(:,1) = feval(ff3,1:I3ox);

% Figure: Fitted exponential curve on occlusion 3 data
figure;plot(f3,'k',timeexp3,pmeoccl3ctrl,'.k');
hold on
plot(ff3,'r',timeexp3,pmeoccl3plx,'.r');

% Fitting of polynomial curve for the slower changes
% Define time by occlusions
timepoccl1 = (I1ox:1:1000)'; % occl1
timepoccl2 = (I2ox:1:1000)'; % occl2
timepoccl3 = (I3ox:1:1000)'; % occl3

%Define groups per occlusions
poccl1ctrl = corrctrl(occl1_ctrl_inx,I1ox:1000)';
poccl2ctrl = corrctrl(occl2_ctrl_inx,I2ox:1000)';
poccl3ctrl = corrctrl(occl3_ctrl_inx,I3ox:1000)';
poccl1plx = corrplx(occl1_plx_inx,I1ox:1000)';
poccl2plx = corrplx(occl2_plx_inx,I2ox:1000)';
poccl3plx = corrplx(occl3_plx_inx,I3ox:1000)';

% Calculate avg of the second/slower changed section (per occlusions)
mpoccl1ctrl = (mean(poccl1ctrl,2));
mpoccl2ctrl = (mean(poccl2ctrl,2));
mpoccl3ctrl = (mean(poccl3ctrl,2));

mpoccl1plx = (mean(poccl1plx,2));
mpoccl2plx = (mean(poccl2plx,2));
mpoccl3plx = (mean(poccl3plx,2));

% Fit polynomial curve
% Ctrl
pctrl1 = nan(size(poccl1ctrl,2),4);
pctrl2 = nan(size(poccl2ctrl,2),4);
pctrl3 = nan(size(poccl3ctrl,2),4);
for n = 1:NumCtrl
    fc1 = fit(timepoccl1,poccl1ctrl(:,n),'poly3'); % occl1
    pctrl1(n,:) = coeffvalues(fc1);

    fc2 = fit(timepoccl2,poccl2ctrl(:,n),'poly3'); % occl2
    pctrl2(n,:) = coeffvalues(fc2);
    
    fc3 = fit(timepoccl3,poccl3ctrl(:,n),'poly3'); % occl3
    pctrl3(n,:) = coeffvalues(fc3);
end

% Plx
pplx1 = nan(size(poccl1plx,2),4);
pplx2 = nan(size(poccl2plx,2),4);
pplx3 = nan(size(poccl3plx,2),4);
for n = 1:NumPlx
    fp1 = fit(timepoccl1,poccl1plx(:,n),'poly3'); % occl1
    pplx1(n,:) = coeffvalues(fp1);

    fp2 = fit(timepoccl2,poccl2plx(:,n),'poly3'); % occl2
    pplx2(n,:) = coeffvalues(fp2);
    
    fp3 = fit(timepoccl3,poccl3plx(:,n),'poly3'); % occl3
    pplx3(n,:) = coeffvalues(fp3);
end

% Figure: Fitted polynomial curves on mean of occl 3 data
f3 = fit(timepoccl3,mpoccl3ctrl(:,1),'poly3');
ff3 = fit(timepoccl3,mpoccl3plx(:,1),'poly3');
figure
plot(f3,'k',timepoccl3,mpoccl3ctrl,'.k');
hold on
plot(ff3,'r',timepoccl3,mpoccl3plx,'.r');
end