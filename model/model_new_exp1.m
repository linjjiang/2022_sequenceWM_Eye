%% Mixture modeling for experiment 1
% By Linjing Jiang
% 05/02/2023
% Perform mixture modeling for experiment 1, using MemFit2D

%% Load directories
clear
close all
clc

% output directory
output = './analysis/exp1';

% add scripts to the path
addpath(genpath('./scripts/Cloned-Repo/'))

% load stats 
load(fullfile(output,'stats_exp1_longdelay.mat')) % stats1
load(fullfile(output,'stats_exp1_shortdelay.mat')) % stats2


%% Standard model: the v3_short, long delay condition - serial position by cue
% only guessing and memory imprecision

%stats1.model_standard.mux = nan(length(stats1.subs_id),2,4);
%stats1.model_standard.muy = nan(length(stats1.subs_id),2,4);
%stats1.model_standard.mu = nan(length(stats1.subs_id),2,4);
stats1.model_standard.sd = nan(length(stats1.subs_id),2,4);
%stats1.model_standard.beta = nan(length(stats1.subs_id),2,4);
stats1.model_standard.gamma = nan(length(stats1.subs_id),2,4);
stats1.model_standard.alpha = nan(length(stats1.subs_id),2,4);
stats1.model_standard.llh = nan(length(stats1.subs_id),2,4);
stats1.model_standard.model = cell(length(stats1.subs_id),2,4);
stats1.model_standard.data = cell(length(stats1.subs_id),2,4);
stats1.model_standard.AIC = nan(length(stats1.subs_id),2,4);
stats1.model_standard.AICc = nan(length(stats1.subs_id),2,4);
stats1.model_standard.BIC = nan(length(stats1.subs_id),2,4);

for ss = 1:length(stats1.subs_id) % for each participant
    % Let's model quadrant and order cue differently
    for cc = 1:2 % 1, quadrant cue (0); 2, order cue (1)
        for pp = 1:4 % 4 serial positions

% row (trial) index of the specific type of cue
idx = find(stats1.cue(ss,:) == cc-1 & ...
    stats1.kept_trial(ss,:) == 1 & ...
    stats1.ord(ss,:) == pp);

if isempty(idx) | length(idx) == 1 % if only 1 trial or if there is no trial
    continue; % we do not do model
end

% what we need:
% X and Y of primary saccade endpoints
% X and Y of targets
% X and Y of non-targets

% Input:
% testdata.targets: 2 x N matrix, where N is the number of trials
    % two rows: X1, Y1

tarx = stats1.xtar_pix(ss,idx); % 1 x N trials
tary = stats1.ytar_pix(ss,idx); % 1 x N trials
testdata.targets = [tarx;tary];

% testdata.distractors: 6 x N matrix, if we have 3 distractors (3 x 2 = 6)
    % six rows: X1, Y1, X2, Y2, X3, Y3
    % distance between the non-target and target, e.g., distractor minus
    % target (not vice versa!!!!!)

% distractor location
ntarx = squeeze(stats1.xntar_pix(ss,idx,:))'; % 3 nontargets x N trials
ntary = squeeze(stats1.yntar_pix(ss,idx,:))';
% distractor location minus target location
ntarx = ntarx - tarx; ntary = ntary - tary;
% concatenate distractor locations
if length(idx) > 1 % if more than 1 trial
testdata.distractors = [ntarx(1,:);ntary(1,:);ntarx(2,:);ntary(2,:);ntarx(3,:);ntary(3,:)];
else
    testdata.distractors = [ntarx(1);ntary(1);ntarx(2);ntary(2);ntarx(3);ntary(3)];
end
% testdata.errors: 2 x N matrix
    % two rows: X1, Y1
    % distance between the response and target, e.g., response minus target
    % (not vice versa!!!!!)

    % get endpoint fixation position across trials
respx = stats1.p_xend_pix(ss,idx); % 1 x N trials
respy = stats1.p_yend_pix(ss,idx);
% response minus target location
respx = respx - tarx; respy = respy - tary;
% concatenate responses
testdata.errors = [respx;respy];

% testdata.dimensions = [x_resolution y_resolution]; % in pixel
testdata.dimensions = [1280 1024]; % adjust this if needed

% Initialize the model
%model = WithBias2D(SwapModel2D());
%model = SwapModel2D();
model = StandardMixtureModel2D();

SD = mean(std(testdata.errors,0,2));
MU = mean(testdata.errors,2);
MUX= MU(1);
MUY = MU(2);

model_init = AddInit(model,SD,MUX,MUY); % add initial parameters
%model_init = model;

% fitting results
% Output:
% fit: {'g', 'sd'}
% llh: log likelihood value
[fit,llh] = MLE(testdata, model);
% stats.model.mu = [fit(4) nan;fit(5) nan]; % 2 dimenions (x and y) x 2 types of items (target and nontarget)
% stats.model.Sigma(:,:,1) = [fit(3) 0;0 fit(3)]; % sigma for target
% stats.model.Sigma(:,:,2) = stats.model.Sigma(:,:,1); % sigma for non-target - they are the same
% stats.model.w = [1-fit(1)-fit(2) fit(2) fit(1)]; % alpha (target prob), beta (swap rate), gamma (guess rate)
% stats.model.range = testdata.dimensions;
% stats.model.llh = llh;
% stats.model.MemModel = model;
%stats1.model_standard.mux(ss,cc,pp) = fit(4);
%stats1.model_standard.muy(ss,cc,pp)= fit(5);
%stats1.model_standard.mu(ss,cc,pp) = sqrt(fit(4)^2 + fit(5)^2);
stats1.model_standard.sd(ss,cc,pp) = fit(2); % standard deviation
stats1.model_standard.alpha(ss,cc,pp) = 1-fit(1);
%stats1.model_standard.beta(ss,cc,pp) = fit(2);
stats1.model_standard.gamma(ss,cc,pp) = fit(1);
stats1.model_standard.llh(ss,cc,pp) = llh;
stats1.model_standard.model{ss,cc,pp} = model;
stats1.model_standard.data{ss,cc,pp} = testdata;
% Calculate goodness of fit
k = length(model.upperbound); % how many parameters
dataLen = size(respx,2); % how many trials/data points
stats1.model_standard.AIC(ss,cc,pp) = -2*llh + 2 * k; % 2 parameters
stats1.model_standard.AICc(ss,cc,pp) = -2*llh + 2 * k * (dataLen / (dataLen - k - 1)); % 2 parameters
stats1.model_standard.BIC(ss,cc,pp) = -2*llh + (log(dataLen)+log(2*pi))*k; % 2 parameters
        end
    end
end

%% Standard mixture model: the v3_short_DELAY, short delay condition - serial position by cue
% only guessing and memory imprecision

%stats2.model_standard.mux = nan(length(stats2.subs_id),2,4);
%stats2.model_standard.muy = nan(length(stats2.subs_id),2,4);
%stats2.model_standard.mu = nan(length(stats2.subs_id),2,4);
stats2.model_standard.sd = nan(length(stats2.subs_id),2,4);
%stats2.model_standard.beta = nan(length(stats2.subs_id),2,4);
stats2.model_standard.gamma = nan(length(stats2.subs_id),2,4);
stats2.model_standard.alpha = nan(length(stats2.subs_id),2,4);
stats2.model_standard.llh = nan(length(stats2.subs_id),2,4);
stats2.model_standard.model = cell(length(stats2.subs_id),2,4);
stats2.model_standard.data = cell(length(stats2.subs_id),2,4);
stats2.model_standard.AIC = nan(length(stats2.subs_id),2,4);
stats2.model_standard.AICc = nan(length(stats2.subs_id),2,4);
stats2.model_standard.BIC = nan(length(stats2.subs_id),2,4);

for ss = 1:length(stats2.subs_id) % for each participant
    % Let's model quadrant and order cue differently
    for cc = 1:2 % 1, quadrant cue (0); 2, order cue (1)
        for pp = 1:4 % 4 serial positions

% row (trial) index of the specific type of cue
idx = find(stats2.cue(ss,:) == cc-1 & ...
    stats2.kept_trial(ss,:) == 1 & ...
    stats2.ord(ss,:) == pp);

if isempty(idx) | length(idx) == 1 % if only 1 trial or if there is no trial
    continue; % we do not do model
end


% what we need:
% X and Y of primary saccade endpoints
% X and Y of targets
% X and Y of non-targets

% Input:
% testdata.targets: 2 x N matrix, where N is the number of trials
    % two rows: X1, Y1

tarx = stats2.xtar_pix(ss,idx); % 1 x N trials
tary = stats2.ytar_pix(ss,idx); % 1 x N trials
testdata.targets = [tarx;tary];

% testdata.distractors: 6 x N matrix, if we have 3 distractors (3 x 2 = 6)
    % six rows: X1, Y1, X2, Y2, X3, Y3
    % distance between the non-target and target, e.g., distractor minus
    % target (not vice versa!!!!!)

% distractor location
ntarx = squeeze(stats2.xntar_pix(ss,idx,:))'; % 3 nontargets x N trials
ntary = squeeze(stats2.yntar_pix(ss,idx,:))';
% distractor location minus target location
ntarx = ntarx - tarx; ntary = ntary - tary;
% concatenate distractor locations
testdata.distractors = [ntarx(1,:);ntary(1,:);ntarx(2,:);ntary(2,:);ntarx(3,:);ntary(3,:)];

% testdata.errors: 2 x N matrix
    % two rows: X1, Y1
    % distance between the response and target, e.g., response minus target
    % (not vice versa!!!!!)

    % get endpoint fixation position across trials
respx = stats2.p_xend_pix(ss,idx); % 1 x N trials
respy = stats2.p_yend_pix(ss,idx);
% response minus target location
respx = respx - tarx; respy = respy - tary;
% concatenate responses
testdata.errors = [respx;respy];

% testdata.dimensions = [x_resolution y_resolution]; % in pixel
testdata.dimensions = [1280 1024]; % adjust this if needed

% Initialize the model
%model = WithBias2D(SwapModel2D());
%model = SwapModel2D();
model = StandardMixtureModel2D();

SD = mean(std(testdata.errors,0,2));
MU = mean(testdata.errors,2);
MUX= MU(1);
MUY = MU(2);

model_init = AddInit(model,SD,MUX,MUY); % add initial parameters
%model_init = model;

% fitting results
% Output:
% fit: {'g', 'B', 'sd', ,'muX', 'muY'}
% llh: log likelihood value
[fit,llh] = MLE(testdata, model);
% stats.model.mu = [fit(4) nan;fit(5) nan]; % 2 dimenions (x and y) x 2 types of items (target and nontarget)
% stats.model.Sigma(:,:,1) = [fit(3) 0;0 fit(3)]; % sigma for target
% stats.model.Sigma(:,:,2) = stats.model.Sigma(:,:,1); % sigma for non-target - they are the same
% stats.model.w = [1-fit(1)-fit(2) fit(2) fit(1)]; % alpha (target prob), beta (swap rate), gamma (guess rate)
% stats.model.range = testdata.dimensions;
% stats.model.llh = llh;
% stats.model.MemModel = model;
%stats2.model_standard.mux(ss,cc,pp) = fit(4);
%stats2.model_standard.muy(ss,cc,pp)= fit(5);
%stats2.model_standard.mu(ss,cc,pp) = sqrt(fit(4)^2 + fit(5)^2);
stats2.model_standard.sd(ss,cc,pp) = fit(2); % standard deviation
stats2.model_standard.alpha(ss,cc,pp) = 1-fit(1);
%stats2.model_standard.beta(ss,cc,pp) = fit(2);
stats2.model_standard.gamma(ss,cc,pp) = fit(1);
stats2.model_standard.llh(ss,cc,pp) = llh;
stats2.model_standard.model{ss,cc,pp} = model;
stats2.model_standard.data{ss,cc,pp} = testdata;
% Calculate goodness of fit
k = length(model.upperbound); % how many parameters
dataLen = size(respx,2); % how many trials/data points
stats2.model_standard.AIC(ss,cc,pp) = -2*llh + 2 * k; % 2 parameters
stats2.model_standard.AICc(ss,cc,pp) = -2*llh + 2 * k * (dataLen / (dataLen - k - 1)); % 2 parameters
stats2.model_standard.BIC(ss,cc,pp) = -2*llh + (log(dataLen)+log(2*pi))*k; % 2 parameters
        end
    end
end

%% Swap model: the v3_short, long delay condition - serial position by cue
% only guessing, swapping, and memory imprecision

%stats1.model_swap.mux = nan(length(stats1.subs_id),2,4);
%stats1.model_swap.muy = nan(length(stats1.subs_id),2,4);
%stats1.model_swap.mu = nan(length(stats1.subs_id),2,4);
stats1.model_swap.sd = nan(length(stats1.subs_id),2,4);
stats1.model_swap.beta = nan(length(stats1.subs_id),2,4);
stats1.model_swap.gamma = nan(length(stats1.subs_id),2,4);
stats1.model_swap.alpha = nan(length(stats1.subs_id),2,4);
stats1.model_swap.llh = nan(length(stats1.subs_id),2,4);
stats1.model_swap.model = cell(length(stats1.subs_id),2,4);
stats1.model_swap.data = cell(length(stats1.subs_id),2,4);
stats1.model_swap.AIC = nan(length(stats1.subs_id),2,4);
stats1.model_swap.AICc = nan(length(stats1.subs_id),2,4);
stats1.model_swap.BIC = nan(length(stats1.subs_id),2,4);

for ss = 1:length(stats1.subs_id) % for each participant
    % Let's model quadrant and order cue differently
    for cc = 1:2 % 1, quadrant cue (0); 2, order cue (1)
        for pp = 1:4 % 4 serial positions

% row (trial) index of the specific type of cue
idx = find(stats1.cue(ss,:) == cc-1 & ...
    stats1.kept_trial(ss,:) == 1 & ...
    stats1.ord(ss,:) == pp);

if isempty(idx) | length(idx) == 1 % if only 1 trial or if there is no trial
    continue; % we do not do model
end

% what we need:
% X and Y of primary saccade endpoints
% X and Y of targets
% X and Y of non-targets

% Input:
% testdata.targets: 2 x N matrix, where N is the number of trials
    % two rows: X1, Y1

tarx = stats1.xtar_pix(ss,idx); % 1 x N trials
tary = stats1.ytar_pix(ss,idx); % 1 x N trials
testdata.targets = [tarx;tary];

% testdata.distractors: 6 x N matrix, if we have 3 distractors (3 x 2 = 6)
    % six rows: X1, Y1, X2, Y2, X3, Y3
    % distance between the non-target and target, e.g., distractor minus
    % target (not vice versa!!!!!)

% distractor location
ntarx = squeeze(stats1.xntar_pix(ss,idx,:))'; % 3 nontargets x N trials
ntary = squeeze(stats1.yntar_pix(ss,idx,:))';
% distractor location minus target location
ntarx = ntarx - tarx; ntary = ntary - tary;
% concatenate distractor locations
if length(idx) > 1 % if more than 1 trial
testdata.distractors = [ntarx(1,:);ntary(1,:);ntarx(2,:);ntary(2,:);ntarx(3,:);ntary(3,:)];
else
    testdata.distractors = [ntarx(1);ntary(1);ntarx(2);ntary(2);ntarx(3);ntary(3)];
end
% testdata.errors: 2 x N matrix
    % two rows: X1, Y1
    % distance between the response and target, e.g., response minus target
    % (not vice versa!!!!!)

    % get endpoint fixation position across trials
respx = stats1.p_xend_pix(ss,idx); % 1 x N trials
respy = stats1.p_yend_pix(ss,idx);
% response minus target location
respx = respx - tarx; respy = respy - tary;
% concatenate responses
testdata.errors = [respx;respy];

% testdata.dimensions = [x_resolution y_resolution]; % in pixel
testdata.dimensions = [1280 1024]; % adjust this if needed

% Initialize the model
%model = WithBias2D(SwapModel2D());
model = SwapModel2D();

SD = mean(std(testdata.errors,0,2));
MU = mean(testdata.errors,2);
MUX= MU(1);
MUY = MU(2);

model_init = AddInit(model,SD,MUX,MUY); % add initial parameters
%model_init = model;

% fitting results
% Output:
% fit: {'g', 'B', 'sd', ,'muX', 'muY'}
% llh: log likelihood value
[fit,llh] = MLE(testdata, model);
% stats.model.mu = [fit(4) nan;fit(5) nan]; % 2 dimenions (x and y) x 2 types of items (target and nontarget)
% stats.model.Sigma(:,:,1) = [fit(3) 0;0 fit(3)]; % sigma for target
% stats.model.Sigma(:,:,2) = stats.model.Sigma(:,:,1); % sigma for non-target - they are the same
% stats.model.w = [1-fit(1)-fit(2) fit(2) fit(1)]; % alpha (target prob), beta (swap rate), gamma (guess rate)
% stats.model.range = testdata.dimensions;
% stats.model.llh = llh;
% stats.model.MemModel = model;
%stats1.model_swap.mux(ss,cc,pp) = fit(4);
%stats1.model_swap.muy(ss,cc,pp)= fit(5);
%stats1.model_swap.mu(ss,cc,pp) = sqrt(fit(4)^2 + fit(5)^2);
stats1.model_swap.sd(ss,cc,pp) = fit(3); % standard deviation
stats1.model_swap.alpha(ss,cc,pp) = 1-fit(1)-fit(2);
stats1.model_swap.beta(ss,cc,pp) = fit(2);
stats1.model_swap.gamma(ss,cc,pp) = fit(1);
stats1.model_swap.llh(ss,cc,pp) = llh;
stats1.model_swap.model{ss,cc,pp} = model;
stats1.model_swap.data{ss,cc,pp} = testdata;
% Calculate goodness of fit
k = length(model.upperbound); % how many parameters
dataLen = size(respx,2); % how many trials/data points
stats1.model_swap.AIC(ss,cc,pp) = -2*llh + 2 * k; % 2 parameters
stats1.model_swap.AICc(ss,cc,pp) = -2*llh + 2 * k * (dataLen / (dataLen - k - 1)); % 2 parameters
stats1.model_swap.BIC(ss,cc,pp) = -2*llh + (log(dataLen)+log(2*pi))*k; % 2 parameters

        end
    end
end

%% Swap model: the v3_short_DELAY, short delay condition - serial position by cue
% only guessing, swapping, and memory imprecision

%stats2.model_swap.mux = nan(length(stats2.subs_id),2,4);
%stats2.model_swap.muy = nan(length(stats2.subs_id),2,4);
%stats2.model_swap.mu = nan(length(stats2.subs_id),2,4);
stats2.model_swap.sd = nan(length(stats2.subs_id),2,4);
stats2.model_swap.beta = nan(length(stats2.subs_id),2,4);
stats2.model_swap.gamma = nan(length(stats2.subs_id),2,4);
stats2.model_swap.alpha = nan(length(stats2.subs_id),2,4);
stats2.model_swap.llh = nan(length(stats2.subs_id),2,4);
stats2.model_swap.model = cell(length(stats2.subs_id),2,4);
stats2.model_swap.data = cell(length(stats2.subs_id),2,4);
stats2.model_swap.AIC = nan(length(stats2.subs_id),2,4);
stats2.model_swap.AICc = nan(length(stats2.subs_id),2,4);
stats2.model_swap.BIC = nan(length(stats2.subs_id),2,4);

for ss = 1:length(stats2.subs_id) % for each participant
    % Let's model quadrant and order cue differently
    for cc = 1:2 % 1, quadrant cue (0); 2, order cue (1)
        for pp = 1:4 % 4 serial positions

% row (trial) index of the specific type of cue
idx = find(stats2.cue(ss,:) == cc-1 & ...
    stats2.kept_trial(ss,:) == 1 & ...
    stats2.ord(ss,:) == pp);

if isempty(idx) | length(idx) == 1 % if only 1 trial or if there is no trial
    continue; % we do not do model
end


% what we need:
% X and Y of primary saccade endpoints
% X and Y of targets
% X and Y of non-targets

% Input:
% testdata.targets: 2 x N matrix, where N is the number of trials
    % two rows: X1, Y1

tarx = stats2.xtar_pix(ss,idx); % 1 x N trials
tary = stats2.ytar_pix(ss,idx); % 1 x N trials
testdata.targets = [tarx;tary];

% testdata.distractors: 6 x N matrix, if we have 3 distractors (3 x 2 = 6)
    % six rows: X1, Y1, X2, Y2, X3, Y3
    % distance between the non-target and target, e.g., distractor minus
    % target (not vice versa!!!!!)

% distractor location
ntarx = squeeze(stats2.xntar_pix(ss,idx,:))'; % 3 nontargets x N trials
ntary = squeeze(stats2.yntar_pix(ss,idx,:))';
% distractor location minus target location
ntarx = ntarx - tarx; ntary = ntary - tary;
% concatenate distractor locations
testdata.distractors = [ntarx(1,:);ntary(1,:);ntarx(2,:);ntary(2,:);ntarx(3,:);ntary(3,:)];

% testdata.errors: 2 x N matrix
    % two rows: X1, Y1
    % distance between the response and target, e.g., response minus target
    % (not vice versa!!!!!)

    % get endpoint fixation position across trials
respx = stats2.p_xend_pix(ss,idx); % 1 x N trials
respy = stats2.p_yend_pix(ss,idx);
% response minus target location
respx = respx - tarx; respy = respy - tary;
% concatenate responses
testdata.errors = [respx;respy];

% testdata.dimensions = [x_resolution y_resolution]; % in pixel
testdata.dimensions = [1280 1024]; % adjust this if needed

% Initialize the model
%model = WithBias2D(SwapModel2D());
model = SwapModel2D();

SD = mean(std(testdata.errors,0,2));
MU = mean(testdata.errors,2);
MUX= MU(1);
MUY = MU(2);

model_init = AddInit(model,SD,MUX,MUY); % add initial parameters
%model_init = model;

% fitting results
% Output:
% fit: {'g', 'B', 'sd', ,'muX', 'muY'}
% llh: log likelihood value
[fit,llh] = MLE(testdata, model);
% stats.model.mu = [fit(4) nan;fit(5) nan]; % 2 dimenions (x and y) x 2 types of items (target and nontarget)
% stats.model.Sigma(:,:,1) = [fit(3) 0;0 fit(3)]; % sigma for target
% stats.model.Sigma(:,:,2) = stats.model.Sigma(:,:,1); % sigma for non-target - they are the same
% stats.model.w = [1-fit(1)-fit(2) fit(2) fit(1)]; % alpha (target prob), beta (swap rate), gamma (guess rate)
% stats.model.range = testdata.dimensions;
% stats.model.llh = llh;
% stats.model.MemModel = model;
%stats2.model_swap.mux(ss,cc,pp) = fit(4);
%stats2.model_swap.muy(ss,cc,pp)= fit(5);
%stats2.model_swap.mu(ss,cc,pp) = sqrt(fit(4)^2 + fit(5)^2);
stats2.model_swap.sd(ss,cc,pp) = fit(3); % standard deviation
stats2.model_swap.alpha(ss,cc,pp) = 1-fit(1)-fit(2);
stats2.model_swap.beta(ss,cc,pp) = fit(2);
stats2.model_swap.gamma(ss,cc,pp) = fit(1);
stats2.model_swap.llh(ss,cc,pp) = llh;
stats2.model_swap.model{ss,cc,pp} = model;
stats2.model_swap.data{ss,cc,pp} = testdata;
% Calculate goodness of fit
k = length(model.upperbound); % how many parameters
dataLen = size(respx,2); % how many trials/data points
stats2.model_swap.AIC(ss,cc,pp) = -2*llh + 2 * k; % 2 parameters
stats2.model_swap.AICc(ss,cc,pp) = -2*llh + 2 * k * (dataLen / (dataLen - k - 1)); % 2 parameters
stats2.model_swap.BIC(ss,cc,pp) = -2*llh + (log(dataLen)+log(2*pi))*k; % 2 parameters

        end
    end
end

% %% Swap with response sampling: the v3_short, long delay condition - serial position by cue
% % only guessing, swapping, and memory imprecision
% 
% %stats1.model_swap_samp.mux = nan(length(stats1.subs_id),2,4);
% %stats1.model_swap_samp.muy = nan(length(stats1.subs_id),2,4);
% %stats1.model_swap_samp.mu = nan(length(stats1.subs_id),2,4);
% stats1.model_swap_samp.sd = nan(length(stats1.subs_id),2,4);
% stats1.model_swap_samp.beta = nan(length(stats1.subs_id),2,4);
% stats1.model_swap_samp.gamma = nan(length(stats1.subs_id),2,4);
% stats1.model_swap_samp.alpha = nan(length(stats1.subs_id),2,4);
% stats1.model_swap_samp.llh = nan(length(stats1.subs_id),2,4);
% 
% for ss = 1:length(stats1.subs_id) % for each participant
%     % Let's model quadrant and order cue differently
%     for cc = 1:2 % 1, quadrant cue (0); 2, order cue (1)
%         for pp = 1:4 % 4 serial positions
% 
% % row (trial) index of the specific type of cue
% idx = find(stats1.cue(ss,[1:64]+64*(cc-1)) == cc-1 & ...
%     stats1.kept_trial(ss,[1:64]+64*(cc-1)) == 1 & ...
%     stats1.ord(ss,[1:64]+64*(cc-1)) == pp);
% 
% if isempty(idx) | length(idx) == 1 % if only 1 trial or if there is no trial
%     continue; % we do not do model
% end
% 
% % what we need:
% % X and Y of primary saccade endpoints
% % X and Y of targets
% % X and Y of non-targets
% 
% % Input:
% % testdata.targets: 2 x N matrix, where N is the number of trials
%     % two rows: X1, Y1
% 
% tarx = stats1.xtar_pix(ss,idx); % 1 x N trials
% tary = stats1.ytar_pix(ss,idx); % 1 x N trials
% testdata.targets = [tarx;tary];
% 
% % testdata.distractors: 6 x N matrix, if we have 3 distractors (3 x 2 = 6)
%     % six rows: X1, Y1, X2, Y2, X3, Y3
%     % distance between the non-target and target, e.g., distractor minus
%     % target (not vice versa!!!!!)
% 
% % distractor location
% ntarx = squeeze(stats1.xntar_pix(ss,idx,:))'; % 3 nontargets x N trials
% ntary = squeeze(stats1.yntar_pix(ss,idx,:))';
% % distractor location minus target location
% ntarx = ntarx - tarx; ntary = ntary - tary;
% % concatenate distractor locations
% if length(idx) > 1 % if more than 1 trial
% testdata.distractors = [ntarx(1,:);ntary(1,:);ntarx(2,:);ntary(2,:);ntarx(3,:);ntary(3,:)];
% else
%     testdata.distractors = [ntarx(1);ntary(1);ntarx(2);ntary(2);ntarx(3);ntary(3)];
% end
% % testdata.errors: 2 x N matrix
%     % two rows: X1, Y1
%     % distance between the response and target, e.g., response minus target
%     % (not vice versa!!!!!)
% 
%     % get endpoint fixation position across trials
% respx = stats1.p_xend_pix(ss,idx); % 1 x N trials
% respy = stats1.p_yend_pix(ss,idx);
% % response minus target location
% respx = respx - tarx; respy = respy - tary;
% % concatenate responses
% testdata.errors = [respx;respy];
% 
% % testdata.dimensions = [x_resolution y_resolution]; % in pixel
% testdata.dimensions = [1280 1024]; % adjust this if needed
% 
% % Initialize the model
% %model = WithBias2D(SwapModel2D());
% %model = SwapModel2D();
% model = WithResponseSampling2D(SwapModel2D());
% 
% SD = mean(std(testdata.errors,0,2));
% MU = mean(testdata.errors,2);
% MUX= MU(1);
% MUY = MU(2);
% 
% model_init = AddInit(model,SD,MUX,MUY); % add initial parameters
% %model_init = model;
% 
% % fitting results
% % Output:
% % fit: {'g', 'B', 'sd', ,'muX', 'muY'}
% % llh: log likelihood value
% [fit,llh] = MLE(testdata, model);
% % stats.model.mu = [fit(4) nan;fit(5) nan]; % 2 dimenions (x and y) x 2 types of items (target and nontarget)
% % stats.model.Sigma(:,:,1) = [fit(3) 0;0 fit(3)]; % sigma for target
% % stats.model.Sigma(:,:,2) = stats.model.Sigma(:,:,1); % sigma for non-target - they are the same
% % stats.model.w = [1-fit(1)-fit(2) fit(2) fit(1)]; % alpha (target prob), beta (swap rate), gamma (guess rate)
% % stats.model.range = testdata.dimensions;
% % stats.model.llh = llh;
% % stats.model.MemModel = model;
% %stats1.model_swap_samp.mux(ss,cc,pp) = fit(4);
% %stats1.model_swap_samp.muy(ss,cc,pp)= fit(5);
% %stats1.model_swap_samp.mu(ss,cc,pp) = sqrt(fit(4)^2 + fit(5)^2);
% stats1.model_swap_samp.sd(ss,cc,pp) = fit(3); % standard deviation
% stats1.model_swap_samp.alpha(ss,cc,pp) = 1-fit(1)-fit(2);
% stats1.model_swap_samp.beta(ss,cc,pp) = fit(2);
% stats1.model_swap_samp.gamma(ss,cc,pp) = fit(1);
% stats1.model_swap_samp.llh(ss,cc,pp) = llh;
%         end
%     end
% end
% 
% %% Swap with response sampling: the v3_short_DELAY, short delay condition - serial position by cue
% % only guessing, swapping, and memory imprecision
% 
% %stats2.model_swap_samp.mux = nan(length(stats2.subs_id),2,4);
% %stats2.model_swap_samp.muy = nan(length(stats2.subs_id),2,4);
% %stats2.model_swap_samp.mu = nan(length(stats2.subs_id),2,4);
% stats2.model_swap_samp.sd = nan(length(stats2.subs_id),2,4);
% stats2.model_swap_samp.beta = nan(length(stats2.subs_id),2,4);
% stats2.model_swap_samp.gamma = nan(length(stats2.subs_id),2,4);
% stats2.model_swap_samp.alpha = nan(length(stats2.subs_id),2,4);
% stats2.model_swap_samp.llh = nan(length(stats2.subs_id),2,4);
% 
% for ss = 1:length(stats2.subs_id) % for each participant
%     % Let's model quadrant and order cue differently
%     for cc = 1:2 % 1, quadrant cue (0); 2, order cue (1)
%         for pp = 1:4 % 4 serial positions
% 
% % row (trial) index of the specific type of cue
% idx = find(stats2.cue(ss,[1:64]+64*(cc-1)) == cc-1 & ...
%     stats2.kept_trial(ss,[1:64]+64*(cc-1)) == 1 & ...
%     stats2.ord(ss,[1:64]+64*(cc-1)) == pp);
% 
% if isempty(idx) | length(idx) == 1 % if only 1 trial or if there is no trial
%     continue; % we do not do model
% end
% 
% 
% % what we need:
% % X and Y of primary saccade endpoints
% % X and Y of targets
% % X and Y of non-targets
% 
% % Input:
% % testdata.targets: 2 x N matrix, where N is the number of trials
%     % two rows: X1, Y1
% 
% tarx = stats2.xtar_pix(ss,idx); % 1 x N trials
% tary = stats2.ytar_pix(ss,idx); % 1 x N trials
% testdata.targets = [tarx;tary];
% 
% % testdata.distractors: 6 x N matrix, if we have 3 distractors (3 x 2 = 6)
%     % six rows: X1, Y1, X2, Y2, X3, Y3
%     % distance between the non-target and target, e.g., distractor minus
%     % target (not vice versa!!!!!)
% 
% % distractor location
% ntarx = squeeze(stats2.xntar_pix(ss,idx,:))'; % 3 nontargets x N trials
% ntary = squeeze(stats2.yntar_pix(ss,idx,:))';
% % distractor location minus target location
% ntarx = ntarx - tarx; ntary = ntary - tary;
% % concatenate distractor locations
% testdata.distractors = [ntarx(1,:);ntary(1,:);ntarx(2,:);ntary(2,:);ntarx(3,:);ntary(3,:)];
% 
% % testdata.errors: 2 x N matrix
%     % two rows: X1, Y1
%     % distance between the response and target, e.g., response minus target
%     % (not vice versa!!!!!)
% 
%     % get endpoint fixation position across trials
% respx = stats2.p_xend_pix(ss,idx); % 1 x N trials
% respy = stats2.p_yend_pix(ss,idx);
% % response minus target location
% respx = respx - tarx; respy = respy - tary;
% % concatenate responses
% testdata.errors = [respx;respy];
% 
% % testdata.dimensions = [x_resolution y_resolution]; % in pixel
% testdata.dimensions = [1280 1024]; % adjust this if needed
% 
% % Initialize the model
% %model = WithBias2D(SwapModel2D());
% %model = SwapModel2D();
% model = WithResponseSampling2D(SwapModel2D());
% 
% SD = mean(std(testdata.errors,0,2));
% MU = mean(testdata.errors,2);
% MUX= MU(1);
% MUY = MU(2);
% 
% model_init = AddInit(model,SD,MUX,MUY); % add initial parameters
% %model_init = model;
% 
% % fitting results
% % Output:
% % fit: {'g', 'B', 'sd', ,'muX', 'muY'}
% % llh: log likelihood value
% [fit,llh] = MLE(testdata, model);
% % stats.model.mu = [fit(4) nan;fit(5) nan]; % 2 dimenions (x and y) x 2 types of items (target and nontarget)
% % stats.model.Sigma(:,:,1) = [fit(3) 0;0 fit(3)]; % sigma for target
% % stats.model.Sigma(:,:,2) = stats.model.Sigma(:,:,1); % sigma for non-target - they are the same
% % stats.model.w = [1-fit(1)-fit(2) fit(2) fit(1)]; % alpha (target prob), beta (swap rate), gamma (guess rate)
% % stats.model.range = testdata.dimensions;
% % stats.model.llh = llh;
% % stats.model.MemModel = model;
% %stats2.model_swap_samp.mux(ss,cc,pp) = fit(4);
% %stats2.model_swap_samp.muy(ss,cc,pp)= fit(5);
% %stats2.model_swap_samp.mu(ss,cc,pp) = sqrt(fit(4)^2 + fit(5)^2);
% stats2.model_swap_samp.sd(ss,cc,pp) = fit(3); % standard deviation
% stats2.model_swap_samp.alpha(ss,cc,pp) = 1-fit(1)-fit(2);
% stats2.model_swap_samp.beta(ss,cc,pp) = fit(2);
% stats2.model_swap_samp.gamma(ss,cc,pp) = fit(1);
% stats2.model_swap_samp.llh(ss,cc,pp) = llh;
%         end
%     end
% end
% 
%% Swap with systematic bias: the v3_short, long delay condition - serial position by cue
% only guessing, swapping, and memory imprecision

stats1.model_swap_sys.mux = nan(length(stats1.subs_id),2,4);
stats1.model_swap_sys.muy = nan(length(stats1.subs_id),2,4);
stats1.model_swap_sys.mu = nan(length(stats1.subs_id),2,4);
stats1.model_swap_sys.sd = nan(length(stats1.subs_id),2,4);
stats1.model_swap_sys.beta = nan(length(stats1.subs_id),2,4);
stats1.model_swap_sys.gamma = nan(length(stats1.subs_id),2,4);
stats1.model_swap_sys.alpha = nan(length(stats1.subs_id),2,4);
stats1.model_swap_sys.llh = nan(length(stats1.subs_id),2,4);
stats1.model_swap_sys.model = cell(length(stats1.subs_id),2,4);
stats1.model_swap_sys.data = cell(length(stats1.subs_id),2,4);
stats1.model_swap_sys.AIC = nan(length(stats1.subs_id),2,4);
stats1.model_swap_sys.AICc = nan(length(stats1.subs_id),2,4);
stats1.model_swap_sys.BIC = nan(length(stats1.subs_id),2,4);

for ss = 1:length(stats1.subs_id) % for each participant
    % Let's model quadrant and order cue differently
    for cc = 1:2 % 1, quadrant cue (0); 2, order cue (1)
        for pp = 1:4 % 4 serial positions

% row (trial) index of the specific type of cue
idx = find(stats1.cue(ss,:) == cc-1 & ...
    stats1.kept_trial(ss,:) == 1 & ...
    stats1.ord(ss,:) == pp);

if isempty(idx) | length(idx) == 1 % if only 1 trial or if there is no trial
    continue; % we do not do model
end

% what we need:
% X and Y of primary saccade endpoints
% X and Y of targets
% X and Y of non-targets

% Input:
% testdata.targets: 2 x N matrix, where N is the number of trials
    % two rows: X1, Y1

tarx = stats1.xtar_pix(ss,idx); % 1 x N trials
tary = stats1.ytar_pix(ss,idx); % 1 x N trials
testdata.targets = [tarx;tary];

% testdata.distractors: 6 x N matrix, if we have 3 distractors (3 x 2 = 6)
    % six rows: X1, Y1, X2, Y2, X3, Y3
    % distance between the non-target and target, e.g., distractor minus
    % target (not vice versa!!!!!)

% distractor location
ntarx = squeeze(stats1.xntar_pix(ss,idx,:))'; % 3 nontargets x N trials
ntary = squeeze(stats1.yntar_pix(ss,idx,:))';
% distractor location minus target location
ntarx = ntarx - tarx; ntary = ntary - tary;
% concatenate distractor locations
if length(idx) > 1 % if more than 1 trial
testdata.distractors = [ntarx(1,:);ntary(1,:);ntarx(2,:);ntary(2,:);ntarx(3,:);ntary(3,:)];
else
    testdata.distractors = [ntarx(1);ntary(1);ntarx(2);ntary(2);ntarx(3);ntary(3)];
end
% testdata.errors: 2 x N matrix
    % two rows: X1, Y1
    % distance between the response and target, e.g., response minus target
    % (not vice versa!!!!!)

    % get endpoint fixation position across trials
respx = stats1.p_xend_pix(ss,idx); % 1 x N trials
respy = stats1.p_yend_pix(ss,idx);
% response minus target location
respx = respx - tarx; respy = respy - tary;
% concatenate responses
testdata.errors = [respx;respy];

% testdata.dimensions = [x_resolution y_resolution]; % in pixel
testdata.dimensions = [1280 1024]; % adjust this if needed

% Initialize the model
model = WithBias2D(SwapModel2D());
%model = SwapModel2D();
%model = WithResponseSampling2D(SwapModel2D());

SD = mean(std(testdata.errors,0,2));
MU = mean(testdata.errors,2);
MUX= MU(1);
MUY = MU(2);

model_init = AddInit(model,SD,MUX,MUY); % add initial parameters
%model_init = model;

% fitting results
% Output:
% fit: {'g', 'B', 'sd', ,'muX', 'muY'}
% llh: log likelihood value
[fit,llh] = MLE(testdata, model);
% stats.model.mu = [fit(4) nan;fit(5) nan]; % 2 dimenions (x and y) x 2 types of items (target and nontarget)
% stats.model.Sigma(:,:,1) = [fit(3) 0;0 fit(3)]; % sigma for target
% stats.model.Sigma(:,:,2) = stats.model.Sigma(:,:,1); % sigma for non-target - they are the same
% stats.model.w = [1-fit(1)-fit(2) fit(2) fit(1)]; % alpha (target prob), beta (swap rate), gamma (guess rate)
% stats.model.range = testdata.dimensions;
% stats.model.llh = llh;
% stats.model.MemModel = model;
stats1.model_swap_sys.mux(ss,cc,pp) = fit(4);
stats1.model_swap_sys.muy(ss,cc,pp)= fit(5);
stats1.model_swap_sys.mu(ss,cc,pp) = sqrt(fit(4)^2 + fit(5)^2);
stats1.model_swap_sys.sd(ss,cc,pp) = fit(3); % standard deviation
stats1.model_swap_sys.alpha(ss,cc,pp) = 1-fit(1)-fit(2);
stats1.model_swap_sys.beta(ss,cc,pp) = fit(2);
stats1.model_swap_sys.gamma(ss,cc,pp) = fit(1);
stats1.model_swap_sys.llh(ss,cc,pp) = llh;
stats1.model_swap_sys.model{ss,cc,pp} = model;
stats1.model_swap_sys.data{ss,cc,pp} = testdata;
% Calculate goodness of fit
k = length(model.upperbound); % how many parameters
dataLen = size(respx,2); % how many trials/data points
stats1.model_swap_sys.AIC(ss,cc,pp) = -2*llh + 2 * k; % 2 parameters
stats1.model_swap_sys.AICc(ss,cc,pp) = -2*llh + 2 * k * (dataLen / (dataLen - k - 1)); % 2 parameters
stats1.model_swap_sys.BIC(ss,cc,pp) = -2*llh + (log(dataLen)+log(2*pi))*k; % 2 parameters

        end
    end
end

%% Swap with systematic bias: the v3_short_DELAY, short delay condition - serial position by cue
% only guessing, swapping, and memory imprecision

stats2.model_swap_sys.mux = nan(length(stats2.subs_id),2,4);
stats2.model_swap_sys.muy = nan(length(stats2.subs_id),2,4);
stats2.model_swap_sys.mu = nan(length(stats2.subs_id),2,4);
stats2.model_swap_sys.sd = nan(length(stats2.subs_id),2,4);
stats2.model_swap_sys.beta = nan(length(stats2.subs_id),2,4);
stats2.model_swap_sys.gamma = nan(length(stats2.subs_id),2,4);
stats2.model_swap_sys.alpha = nan(length(stats2.subs_id),2,4);
stats2.model_swap_sys.llh = nan(length(stats2.subs_id),2,4);
stats2.model_swap_sys.model = cell(length(stats2.subs_id),2,4);
stats2.model_swap_sys.data = cell(length(stats2.subs_id),2,4);
stats2.model_swap_sys.AIC = nan(length(stats2.subs_id),2,4);
stats2.model_swap_sys.AICc = nan(length(stats2.subs_id),2,4);
stats2.model_swap_sys.BIC = nan(length(stats2.subs_id),2,4);

for ss = 1:length(stats2.subs_id) % for each participant
    % Let's model quadrant and order cue differently
    for cc = 1:2 % 1, quadrant cue (0); 2, order cue (1)
        for pp = 1:4 % 4 serial positions

% row (trial) index of the specific type of cue
idx = find(stats2.cue(ss,:) == cc-1 & ...
    stats2.kept_trial(ss,:) == 1 & ...
    stats2.ord(ss,:) == pp);

if isempty(idx) | length(idx) == 1 % if only 1 trial or if there is no trial
    continue; % we do not do model
end


% what we need:
% X and Y of primary saccade endpoints
% X and Y of targets
% X and Y of non-targets

% Input:
% testdata.targets: 2 x N matrix, where N is the number of trials
    % two rows: X1, Y1

tarx = stats2.xtar_pix(ss,idx); % 1 x N trials
tary = stats2.ytar_pix(ss,idx); % 1 x N trials
testdata.targets = [tarx;tary];

% testdata.distractors: 6 x N matrix, if we have 3 distractors (3 x 2 = 6)
    % six rows: X1, Y1, X2, Y2, X3, Y3
    % distance between the non-target and target, e.g., distractor minus
    % target (not vice versa!!!!!)

% distractor location
ntarx = squeeze(stats2.xntar_pix(ss,idx,:))'; % 3 nontargets x N trials
ntary = squeeze(stats2.yntar_pix(ss,idx,:))';
% distractor location minus target location
ntarx = ntarx - tarx; ntary = ntary - tary;
% concatenate distractor locations
testdata.distractors = [ntarx(1,:);ntary(1,:);ntarx(2,:);ntary(2,:);ntarx(3,:);ntary(3,:)];

% testdata.errors: 2 x N matrix
    % two rows: X1, Y1
    % distance between the response and target, e.g., response minus target
    % (not vice versa!!!!!)

    % get endpoint fixation position across trials
respx = stats2.p_xend_pix(ss,idx); % 1 x N trials
respy = stats2.p_yend_pix(ss,idx);
% response minus target location
respx = respx - tarx; respy = respy - tary;
% concatenate responses
testdata.errors = [respx;respy];

% testdata.dimensions = [x_resolution y_resolution]; % in pixel
testdata.dimensions = [1280 1024]; % adjust this if needed

% Initialize the model
model = WithBias2D(SwapModel2D());
%model = SwapModel2D();
%model = WithResponseSampling2D(SwapModel2D());

SD = mean(std(testdata.errors,0,2));
MU = mean(testdata.errors,2);
MUX= MU(1);
MUY = MU(2);

model_init = AddInit(model,SD,MUX,MUY); % add initial parameters
%model_init = model;

% fitting results
% Output:
% fit: {'g', 'B', 'sd', ,'muX', 'muY'}
% llh: log likelihood value
[fit,llh] = MLE(testdata, model);
% stats.model.mu = [fit(4) nan;fit(5) nan]; % 2 dimenions (x and y) x 2 types of items (target and nontarget)
% stats.model.Sigma(:,:,1) = [fit(3) 0;0 fit(3)]; % sigma for target
% stats.model.Sigma(:,:,2) = stats.model.Sigma(:,:,1); % sigma for non-target - they are the same
% stats.model.w = [1-fit(1)-fit(2) fit(2) fit(1)]; % alpha (target prob), beta (swap rate), gamma (guess rate)
% stats.model.range = testdata.dimensions;
% stats.model.llh = llh;
% stats.model.MemModel = model;
stats2.model_swap_sys.mux(ss,cc,pp) = fit(4);
stats2.model_swap_sys.muy(ss,cc,pp)= fit(5);
stats2.model_swap_sys.mu(ss,cc,pp) = sqrt(fit(4)^2 + fit(5)^2);
stats2.model_swap_sys.sd(ss,cc,pp) = fit(3); % standard deviation
stats2.model_swap_sys.alpha(ss,cc,pp) = 1-fit(1)-fit(2);
stats2.model_swap_sys.beta(ss,cc,pp) = fit(2);
stats2.model_swap_sys.gamma(ss,cc,pp) = fit(1);
stats2.model_swap_sys.llh(ss,cc,pp) = llh;
stats2.model_swap_sys.model{ss,cc,pp} = model;
stats2.model_swap_sys.data{ss,cc,pp} = testdata;
% Calculate goodness of fit
k = length(model.upperbound); % how many parameters
dataLen = size(respx,2); % how many trials/data points
stats2.model_swap_sys.AIC(ss,cc,pp) = -2*llh + 2 * k; % 2 parameters
stats2.model_swap_sys.AICc(ss,cc,pp) = -2*llh + 2 * k * (dataLen / (dataLen - k - 1)); % 2 parameters
stats2.model_swap_sys.BIC(ss,cc,pp) = -2*llh + (log(dataLen)+log(2*pi))*k; % 2 parameters

        end
    end
end

% %% Swap with systematic bias and response sampling the v3_short, long delay condition - serial position by cue
% % only guessing, swapping, and memory imprecision
% 
% stats1.model_swap_sys_samp.mux = nan(length(stats1.subs_id),2,4);
% stats1.model_swap_sys_samp.muy = nan(length(stats1.subs_id),2,4);
% stats1.model_swap_sys_samp.mu = nan(length(stats1.subs_id),2,4);
% stats1.model_swap_sys_samp.sd = nan(length(stats1.subs_id),2,4);
% stats1.model_swap_sys_samp.beta = nan(length(stats1.subs_id),2,4);
% stats1.model_swap_sys_samp.gamma = nan(length(stats1.subs_id),2,4);
% stats1.model_swap_sys_samp.alpha = nan(length(stats1.subs_id),2,4);
% stats1.model_swap_sys_samp.llh = nan(length(stats1.subs_id),2,4);
% 
% for ss = 1:length(stats1.subs_id) % for each participant
%     % Let's model quadrant and order cue differently
%     for cc = 1:2 % 1, quadrant cue (0); 2, order cue (1)
%         for pp = 1:4 % 4 serial positions
% 
% % row (trial) index of the specific type of cue
% idx = find(stats1.cue(ss,[1:64]+64*(cc-1)) == cc-1 & ...
%     stats1.kept_trial(ss,[1:64]+64*(cc-1)) == 1 & ...
%     stats1.ord(ss,[1:64]+64*(cc-1)) == pp);
% 
% if isempty(idx) | length(idx) == 1 % if only 1 trial or if there is no trial
%     continue; % we do not do model
% end
% 
% % what we need:
% % X and Y of primary saccade endpoints
% % X and Y of targets
% % X and Y of non-targets
% 
% % Input:
% % testdata.targets: 2 x N matrix, where N is the number of trials
%     % two rows: X1, Y1
% 
% tarx = stats1.xtar_pix(ss,idx); % 1 x N trials
% tary = stats1.ytar_pix(ss,idx); % 1 x N trials
% testdata.targets = [tarx;tary];
% 
% % testdata.distractors: 6 x N matrix, if we have 3 distractors (3 x 2 = 6)
%     % six rows: X1, Y1, X2, Y2, X3, Y3
%     % distance between the non-target and target, e.g., distractor minus
%     % target (not vice versa!!!!!)
% 
% % distractor location
% ntarx = squeeze(stats1.xntar_pix(ss,idx,:))'; % 3 nontargets x N trials
% ntary = squeeze(stats1.yntar_pix(ss,idx,:))';
% % distractor location minus target location
% ntarx = ntarx - tarx; ntary = ntary - tary;
% % concatenate distractor locations
% if length(idx) > 1 % if more than 1 trial
% testdata.distractors = [ntarx(1,:);ntary(1,:);ntarx(2,:);ntary(2,:);ntarx(3,:);ntary(3,:)];
% else
%     testdata.distractors = [ntarx(1);ntary(1);ntarx(2);ntary(2);ntarx(3);ntary(3)];
% end
% % testdata.errors: 2 x N matrix
%     % two rows: X1, Y1
%     % distance between the response and target, e.g., response minus target
%     % (not vice versa!!!!!)
% 
%     % get endpoint fixation position across trials
% respx = stats1.p_xend_pix(ss,idx); % 1 x N trials
% respy = stats1.p_yend_pix(ss,idx);
% % response minus target location
% respx = respx - tarx; respy = respy - tary;
% % concatenate responses
% testdata.errors = [respx;respy];
% 
% % testdata.dimensions = [x_resolution y_resolution]; % in pixel
% testdata.dimensions = [1280 1024]; % adjust this if needed
% 
% % Initialize the model
% model = WithResponseSampling2D(WithBias2D(SwapModel2D()));
% %model = SwapModel2D();
% %model = WithResponseSampling2D(SwapModel2D());
% 
% SD = mean(std(testdata.errors,0,2));
% MU = mean(testdata.errors,2);
% MUX= MU(1);
% MUY = MU(2);
% 
% model_init = AddInit(model,SD,MUX,MUY); % add initial parameters
% %model_init = model;
% 
% % fitting results
% % Output:
% % fit: {'g', 'B', 'sd', ,'muX', 'muY'}
% % llh: log likelihood value
% [fit,llh] = MLE(testdata, model);
% % stats.model.mu = [fit(4) nan;fit(5) nan]; % 2 dimenions (x and y) x 2 types of items (target and nontarget)
% % stats.model.Sigma(:,:,1) = [fit(3) 0;0 fit(3)]; % sigma for target
% % stats.model.Sigma(:,:,2) = stats.model.Sigma(:,:,1); % sigma for non-target - they are the same
% % stats.model.w = [1-fit(1)-fit(2) fit(2) fit(1)]; % alpha (target prob), beta (swap rate), gamma (guess rate)
% % stats.model.range = testdata.dimensions;
% % stats.model.llh = llh;
% % stats.model.MemModel = model;
% stats1.model_swap_sys_samp.mux(ss,cc,pp) = fit(4);
% stats1.model_swap_sys_samp.muy(ss,cc,pp)= fit(5);
% stats1.model_swap_sys_samp.mu(ss,cc,pp) = sqrt(fit(4)^2 + fit(5)^2);
% stats1.model_swap_sys_samp.sd(ss,cc,pp) = fit(3); % standard deviation
% stats1.model_swap_sys_samp.alpha(ss,cc,pp) = 1-fit(1)-fit(2);
% stats1.model_swap_sys_samp.beta(ss,cc,pp) = fit(2);
% stats1.model_swap_sys_samp.gamma(ss,cc,pp) = fit(1);
% stats1.model_swap_sys_samp.llh(ss,cc,pp) = llh;
%         end
%     end
% end
% 
% %% Swap with systematic bias and response sampling: the v3_short_DELAY, short delay condition - serial position by cue
% % only guessing, swapping, and memory imprecision
% 
% stats2.model_swap_sys_samp.mux = nan(length(stats2.subs_id),2,4);
% stats2.model_swap_sys_samp.muy = nan(length(stats2.subs_id),2,4);
% stats2.model_swap_sys_samp.mu = nan(length(stats2.subs_id),2,4);
% stats2.model_swap_sys_samp.sd = nan(length(stats2.subs_id),2,4);
% stats2.model_swap_sys_samp.beta = nan(length(stats2.subs_id),2,4);
% stats2.model_swap_sys_samp.gamma = nan(length(stats2.subs_id),2,4);
% stats2.model_swap_sys_samp.alpha = nan(length(stats2.subs_id),2,4);
% stats2.model_swap_sys_samp.llh = nan(length(stats2.subs_id),2,4);
% 
% for ss = 1:length(stats2.subs_id) % for each participant
%     % Let's model quadrant and order cue differently
%     for cc = 1:2 % 1, quadrant cue (0); 2, order cue (1)
%         for pp = 1:4 % 4 serial positions
% 
% % row (trial) index of the specific type of cue
% idx = find(stats2.cue(ss,[1:64]+64*(cc-1)) == cc-1 & ...
%     stats2.kept_trial(ss,[1:64]+64*(cc-1)) == 1 & ...
%     stats2.ord(ss,[1:64]+64*(cc-1)) == pp);
% 
% if isempty(idx) | length(idx) == 1 % if only 1 trial or if there is no trial
%     continue; % we do not do model
% end
% 
% 
% % what we need:
% % X and Y of primary saccade endpoints
% % X and Y of targets
% % X and Y of non-targets
% 
% % Input:
% % testdata.targets: 2 x N matrix, where N is the number of trials
%     % two rows: X1, Y1
% 
% tarx = stats2.xtar_pix(ss,idx); % 1 x N trials
% tary = stats2.ytar_pix(ss,idx); % 1 x N trials
% testdata.targets = [tarx;tary];
% 
% % testdata.distractors: 6 x N matrix, if we have 3 distractors (3 x 2 = 6)
%     % six rows: X1, Y1, X2, Y2, X3, Y3
%     % distance between the non-target and target, e.g., distractor minus
%     % target (not vice versa!!!!!)
% 
% % distractor location
% ntarx = squeeze(stats2.xntar_pix(ss,idx,:))'; % 3 nontargets x N trials
% ntary = squeeze(stats2.yntar_pix(ss,idx,:))';
% % distractor location minus target location
% ntarx = ntarx - tarx; ntary = ntary - tary;
% % concatenate distractor locations
% testdata.distractors = [ntarx(1,:);ntary(1,:);ntarx(2,:);ntary(2,:);ntarx(3,:);ntary(3,:)];
% 
% % testdata.errors: 2 x N matrix
%     % two rows: X1, Y1
%     % distance between the response and target, e.g., response minus target
%     % (not vice versa!!!!!)
% 
%     % get endpoint fixation position across trials
% respx = stats2.p_xend_pix(ss,idx); % 1 x N trials
% respy = stats2.p_yend_pix(ss,idx);
% % response minus target location
% respx = respx - tarx; respy = respy - tary;
% % concatenate responses
% testdata.errors = [respx;respy];
% 
% % testdata.dimensions = [x_resolution y_resolution]; % in pixel
% testdata.dimensions = [1280 1024]; % adjust this if needed
% 
% % Initialize the model
% %model = WithBias2D(SwapModel2D());
% %model = SwapModel2D();
% %model = WithResponseSampling2D(SwapModel2D());
% model = WithResponseSampling2D(WithBias2D(SwapModel2D()));
% 
% SD = mean(std(testdata.errors,0,2));
% MU = mean(testdata.errors,2);
% MUX= MU(1);
% MUY = MU(2);
% 
% model_init = AddInit(model,SD,MUX,MUY); % add initial parameters
% %model_init = model;
% 
% % fitting results
% % Output:
% % fit: {'g', 'B', 'sd', ,'muX', 'muY'}
% % llh: log likelihood value
% [fit,llh] = MLE(testdata, model);
% % stats.model.mu = [fit(4) nan;fit(5) nan]; % 2 dimenions (x and y) x 2 types of items (target and nontarget)
% % stats.model.Sigma(:,:,1) = [fit(3) 0;0 fit(3)]; % sigma for target
% % stats.model.Sigma(:,:,2) = stats.model.Sigma(:,:,1); % sigma for non-target - they are the same
% % stats.model.w = [1-fit(1)-fit(2) fit(2) fit(1)]; % alpha (target prob), beta (swap rate), gamma (guess rate)
% % stats.model.range = testdata.dimensions;
% % stats.model.llh = llh;
% % stats.model.MemModel = model;
% stats2.model_swap_sys_samp.mux(ss,cc,pp) = fit(4);
% stats2.model_swap_sys_samp.muy(ss,cc,pp)= fit(5);
% stats2.model_swap_sys_samp.mu(ss,cc,pp) = sqrt(fit(4)^2 + fit(5)^2);
% stats2.model_swap_sys_samp.sd(ss,cc,pp) = fit(3); % standard deviation
% stats2.model_swap_sys_samp.alpha(ss,cc,pp) = 1-fit(1)-fit(2);
% stats2.model_swap_sys_samp.beta(ss,cc,pp) = fit(2);
% stats2.model_swap_sys_samp.gamma(ss,cc,pp) = fit(1);
% stats2.model_swap_sys_samp.llh(ss,cc,pp) = llh;
%         end
%     end
% end

%% save the stats
save(fullfile(output,'stats_exp1_longdelay.mat'),'stats1')
save(fullfile(output,'stats_exp1_shortdelay.mat'),'stats2')
