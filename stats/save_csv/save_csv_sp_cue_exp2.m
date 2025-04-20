% Save the data into csv files

% By: Linjing Jiang
% Date: 04/03/2024

%% 
clear
close all
clc

% output directory
output = './analysis/exp2';


% Load data
cd(output)
%load('data1.mat') % v3_ord, long delay
%load('data2.mat') % v3_ord_delay, short delay
%qa1 = load('qa_exp1_quaddelay.mat','qa'); qa1 = qa1.qa;
%qa2 = load('qa_exp1_orddelay.mat','qa'); qa2 = qa2.qa;
load('stats_exp2_quad.mat')
load('stats_exp2_ord.mat')

%% Primary saccade error (of ANY stimuli)

clearvars -except output stats1 stats2 qa1 qa2

data_quad = squeeze(mean(stats1.p_err_ord,2,'omitnan'));
data_ord = squeeze(mean(stats2.p_err_ord,2,'omitnan'));
data = [data_quad(:,1:4);data_ord(:,5:8)];

id = [stats1.subs_id';stats2.subs_id'];
cue = [zeros(size(stats1.subs_id))';ones(size(stats2.subs_id))']; % 0 - quadrant cue; 1 - order cue
final_sind = [stats1.final_sind;stats2.final_sind+length(stats1.subs_id)];

data = data(final_sind,:);
id = id(final_sind,:);
cue = cue(final_sind,:);

tbl = array2table(data,'VariableNames',{'o1', 'o2',...
                                        'o3', 'o4'});
tbl.cue = cue;
tbl.id = id;
writetable(tbl,fullfile(output,'p_err.csv'));


%% Primary saccade latency (of ANY stimuli)

clearvars -except output stats1 stats2 qa1 qa2
data_quad = squeeze(mean(stats1.p_rt_ord,2,'omitnan'));
data_ord = squeeze(mean(stats2.p_rt_ord,2,'omitnan'));
data = [data_quad(:,1:4);data_ord(:,5:8)];

id = [stats1.subs_id';stats2.subs_id'];
cue = [zeros(size(stats1.subs_id))';ones(size(stats2.subs_id))']; % 0 - quadrant cue; 1 - order cue
final_sind = [stats1.final_sind;stats2.final_sind+length(stats1.subs_id)];

data = data(final_sind,:);
id = id(final_sind,:);
cue = cue(final_sind,:);

tbl = array2table(data,'VariableNames',{'o1', 'o2',...
                                        'o3', 'o4'});
tbl.cue = cue;
tbl.id = id;
writetable(tbl,fullfile(output,'p_rt.csv'));

%% Primary saccade latency (of target stimuli)

clearvars -except output stats1 stats2 qa1 qa2
data_quad = squeeze(mean(stats1.pt_rt_ord,2,'omitnan'));
data_ord = squeeze(mean(stats2.pt_rt_ord,2,'omitnan'));
data = [data_quad(:,1:4);data_ord(:,5:8)];

id = [stats1.subs_id';stats2.subs_id'];
cue = [zeros(size(stats1.subs_id))';ones(size(stats2.subs_id))']; % 0 - quadrant cue; 1 - order cue
final_sind = [stats1.final_sind;stats2.final_sind+length(stats1.subs_id)];

data = data(final_sind,:);
id = id(final_sind,:);
cue = cue(final_sind,:);

tbl = array2table(data,'VariableNames',{'o1', 'o2',...
                                        'o3', 'o4'});
tbl.cue = cue;
tbl.id = id;
writetable(tbl,fullfile(output,'pt_rt.csv'));

%% Primary saccade error (of target stimuli)
clearvars -except output stats1 stats2 qa1 qa2
data_quad = squeeze(mean(stats1.pt_err_ord,2,'omitnan'));
data_ord = squeeze(mean(stats2.pt_err_ord,2,'omitnan'));
data = [data_quad(:,1:4);data_ord(:,5:8)];

id = [stats1.subs_id';stats2.subs_id'];
cue = [zeros(size(stats1.subs_id))';ones(size(stats2.subs_id))']; % 0 - quadrant cue; 1 - order cue
final_sind = [stats1.final_sind;stats2.final_sind+length(stats1.subs_id)];

data = data(final_sind,:);
id = id(final_sind,:);
cue = cue(final_sind,:);

tbl = array2table(data,'VariableNames',{'o1', 'o2',...
                                        'o3', 'o4'});
tbl.cue = cue;
tbl.id = id;
writetable(tbl,fullfile(output,'pt_err.csv'));

%% Most accurate saccade error (of target stimuli)

clearvars -except output stats1 stats2 qa1 qa2
data_quad = squeeze(mean(stats1.a_err_ord,2,'omitnan'));
data_ord = squeeze(mean(stats2.a_err_ord,2,'omitnan'));
data = [data_quad(:,1:4);data_ord(:,5:8)];

id = [stats1.subs_id';stats2.subs_id'];
cue = [zeros(size(stats1.subs_id))';ones(size(stats2.subs_id))']; % 0 - quadrant cue; 1 - order cue
final_sind = [stats1.final_sind;stats2.final_sind+length(stats1.subs_id)];

data = data(final_sind,:);
id = id(final_sind,:);
cue = cue(final_sind,:);

tbl = array2table(data,'VariableNames',{'o1', 'o2',...
                                        'o3', 'o4'});
tbl.cue = cue;
tbl.id = id;
writetable(tbl,fullfile(output,'a_err.csv'));

%% Most accurate saccade latency (of target stimuli)

clearvars -except output stats1 stats2 qa1 qa2
data_quad = squeeze(mean(stats1.a_rt_ord,2,'omitnan'));
data_ord = squeeze(mean(stats2.a_rt_ord,2,'omitnan'));
data = [data_quad(:,1:4);data_ord(:,5:8)];

id = [stats1.subs_id';stats2.subs_id'];
cue = [zeros(size(stats1.subs_id))';ones(size(stats2.subs_id))']; % 0 - quadrant cue; 1 - order cue
final_sind = [stats1.final_sind;stats2.final_sind+length(stats1.subs_id)];

data = data(final_sind,:);
id = id(final_sind,:);
cue = cue(final_sind,:);

tbl = array2table(data,'VariableNames',{'o1', 'o2',...
                                        'o3', 'o4'});
tbl.cue = cue;
tbl.id = id;
writetable(tbl,fullfile(output,'a_rt.csv'));

%% Order error (of any stimuli)

clearvars -except output stats1 stats2 final_sind1 final_sind2
% transposition error % for all the not nan trials
final_sind1 = stats1.final_sind; final_sind2 = stats2.final_sind;
p_ord_err = cat(1,stats1.p_ord_err_ord(final_sind1,:,1:4),stats2.p_ord_err_ord(final_sind2,:,5:8)); % N participants x M trials x 8 conditions
% transposition error % for all the not nan trials
for ii = 1:numel(p_ord_err)
    curr = p_ord_err(ii);
    if isnan(curr)
    p_trans(ii) = nan;
    elseif curr ~= 0
        p_trans(ii) = 1; % one transposition error
    else
        p_trans(ii) = 0;
    end
end

p_trans = reshape(p_trans,size(p_ord_err,1),size(p_ord_err,2),size(p_ord_err,3));
data = squeeze(sum(p_trans,2,'omitnan')./64*100); %sum(~isnan(p_trans),2))

% Let's store that data
id = [stats1.subs_id';stats2.subs_id'];
cue = [zeros(size(stats1.subs_id))';ones(size(stats2.subs_id))'];
final_sind = [stats1.final_sind;stats2.final_sind+length(stats1.subs_id)];

id = id(final_sind,:);
cue = cue(final_sind,:);

tbl = array2table(data,'VariableNames',{'o1',...
                                        'o2',...
                                        'o3',...
                                        'o4'});
tbl.cue = cue;
tbl.id = id;
writetable(tbl,fullfile(output,'p_orderr.csv'));


%%%%%%%%%%%%%%% GOOD FOR NOW %%%%%%%%%%%%%%%
% %% Modeling - Mu
% clearvars -except output stats1 stats2 final_sind1 final_sind2
% data_quad = squeeze(stats1.model_sp_by_cue.mu(:,1,:));
% data_ord = squeeze(stats2.model_sp_by_cue.mu(:,2,:));
% 
% data = [data_quad;data_ord];
% 
% id = [stats1.subs_id';stats2.subs_id'];
% cue = [zeros(size(stats1.subs_id))';ones(size(stats2.subs_id))'];
% final_sind = [stats1.final_sind;stats2.final_sind+length(stats1.subs_id)];
% 
% data = data(final_sind,:);
% id = id(final_sind,:);
% cue = cue(final_sind,:);
% 
% tbl = array2table(data,'VariableNames',{'o1',...
%                                         'o2',...
%                                         'o3',...
%                                         'o4'});
% tbl.cue = cue;
% tbl.id = id;
% writetable(tbl,fullfile(output,'p_mu_wide.csv'));

%% Modeling - sd
% calculate how many dva per pixel
% Screen info
% distance from the eyes to the screen (in cm)
d = 67.7;
% width (in cm)
w = 38;
% height (in cm)
h = 30.5;
% x resolution (in pixel)
xres = 1280;
% y resolution (in pixel) of the monitor
yres = 1024;

% how many dva per pixel
xdva = atand(w/d)/xres;
ydva = atand(h/d)/yres;
dva_per_pix = (xdva+ydva)/2;

clearvars -except output stats1 stats2 final_sind1 final_sind2 dva_per_pix

data_quad = squeeze(stats1.model_swap.sd(:,1,:));
data_ord = squeeze(stats2.model_swap.sd(:,2,:));

data = [data_quad;data_ord]*dva_per_pix;

id = [stats1.subs_id';stats2.subs_id'];
cue = [zeros(size(stats1.subs_id))';ones(size(stats2.subs_id))'];
final_sind = [stats1.final_sind;stats2.final_sind+length(stats1.subs_id)];

data = data(final_sind,:);
id = id(final_sind,:);
cue = cue(final_sind,:);

tbl = array2table(data,'VariableNames',{'o1',...
                                        'o2',...
                                        'o3',...
                                        'o4'});
tbl.cue = cue;
tbl.id = id;
writetable(tbl,fullfile(output,'p_sd.csv'));

%% Modeling - guessing rate
clearvars -except output stats1 stats2 final_sind1 final_sind2
data_quad = squeeze(stats1.model_swap.gamma(:,1,:));
data_ord = squeeze(stats2.model_swap.gamma(:,2,:));

data = [data_quad;data_ord];

id = [stats1.subs_id';stats2.subs_id'];
cue = [zeros(size(stats1.subs_id))';ones(size(stats2.subs_id))'];
final_sind = [stats1.final_sind;stats2.final_sind+length(stats1.subs_id)];

data = data(final_sind,:);
id = id(final_sind,:);
cue = cue(final_sind,:);

tbl = array2table(data,'VariableNames',{'o1',...
                                        'o2',...
                                        'o3',...
                                        'o4'});
tbl.cue = cue;
tbl.id = id;
writetable(tbl,fullfile(output,'p_gamma.csv'));


%% Modeling - Swap rate
clearvars -except output stats1 stats2 final_sind1 final_sind2
data_quad = squeeze(stats1.model_swap.beta(:,1,:));
data_ord = squeeze(stats2.model_swap.beta(:,2,:));

data = [data_quad;data_ord];

id = [stats1.subs_id';stats2.subs_id'];
cue = [zeros(size(stats1.subs_id))';ones(size(stats2.subs_id))'];
final_sind = [stats1.final_sind;stats2.final_sind+length(stats1.subs_id)];

data = data(final_sind,:);
id = id(final_sind,:);
cue = cue(final_sind,:);

tbl = array2table(data,'VariableNames',{'o1',...
                                        'o2',...
                                        'o3',...
                                        'o4'});
tbl.cue = cue;
tbl.id = id;
writetable(tbl,fullfile(output,'p_beta.csv'));

%% Modeling - Target rate
clearvars -except output stats1 stats2 final_sind1 final_sind2
data_quad = squeeze(stats1.model_swap.alpha(:,1,:));
data_ord = squeeze(stats2.model_swap.alpha(:,2,:));

data = [data_quad;data_ord];

id = [stats1.subs_id';stats2.subs_id'];
cue = [zeros(size(stats1.subs_id))';ones(size(stats2.subs_id))'];
final_sind = [stats1.final_sind;stats2.final_sind+length(stats1.subs_id)];

data = data(final_sind,:);
id = id(final_sind,:);
cue = cue(final_sind,:);

tbl = array2table(data,'VariableNames',{'o1',...
                                        'o2',...
                                        'o3',...
                                        'o4'});
tbl.cue = cue;
tbl.id = id;
writetable(tbl,fullfile(output,'p_alpha.csv'));

%% How many trials used for modeling
clearvars -except output stats1 stats2 final_sind1 final_sind2 dva_per_pix
data_quad = nan(length(stats1.subs_id),4);
data_ord = nan(length(stats2.subs_id),4);
for ss = 1:length(stats1.subs_id)
    for pp = 1:4
clear temp
temp = stats1.model_swap.data(ss,1,pp); 
if ~isempty(temp{1})
temp = temp{1}.errors;
data_quad(ss,pp) = size(temp,2);
end
    end
end

for ss = 1:length(stats2.subs_id)
    for pp = 1:4
clear temp
temp = stats2.model_swap.data(ss,2,pp); 
if ~isempty(temp{1})
temp = temp{1}.errors;
data_ord(ss,pp) = size(temp,2);
end
    end
end

data = [data_quad;data_ord];

id = [stats1.subs_id';stats2.subs_id'];
cue = [zeros(size(stats1.subs_id))';ones(size(stats2.subs_id))'];
final_sind = [stats1.final_sind;stats2.final_sind+length(stats1.subs_id)];

data = data(final_sind,:);
id = id(final_sind,:);
cue = cue(final_sind,:);

tbl = array2table(data,'VariableNames',{'o1',...
                                        'o2',...
                                        'o3',...
                                        'o4'});
tbl.cue = cue;
tbl.id = id;
writetable(tbl,fullfile(output,'p_numtr.csv'));

%% Modeling - aic
clearvars -except output stats1 stats2 final_sind1 final_sind2
data_quad = squeeze(stats1.model_swap.AIC(:,1,:));
data_ord = squeeze(stats2.model_swap.AIC(:,2,:));

data = [data_quad;data_ord];

id = [stats1.subs_id';stats2.subs_id'];
cue = [zeros(size(stats1.subs_id))';ones(size(stats2.subs_id))'];
final_sind = [stats1.final_sind;stats2.final_sind+length(stats1.subs_id)];

data = data(final_sind,:);
id = id(final_sind,:);
cue = cue(final_sind,:);

tbl = array2table(data,'VariableNames',{'o1',...
                                        'o2',...
                                        'o3',...
                                        'o4'});
tbl.cue = cue;
tbl.id = id;
writetable(tbl,fullfile(output,'p_aic.csv'));

%% Modeling - bic
clearvars -except output stats1 stats2 final_sind1 final_sind2
data_quad = squeeze(stats1.model_swap.BIC(:,1,:));
data_ord = squeeze(stats2.model_swap.BIC(:,2,:));

data = [data_quad;data_ord];

id = [stats1.subs_id';stats2.subs_id'];
cue = [zeros(size(stats1.subs_id))';ones(size(stats2.subs_id))'];
final_sind = [stats1.final_sind;stats2.final_sind+length(stats1.subs_id)];

data = data(final_sind,:);
id = id(final_sind,:);
cue = cue(final_sind,:);

tbl = array2table(data,'VariableNames',{'o1',...
                                        'o2',...
                                        'o3',...
                                        'o4'});
tbl.cue = cue;
tbl.id = id;
writetable(tbl,fullfile(output,'p_bic.csv'));

%% Save mean errors and latencies for each participant
clearvars -except output stats1 stats2 final_sind1 final_sind2 dva_per_pix
p_err = readtable(fullfile(output,'p_err.csv'));
pt_err = readtable(fullfile(output,'pt_err.csv'));
a_err = readtable(fullfile(output,'a_err.csv'));
p_rt = readtable(fullfile(output,'p_rt.csv'));
pt_rt = readtable(fullfile(output,'pt_rt.csv'));
a_rt = readtable(fullfile(output,'a_rt.csv'));

tbl = table();
tbl.id = p_err.id;
tbl.cue = p_err.cue;
tbl.p_err = mean(p_err{:,1:4},2,'omitnan');
tbl.pt_err = mean(pt_err{:,1:4},2,'omitnan');
tbl.a_err = mean(a_err{:,1:4},2,'omitnan');
tbl.p_rt = mean(p_rt{:,1:4},2,'omitnan');
tbl.pt_rt = mean(pt_rt{:,1:4},2,'omitnan');
tbl.a_rt = mean(a_rt{:,1:4},2,'omitnan');

writetable(tbl,fullfile(output,'all_err_rt.csv'));