% Check data quality and determine trial participant exclusion in
% experiment 1
% By: Linjing Jiang
% Date: 5/1/2024

%% Let's first check data quality
clear
close all
clc

% output directory
output = './analysis/exp2';


%% Load data
cd(output)
%load('data1.mat') % v3_short, long delay
%load('data2.mat') % v3_short_delay, short delay
%qa1 = load('qa_exp1_longdelay.mat','qa'); qa1 = qa1.qa;
%qa2 = load('qa_exp1_shortdelay.mat','qa'); qa2 = qa2.qa;
load('stats_exp2_quad.mat')
load('stats_exp2_ord.mat')

% uncomment this if you run this script for the first time
final_sind1 = stats1.final_sind;
final_sind2 = stats2.final_sind;
final_sind = [final_sind1;final_sind2+length(stats1.subs_id)];

% number of participants in each cue condition
n_quad = length(stats1.subs_id);
n_ord = length(stats2.subs_id);
n_all = n_quad + n_ord;
n_tr = size(stats1.p_err,2);

%% Double check the trial identity
% which type of trial it is
final_sind = 1:60;
resp_type = [stats1.resp_type;stats2.resp_type];
no_data = sum(isnan(resp_type),2)/n_tr*100;
artifact = sum(resp_type(final_sind,:) == 0,2)/n_tr*100;
nobaseline = sum(resp_type(final_sind,:)==1,2)/n_tr*100;
miss_noresp = sum(resp_type(final_sind,:)==2,2)/n_tr*100;
miss_undefined = sum(resp_type(final_sind,:)==3,2)/n_tr*100;
early = sum(resp_type(final_sind,:)==4,2)/n_tr*100;
large_baseline = sum(resp_type(final_sind,:)==5,2)/n_tr*100;
tar = sum(resp_type(final_sind,:)==6,2)/n_tr*100;
nontar = sum(resp_type(final_sind,:)==7,2)/n_tr*100;
outside = sum(resp_type(final_sind,:)==8,2)/n_tr*100;
no_endfix = sum(resp_type(final_sind,:)==9,2)/n_tr*100;
others = sum(resp_type(final_sind,:)==10,2)/n_tr*100;
undefined_all = sum(resp_type(final_sind,:) == 3 | resp_type(final_sind,:)==10,2)/n_tr*100;


% exclude the following types of trials:
% 1. No data - isnan
% 2. Artifact during 500ms of the response window - 0
% 3. No baseline - 1
% 4. Miss trials - no response detected - 2, 3
% 5. Early/premature responses - 4
% 6. Primary saccade with no endpoint fixation - 9
% 7. Others - 10
exclude = no_data+artifact+nobaseline;
include = tar+nontar+outside;
figure;
histogram(include,'NumBins',10)
xlabel('Number of trials included in the final analysis')
ylabel('Number of participants')


% included trials per serial position
ord_all = [stats1.ord;stats2.ord];
for ii = 1:60 % for each participant
    ind1 = find(ord_all(ii,:)==1);
    ind2 = find(ord_all(ii,:)==2);
    ind3 = find(ord_all(ii,:)==3);
    ind4 = find(ord_all(ii,:)==4);
    included_sp(ii,:) = [sum(resp_type(ii,ind1)==6 | resp_type(ii,ind1)==7 | resp_type(ii,ind1)==8,2),...
        sum(resp_type(ii,ind2)==6 | resp_type(ii,ind2)==7 | resp_type(ii,ind2)==8,2),...
        sum(resp_type(ii,ind3)==6 | resp_type(ii,ind3)==7 | resp_type(ii,ind3)==8,2),...
        sum(resp_type(ii,ind4)==6 | resp_type(ii,ind4)==7 | resp_type(ii,ind4)==8,2)]/(n_tr/4)*100;
end

%numtr = 32; % previously used 24
perctr = 32/64*100;
sum(included_sp(:,1) >= perctr & included_sp(:,2) >= perctr & included_sp(:,3) >= perctr & included_sp(:,4) >= perctr)
final_sind = find(included_sp(:,1) >= perctr & included_sp(:,2) >= perctr & included_sp(:,3) >= perctr & included_sp(:,4) >= perctr);

% %sum(isoutlier(exclude))
% %sum(exclude >= 30) % if more than 40% of the trials are excluded
% sum(include >= 50)
% % final participant
% %final_sind = find(exclude < 30);
% %include_final = include(exclude < 50);
% final_sind = find(include >= 50);
% % for the stats1 (quad)
final_sind1 = final_sind(final_sind<=29);
% % for the stats2 (ord)
final_sind2 = final_sind(final_sind>29)-29;
% 
stats1.final_sind = final_sind1;
stats2.final_sind = final_sind2;

% store these ids
save('stats_exp2_quad.mat','stats1')
save('stats_exp2_ord.mat','stats2')


%%
% average percentage of each type of trials
mean(no_data)
std(no_data,0,1)

mean(artifact)
std(artifact,0,1)

mean(nobaseline)
std(nobaseline,0,1)

mean(miss_noresp)
std(miss_noresp,0,1)

mean(miss_undefined)
std(miss_undefined,0,1)

% miss
miss = miss_noresp+miss_undefined;
mean(miss)
std(miss,0,1)

mean(early)
std(early,0,1)

mean(large_baseline)
std(large_baseline,0,1)

mean(tar)
std(tar,0,1)

mean(nontar)
std(nontar,0,1)

mean(outside)
std(outside,0,1)

mean(no_endfix)
std(no_endfix,0,1)

mean(others)
std(others,0,1)

all_others = others+no_data;
mean(all_others)
std(all_others,0,1)

mean(tar+nontar+outside)
std(tar+nontar+outside,0,1)
min(tar+nontar+outside)
max(tar+nontar+outside)

% which participant with most miss undefined trials
clc
all_ids = [stats1.subs_id stats2.subs_id];
% 0 - artifact
[val,id] = max(artifact)
all_ids(id)
% 1 - no baseline
[val,id] = max(nobaseline)
all_ids(id)
% 2 - no response
[val,id] = max(miss_noresp)
all_ids(id)
% 3 - miss trial, undefined
[val,id] = max(miss_undefined)
all_ids(id)
% 4 - premature, early saccade
[val,id] = max(early)
all_ids(id)
% 5 - starting position more than 1.5 dva
[val,id] = max(large_baseline)
all_ids(id)
% 6 - target
[val,id] = max(tar)
all_ids(id)
% 7 - nontarget
[val,id] = max(nontar)
all_ids(id)
% 8 - outside 
[val,id] = max(outside)
all_ids(id)
% 9 - no endpoint fixation
[val,id] = max(no_endfix)
all_ids(id)
% 10 - others
[val,id] = max(others)
all_ids(id)

figure(3);clf
%colors = mat2cell(distinguishable_colors(8),ones(8,1),3);
bh = bar([no_data artifact nobaseline miss_noresp miss_undefined others no_endfix early large_baseline outside nontar tar],'Stacked')
set(bh,{'FaceColor'},repmat({'flat'},numel(bh),1))
colors =  mat2cell(jet(numel(bh)),ones(numel(bh),1), 3); 
set(bh,{'CData'},colors)
set(gca,'FontSize',14,'XTick',1:length(artifact),'XTickLabelRotation',-90)
xlabel('Participant')
ylabel('Percentage of trials excluded (%)')
legend({'no data','artifact','no baseline','no response','no saccade - undefined',...
    'early','large baseline','other undefined','outside of all ROIs',...
    'non-target','target'},'box','off','location','northeastoutside')
box off
saveas(gcf,'response_type.jpg')

%% Quadrant versus order cue trials - before exclusion
tar = [[sum(resp_type(1:29,:)==6,2);nan(2,1)] sum(resp_type(30:60,:)==6,2)]/256*100;
nontar = [[sum(resp_type(1:29,:)==7,2);nan(2,1)] sum(resp_type(30:60,:)==7,2)]/256*100;
outside = [[sum(resp_type(1:29,:)==8,2);nan(2,1)] sum(resp_type(30:60,:)==8,2)]/256*100;

include = tar+nontar+outside;

% which type of trials are included in each participant
[h,p,~,stat] = ttest(include(:,1),include(:,2));
stat
p

% plot
figure(1);clf
set(gcf, 'Position',  [0, 0, 800, 600])
boxWithDots3([include(:,1) include(:,2)],...
    {'Quadrant Cue','Order Cue'},0, ...
    [[0 0.4470 0.7410];[0.8500 0.3250 0.0980]],[p],{[1,2]})
ylabel('Number of trials included)')
box off
saveas(gcf,'rmtrial_by_cue_before.jpg')

%% Quadrant versus order cue trials - after exclusion
tar = [sum(resp_type(final_sind1,:)==6,2) [sum(resp_type(final_sind2+29,:)==6,2);nan(3,1)]]/256*100;
nontar = [sum(resp_type(final_sind1,:)==7,2) [sum(resp_type(final_sind2+29,:)==7,2);nan(3,1)]]/256*100;
outside = [sum(resp_type(final_sind1,:)==8,2) [sum(resp_type(final_sind2+29,:)==8,2);nan(3,1)]]/256*100;

include = tar+nontar+outside;
mean(include,'omitnan')
std(include,0,1,'omitnan')
min(include)
max(include)
% which type of trials are included in each participant
[h,p,~,stat] = ttest(include(:,1),include(:,2));
p
stat

% plot
figure(1);clf
set(gcf, 'Position',  [0, 0, 800, 600])
boxWithDots3([include(:,1) include(:,2)],...
    {'Quadrant Cue','Order Cue'},0, ...
    [[0 0.4470 0.7410];[0.8500 0.3250 0.0980]],[p],{[1,2]})
ylabel('Percentage of trials included (%)')
box off
saveas(gcf,'rmtrial_by_cue_after.jpg')

%% Check goodness of fit for models
bic_standard_quad = stats1.model_standard.BIC;
bic_standard_ord = stats2.model_standard.BIC;
bic_swap_quad = stats1.model_swap.BIC;
bic_swap_ord = stats2.model_swap.BIC;
bic_swap_sys_quad = stats1.model_swap_sys.BIC; 
bic_swap_sys_ord = stats2.model_swap_sys.BIC;

fprintf('bic standard model quadrant cue: %f, %f \n',...
    mean(bic_standard_quad(final_sind1,:,:),'all','omitnan'),...
    std(bic_standard_quad(final_sind1,:,:),0,'all','omitnan'))
fprintf('bic swap model quadrant cue: %f, %f \n',...
    mean(bic_swap_quad(final_sind1,:,:),'all','omitnan'),...
    std(bic_swap_quad(final_sind1,:,:),0,'all','omitnan'))
fprintf('bic swap + systematic bias model quadrant cue: %f, %f \n',...
    mean(bic_swap_sys_quad(final_sind1,:,:),'all','omitnan'),...
    std(bic_swap_sys_quad(final_sind1,:,:),0,'all','omitnan'))

fprintf('bic standard model order cue: %f, %f \n',...
    mean(bic_standard_ord(final_sind2,:,:),'all','omitnan'),...
    std(bic_standard_ord(final_sind2,:,:),0,'all','omitnan'))
fprintf('bic swap model order cue: %f, %f \n',...
    mean(bic_swap_ord(final_sind2,:,:),'all','omitnan'),...
    std(bic_swap_ord(final_sind2,:,:),0,'all','omitnan'))
fprintf('bic swap + systematic bias model order cue: %f, %f \n',...
    mean(bic_swap_sys_ord(final_sind2,:,:),'all','omitnan'),...
    std(bic_swap_sys_ord(final_sind2,:,:),0,'all','omitnan'))

figure(1);clf
set(gcf, 'Position',  [0, 0, 1100, 500])
tiledlayout(1,2)
nexttile % quadrant cue
m1 = squeeze(mean(bic_standard_quad(final_sind1,1,:),1,'omitnan')); % quadrant cue, mixture model
sd1 = squeeze(std(bic_standard_quad(final_sind1,1,:),0,1,'omitnan'))/sqrt(length(final_sind1));
m2 = squeeze(mean(bic_swap_quad(final_sind1,1,:),1,'omitnan')); % quadrant cue, swap model
sd2 = squeeze(std(bic_swap_quad(final_sind1,1,:),0,1,'omitnan'))/sqrt(length(final_sind1));
m3 = squeeze(mean(bic_swap_sys_quad(final_sind1,1,:),1,'omitnan')); % quadrant cue, swap + systematic bias model
sd3 = squeeze(std(bic_swap_sys_quad(final_sind1,1,:),0,1,'omitnan'))/sqrt(length(final_sind1));
set(gcf, 'Position',  [0, 0, 1100, 500])
errorbar(m1,sd1,'LineWidth',2,'LineStyle','-','Color',[0 0 0],'Marker','s','MarkerSize',1)
hold on
errorbar(m2,sd2,'LineWidth',2,'LineStyle','--','Color',[0.5 0.5 0.5],'Marker','s','MarkerSize',1)
hold on
errorbar(m3,sd3,'LineWidth',2,'LineStyle',':','Color',[0.8 0.8 0.8],'Marker','s','MarkerSize',1)
xlim([0 5])
ylim([1000 1400])
set(gca,'FontSize',20,'XTick',0:5,'XTickLabel',{'','1','2','3','4',''})
%legend('Standard model','Swap model','Swap model with systematic bias','box','off')
ylabel('BIC')
xlabel('Serial position')
title('Quadrant cue (N = 65)')
box off

nexttile % order cue
m1 = squeeze(mean(bic_standard_ord(final_sind2,2,:),1,'omitnan')); % order cue, mixture model
sd1 = squeeze(std(bic_standard_ord(final_sind2,2,:),0,1,'omitnan'))/sqrt(length(final_sind2));
m2 = squeeze(mean(bic_swap_ord(final_sind2,2,:),1,'omitnan')); % order cue, swap model
sd2 = squeeze(std(bic_swap_ord(final_sind2,2,:),0,1,'omitnan'))/sqrt(length(final_sind2));
m3 = squeeze(mean(bic_swap_sys_ord(final_sind2,2,:),1,'omitnan')); % order cue, swap + systematic bias model
sd3 = squeeze(std(bic_swap_sys_ord(final_sind2,2,:),0,1,'omitnan'))/sqrt(length(final_sind2));
set(gcf, 'Position',  [0, 0, 1100, 500])
errorbar(m1,sd1,'LineWidth',2,'LineStyle','-','Color',[0 0 0],'Marker','s','MarkerSize',1)
hold on
errorbar(m2,sd2,'LineWidth',2,'LineStyle','--','Color',[0.5 0.5 0.5],'Marker','s','MarkerSize',1)
hold on
errorbar(m3,sd3,'LineWidth',2,'LineStyle',':','Color',[0.8 0.8 0.8],'Marker','s','MarkerSize',1)
xlim([0 5])
ylim([1000 1400])
set(gca,'FontSize',20,'XTick',0:5,'XTickLabel',{'','1','2','3','4',''})
legend('Standard model','Swap model','Swap model with systematic bias',...
    'box','off','Location','northeastoutside')
ylabel('BIC')
xlabel('Serial position')
title('Order cue (N = 65)')
box off
saveas(gcf,'bic.jpg')

% Following Jeffrey-Raftery's (1995) guidelines, if the difference in BICs 
% between the two models is 0–2, this constitutes ‘weak’ evidence in favor 
% of the model with the smaller BIC; a difference in BICs between 2 and 6 
% constitutes ‘positive’ evidence; a difference in BICs between 6 and 10 
% constitutes ‘strong’ evidence; and a difference in BICs greater than 10 
% constitutes ‘very strong’ evidence in favor of the model with smaller BIC.
% A lower BIC value indicates a better fit

% Clearly, AIC does not depend directly on sample size. Moreover, generally 
% speaking, AIC presents the danger that it might overfit, whereas BIC 
% presents the danger that it might underfit, simply in virtue of how they 
% penalize free parameters (2*k in AIC; ln(N)*k in BIC). Diachronically, 
% as data is introduced and the scores are recalculated, at relatively 
% low N (7 and less) BIC is more tolerant of free parameters than AIC, but
% less tolerant at higher N (as the natural log of N overcomes 2).

% Additionally, AIC is aimed at finding the best approximating model to the
%unknown data generating process (via minimizing expected estimated K-L divergence). 
% As such, it fails to converge in probability to the true model (assuming 
% one is present in the group evaluated), whereas BIC does converge as N tends to infinity.
aic_standard = [stats1.model_standard.AIC ; stats2.model_standard.AIC];
aic_swap = [stats1.model_swap.AIC ; stats2.model_swap.AIC];
aic_swap_sys = [stats1.model_swap_sys.AIC ; stats2.model_swap_sys.AIC];
mean(aic_standard(final_sind,:,:),'all','omitnan')
mean(aic_swap(final_sind,:,:),'all','omitnan')
mean(aic_swap_sys(final_sind,:,:),'all','omitnan')
std(aic_standard(final_sind,:,:),0,'all','omitnan')
std(aic_swap(final_sind,:,:),0,'all','omitnan')
std(aic_swap_sys(final_sind,:,:),0,'all','omitnan')

% Both AIC and BIC indicates that the standard model fits the data best?
% We will fit the data using the pure swap model

%% Number of trials used for modeling
figure(2)
set(gcf, 'Position',  [0, 0, 550, 500])
tbl = readtable(fullfile(output,'p_numtr.csv'));
% temp = tbl{:,1:16};
% data = [temp(:,1:2:16);temp(:,2:2:16)];
data = tbl{:,1:8};

m = mean(data,'omitnan');
sd = std(data,0,1,'omitnan')/sqrt(size(data,1));

% effects of serial position and cue
% the figure has 8 columns, 1-4 serial positions for quadrant and order cue
% respectively
% scatterplot
xCenter = ones(size(data)).*[1:4,6:9];
xCenter = xCenter + (-1 + 2* rand(size(data)))*0.1;
xCenter(isnan(data)) = NaN;
x1 = reshape(xCenter,1,[]); y1 = reshape(data,1,[]);
%colors = repmat(distinguishable_colors(size(data,1),{'w','k'}),size(data,2),1);
colors = repmat([0 0 0],size(data,2)*size(data,1),1);
%scatter(x1,y1,'lineWidth',0.1,'markerFaceColor',[0.8 0.8 0.8],'markerEdgeColor',[0.2 0.2 0.2])
% for ii = 1:length(x1)
% scatter(x1(ii),y1(ii),'lineWidth',0.1,'MarkerFaceColor',colors(ii,:),'markerEdgeColor',colors(ii,:),'SizeData',5)
% hold on
% end
errorbar(1:4,m(1:4),sd(1:4),'LineWidth',2,'LineStyle','--','Color',[0.5 0.5 0.5],'Marker','s','MarkerSize',10) % [27,158,119]./255
hold on
errorbar(1:4,m(5:8),sd(5:8),'LineWidth',2,'Color',[0 0 0],'Marker','o','MarkerSize',10) % [217,95,2]./255
xlim([0 5])
%ylim([550 860])
set(gca,'FontSize',fz,'XTick',0:5,'XTickLabel',{'','1','2','3','4',''})
ylabel('Number of trials')
xlabel('Serial position')
box off
saveas(gcf,'numtr_swap.jpg')



%% 2. Plot eye movement path
% Plot superimposed gaze path

% First let's load the file directories
% data directory
data_dir1 = '/mnt/Data1/eye_data/vswmda_eye/v3_short/';
file_dirs1 = [dir(fullfile(data_dir1,'sMGS*','p*'))];

% get subjects and run information
clear subs1 runs1 subs2 runs2
for ii = 1:length(file_dirs1)
    subs1{ii} = file_dirs1(ii).name;
    runs1(ii) = str2num(erase(file_dirs1(ii).folder,[data_dir1 'sMGS']));
end
subs_id1 = unique(subs1); % subject id
runs_id1 = unique(runs1); % run id

% Initiate some variables
x_trial_all_mat = []; y_trial_all_mat = []; tarx = []; tary = [];
x_trial_all_flip_mat = []; y_trial_all_flip_mat = [];

% which participant do you want to plot?
ii = final_sind1(10); % the fourth participant in the final list
ind = find(contains(subs1,subs_id1{ii})); % find all the runs for that participant
% run id
run_names = str2double(cellfun(@(x) erase(x,[data_dir1,'sMGS']),{file_dirs1(ind).folder},'UniformOutput',false));

for jj = 1:8 % for each run
        
        [isrun,idx_run] = ismember(jj,run_names);
        if ~isrun % if there is no such run
            continue;
        end
        
        % if there is such run
        % get the file name and directory of preprocessed step6 data for
        % each participant and each run
        file = dir(fullfile(file_dirs1(ind(idx_run)).folder,file_dirs1(ind(idx_run)).name,'result','*step6.mat'));

        if isempty(file) % if there is no preprocessed data *step6.mat
            continue;
        end
        
        % load the preprocessed data
        load(fullfile(file.folder,file.name),'edf');

% 1, start
% 2, fixation
% 3, 5, 7, 9, targets
% 10, delay
% 11, cue
% 13, iti
% 14, end
eye = mod(contains(edf.record.eye,'left')+1,2)+1;
%eye = contains(edf.record.eye,'left');
x = edf.samples.x_deg_clean_dc(:,eye); 
y = edf.samples.y_deg_clean_dc(:,eye);

% event start
ev_srt = edf.events.msg.ind_srt(:,[1 2 3 5 7 9 10 11 13 14]); 

% duration of each task period
ev_len = diff(ev_srt,1,2);

% mainly focus on the delay period
delay_len = round(ev_len(:,7)/1000)*1000;

% which trials need to fill in the na values
delay_len_add = max(delay_len) - delay_len;

% we want to start a trial 250ms before 
ev_srt_new = ev_srt;
ev_srt_new(:,1) = [];
ev_srt_new(:,1) = ev_srt_new(:,2)-250;
ev_srt_new(:,end) = [];

% plotted event timing
ev_plot = [250 1000 1000 1000 500 10000 1500];
ev_plot = cumsum(ev_plot);

% max trial length
tr_len = ev_srt_new(:,end)-ev_srt_new(:,1)+1+delay_len_add;

% get the samples trial by trial
% for ii = 1:16
% x_trial{ii} = [x(ev_srt_new(ii,1):ev_srt_new(ii,7));nan(delay_len_add(ii),1);x((ev_srt_new(ii,7)+1):ev_srt_new(ii,end))];
% y_trial{ii} = [y(ev_srt_new(ii,1):ev_srt_new(ii,7));nan(delay_len_add(ii),1);y((ev_srt_new(ii,7)+1):ev_srt_new(ii,end))];
% tr_len(ii) = length(x_trial{ii});
% end
clear x_trial y_trial
for kk = 1:16
x_trial(:,kk) = [x(ev_srt_new(kk,1):ev_srt_new(kk,7));nan(delay_len_add(kk),1);x((ev_srt_new(kk,7)+1):ev_srt_new(kk,end));nan(max(tr_len)-tr_len(kk),1)];
y_trial(:,kk) = [y(ev_srt_new(kk,1):ev_srt_new(kk,7));nan(delay_len_add(kk),1);y((ev_srt_new(kk,7)+1):ev_srt_new(kk,end));nan(max(tr_len)-tr_len(kk),1)];
end

tr_nodc_err = unique([edf.trackloss.nanerr;edf.trackloss.nodc_trial']);
final_tr =  1:16; final_tr(tr_nodc_err) = []; % valid trials
x_trial_all{jj,1} = x_trial(:,final_tr); % only select subset of trials
y_trial_all{jj,1} = y_trial(:,final_tr);
run_len(jj) = size(x_trial,1);

tarx = [tarx,edf.param.tarx_deg(final_tr)'];
tary = [tary,edf.param.tary_deg(final_tr)'];

xsign = sign(edf.param.tarx_deg(final_tr));
ysign = sign(edf.param.tary_deg(final_tr));

% x will be positive; y will be all negative
x_trial_all_flip{jj,1} = x_trial_all{jj,1}.*xsign';
y_trial_all_flip{jj,1} = y_trial_all{jj,1}.*(-ysign)';

final_tr_all{jj} = final_tr;
end

for jj = 1:8 
% tr_nodc_err = unique([edf.trackloss.nanerr;edf.trackloss.nodc_trial']);
% final_tr =  1:16; final_tr(tr_nodc_err) = []; % valid trials
% 
final_tr = final_tr_all{jj};
x_trial_all_mat = [x_trial_all_mat,[x_trial_all{jj,1};nan(max(run_len)-run_len(jj),length(final_tr))]];
y_trial_all_mat = [y_trial_all_mat,[y_trial_all{jj,1};nan(max(run_len)-run_len(jj),length(final_tr))]];
x_trial_all_flip_mat = [x_trial_all_flip_mat,[x_trial_all_flip{jj,1};nan(max(run_len)-run_len(jj),length(final_tr))]];
y_trial_all_flip_mat = [y_trial_all_flip_mat,[y_trial_all_flip{jj,1};nan(max(run_len)-run_len(jj),length(final_tr))]];

end

tm = (1:max(run_len))*2/1000;
f = figure('Renderer', 'painters', 'Position', [10 10 2000 600])
for ii = 1:size(x_trial_all_flip_mat,2)
plot(tm,x_trial_all_flip_mat(:,ii),'Color',[0 0.4470 0.7410],'LineWidth',1)
hold on
plot(tm,y_trial_all_flip_mat(:,ii),'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
hold on
end

xlim([0 16])
ylim([-11 11])

for ii = 1:length(ev_plot)-1
line([ev_plot(ii) ev_plot(ii)]/1000,[-11 11],'Color','k','LineStyle','--','LineWidth',1.5)
hold on
%text(ev_plot(ii)/1000+0.1,8,task_string{ii})
end

xlabel('Time (s)')
ylabel('Gaze Position (dva)')
box off
legend({'X gaze position (dva)', 'Y gaze position (dva)'},'box','off','location','north')
box off
%title(subs_id1(ii))
set(gca,'FontSize',20)
saveas(f,fullfile(output,['gaze_path_' subs_id1{final_sind1(4)} '.jpg']))
%saveas(f,fullfile(output,['gaze_path_' subs_id1{final_sind1(4)} '.eps']))
print(f, '-depsc2', '-r300', '-tiff', '-loose', fullfile(output,['gaze_path_' subs_id1{final_sind1(4)} '.eps']))