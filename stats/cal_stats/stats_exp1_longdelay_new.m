% Statistical analysis for experiment 1, long delay (v3_short)

% Exp1 consists of v3_short (long delay) and v3_short_delay (short delay)
% Aims to examine serial-order working memory, with ordinal rank, testing
% cue and eccentricity as within-subject factors and delay as
% between-subject factors

% This script does the following:
% 1. Estimate quality of the data and plot the number of trials excluded
% 2. Determine which participant to exclude
% 3. Calculate errors and summarize across trials and participants

% By: Linjing Jiang
% Verdate: 04/13/2024

%% Load directories
clear
close all
clc

% output directory
output = './analysis/exp1';

% add scripts to the path
addpath(genpath('./scripts/Cloned-Repo/'))

% data directory
data_dir1 = './eye_data/v3_short/';
file_dirs1 = [dir(fullfile(data_dir1,'sMGS*','p*'))];

% for the following analysis, all the '1' will refer to long delay
% (v3_short version) and all the '2' variables refer to short delay
% (v3_short_delay)

% get subjects and run information
clear subs1 runs1 subs2 runs2
for ii = 1:length(file_dirs1)
    subs1{ii} = file_dirs1(ii).name;
    runs1(ii) = str2num(erase(file_dirs1(ii).folder,[data_dir1 'sMGS']));
end
subs_id1 = unique(subs1); % subject id
runs_id1 = unique(runs1); % run id

%% Load the preprocessed edf data from each participant and perform further analysis
% for the long delay, v3_short

% Initiate some variables that stores data quality information for all the
% runs and participants
perc = nan(length(subs_id1),length(runs_id1));
nblink = nan(length(subs_id1),length(runs_id1));
ntr_breakfix = nan(length(subs_id1),length(runs_id1));
ntr_err = nan(length(subs_id1),length(runs_id1));
ntr_rm = nan(length(subs_id1),length(runs_id1));
ntr_nodc = nan(length(subs_id1),length(runs_id1));
ntr_err_nodc = nan(length(subs_id1),length(runs_id1));
tr_breakfix = cell(length(subs_id1),length(runs_id1));
tr_err = cell(length(subs_id1),length(runs_id1));
tr_rm = cell(length(subs_id1),length(runs_id1));
tr_nodc = cell(length(subs_id1),length(runs_id1));

% Loop through each participant and run
for ii = 1:length(subs_id1) % for each participant
    ind = find(contains(subs1,subs_id1{ii})); % find all the runs for that participant
    % run id
run_names = str2double(cellfun(@(x) erase(x,[data_dir1,'sMGS']),{file_dirs1(ind).folder},'UniformOutput',false));


                % Let's first initialize some variables
                % supposingly these variables will store error and other
                % data across trials and runs
                % they are all M x 1 vector, where M = number of trials
                % x number of runs, e.g., 16 trials x 8 runs = 128 rows in
                % our current experiment
            p_resp = nan(128,1); % response type of that trial
            p_ord = nan(128,1); % serial position of the target
            p_quad = nan(128,1); % quadrant of the target
            p_img = nan(128,1); % image type of the target
            p_task = nan(128,1); % task type (cue type)
            p_ecc = nan(128,1); % eccentricity of the target
            p_xtar = nan(128,1); p_ytar = nan(128,1); % x and y location of the target in dva
            p_xtar_pix = nan(128,1); p_ytar_pix = nan(128,1); % x and y location of the target in pixel
            p_xntar_pix = nan(128,3); p_yntar_pix = nan(128,3); % x and y location of the three non-targets in pixel
            p_rm = nan(128,1); % trials that DO NOT need to be removed
            

            p_err = nan(128,1); % primary saccade error (that enters AOI of any stimulus)
            a_err = nan(128,1); % closest saccade error
            a_rt = nan(128,1); % closest saccade latency
 
            pt_err = nan(128,1); % primary saccade error (that enters target AOI)
            p_rt = nan(128,1); % primary saccade latency (that enters AOI of any stimulus)
            pt_rt = nan(128,1); % primary saccade latency (that enters target AOI)

            p_xend = nan(128,1); % primary saccade endpoint x location in dva
            p_yend = nan(128,1); % primary saccade endpoint y location in dva
            p_xend_pix = nan(128,1); p_yend_pix = nan(128,1); % primary saccade endpoint locations in pixel
            % primary saccade angular error, quadrant error, and serial
            % position (transposition) error (saccade minus target)
            p_theta_err = nan(128,1); p_quad_err = nan(128,1); p_ord_err = nan(128,1); 
            % the most accurate saccade angular error, quadrant error, and transposition error (saccade minus target) 
            a_theta_err = nan(128,1); a_quad_err = nan(128,1); a_ord_err = nan(128,1);

    for jj = 1:length(ind) % for each run
                % if there are no such run
        [isrun,idx_run] = ismember(jj,run_names);


        if isrun
            
%         % skip 114 - s12, sMGS4
%         if ind(jj) == 114
%             continue
%         end
        
        % get the file name and directory of preprocessed step6 data for
        % each participant and each run
        file = dir(fullfile(file_dirs1(ind(jj)).folder,file_dirs1(ind(jj)).name,'result','*step6.mat'));
        %   fprintf('\n');
        %   fprintf('%s subject %s run %s starts...',subs_id1{ii},erase(file_dirs(ff).folder,[data_dir,'/']));

        if ~isempty(file) % if there is preprocessed data *step6.mat
            % load the preprocessed data
            load(fullfile(file.folder,file.name),'edf','set');

            %% 1. let's calculate and store the quality of the data for that run and participant

            % what is the percentage of samples lost during that run?
            perc(ii,jj) = edf.trackloss.perc;
            % what is the number of blinks during that run?
            nblink(ii,jj) = edf.blink.num;

            % let's store the number of trials breaking fixation, does not have a primary
            % saccade, without drift correction, and removed for participant ii and
            % run jj
            ntr_breakfix(ii,jj) = size(edf.trackloss.breakfix,1);
            ntr_err(ii,jj) = size(edf.trackloss.nanerr,1);
            ntr_rm(ii,jj) = size(edf.trackloss.rmtrial,1);
            ntr_nodc(ii,jj) = length(edf.trackloss.nodc_trial);
            ntr_err_nodc(ii,jj) = ntr_nodc(ii,jj) +  ntr_err(ii,jj); % THIS IS WRONG!!!

            % the actual trial number that breaks fixation, without primary
            % saccade, without drift correction and removed
            tr_breakfix{ii,jj} = edf.trackloss.breakfix;
            tr_err{ii,jj} = edf.trackloss.nanerr;
            tr_rm{ii,jj} = edf.trackloss.rmtrial;
            tr_nodc{ii,jj} = edf.trackloss.nodc_trial;

            %% 2. Let's calculate and summarize er16*(jj-1)+[1:16]rors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Step 1: concatenate errors and measures across trials and
            % runs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            

                % trials that need to be removed (let's only remove trials
                % with big error and without drift correction
                temp1 = ones(16,1);
                temp2 = [edf.trackloss.nanerr;edf.trackloss.nodc_trial'];
                temp1(temp2) = 0;
                p_rm(16*(jj-1)+[1:16]) = temp1;
     
                % response type
                p_resp(16*(jj-1)+[1:16]) = edf.cal.resp_type; 
                
                % cue type (task) - quadrant cue is 0, order cue is 1
                p_task(16*(jj-1)+[1:16]) = floor((jj-1)/4)*ones(16,1);
                % % for experiment 2
                % p_task(16*(jj-1)+[1:16]) = zeros(16,1); % for model_quad
                % p_task(16*(jj-1)+[1:16]) = ones(16,1); % for model_ord

                % target's serial position
                p_ord(16*(jj-1)+[1:16]) = edf.param.probe_order(1:16);

                % target quadrant
                p_quad(16*(jj-1)+[1:16]) = edf.param.probe_quad(1:16);

                % target image
                p_img(16*(jj-1)+[1:16]) = edf.param.probe_img(1:16);

                % target eccentricity
                p_ecc(16*(jj-1)+[1:16]) = edf.param.ecc(1:16);

                % Primary saccade for any stimuuli
                p_err(16*(jj-1)+[1:16]) = edf.cal.primary_sac_err;

                % Primary saccade for target stimulus
                pt_err(16*(jj-1)+[1:16]) = edf.cal.primary_sac_tar_err;

                % Most accurate saccade
                a_err(16*(jj-1)+[1:16]) = edf.cal.acc_sac_tar_err;
                % % Number of saccades
                % n_sac = [n_sac;edf.cal.num_sac(:,1)];
                
                % most accurate saccade latency
                a_rt(16*(jj-1)+[1:16]) = edf.cal.acc_sac_tar_rt;
                
                % latency primary saccade for any stimulus
                p_rt(16*(jj-1)+[1:16]) = edf.cal.primary_sac_rt;

                % latency: primary saccade for target stimulus
                pt_rt(16*(jj-1)+[1:16]) = edf.cal.primary_sac_tar_rt;

                % other type of errors
                % angular error, quadrant disposition error, order disposition error
                p_theta_err(16*(jj-1)+[1:16]) = edf.cal.p_minus_tar_theta;
                p_quad_err(16*(jj-1)+[1:16]) = edf.cal.p_minus_tar_quad;
                p_ord_err(16*(jj-1)+[1:16]) = edf.cal.p_minus_tar_ord;

                a_theta_err(16*(jj-1)+[1:16]) = edf.cal.a_minus_tar_theta;
                a_quad_err(16*(jj-1)+[1:16]) = edf.cal.a_minus_tar_quad;
                a_ord_err(16*(jj-1)+[1:16]) = edf.cal.a_minus_tar_ord;

                % to calculate systematic and unsystematic error, we want
                % to make sure to store the endpoint fixation (averaged)
                % position for each saccade
                sacind = edf.cal.primary_sac_ind; % this is the index relative to sac_dc
                for dd = 1:16
                    if ~isnan(sacind(dd))
                        % endpoint fixation position, not saccade endpoint
                        % position
                        xend_pix(dd) = edf.events.sac_dc.endfix_avg_x_pix(sacind(dd));
                        yend_pix(dd) = edf.events.sac_dc.endfix_avg_y_pix(sacind(dd));
                        % saccade end position in dva
                        xend(dd) = edf.events.sac_dc.endfix_avg_x(sacind(dd));
                        yend(dd) = edf.events.sac_dc.endfix_avg_y(sacind(dd));
                        % stimulus location
                        xtar(dd) = edf.param.tarx_deg((dd));
                        ytar(dd) = edf.param.tary_deg((dd));
                        xtar_pix(dd) = edf.param.tarx((dd));
                        ytar_pix(dd) = edf.param.tary((dd));

                        % distractor location
                        % which items are distractors
                        ntar_idx = 1:4; ntar_idx(edf.param.probe_order(dd)) = [];
                        xntar_pix(dd,:) = edf.param.stimx(dd,ntar_idx); % nontarget x location in pixel
                        yntar_pix(dd,:) = edf.param.stimy(dd,ntar_idx); % nontarget y location in pixel
                    else
                        xend_pix(dd) = nan; yend_pix(dd) = nan;
                        xend(dd) = nan; yend(dd) = nan;
                        xtar(dd) = nan; ytar(dd) = nan;
                        xtar_pix(dd) = nan; ytar_pix(dd) = nan;
                        xntar_pix(dd,1:3) = nan(1,3); yntar_pix(dd,1:3) = nan(1,3);
                    end
                end

                p_xend(16*(jj-1)+[1:16]) = xend;
                p_yend(16*(jj-1)+[1:16]) = yend;
                p_xend_pix(16*(jj-1)+[1:16]) = xend_pix;
                p_yend_pix(16*(jj-1)+[1:16]) = yend_pix;
                p_xtar(16*(jj-1)+[1:16]) = xtar;
                p_ytar(16*(jj-1)+[1:16]) = ytar;
                p_xtar_pix(16*(jj-1)+[1:16]) = xtar_pix;
                p_ytar_pix(16*(jj-1)+[1:16]) = ytar_pix;
                p_xntar_pix(16*(jj-1)+[1:16],:) = xntar_pix;
                p_yntar_pix(16*(jj-1)+[1:16],:) = yntar_pix;    

%         else % if there is no preprocessed data
%             perc(ii,jj) = nan;
%             nblink(ii,jj) = nan;
%             ntr_breakfix(ii,jj) = nan;
%             ntr_err(ii,jj) = nan;
%             ntr_rm(ii,jj) = nan;
%             ntr_nodc(ii,jj) = nan;
%             ntr_err_nodc(ii,jj) = nan;
%             tr_breakfix{ii,jj} = [];
%             tr_err{ii,jj} = [];
%             tr_rm{ii,jj} = [];
%             tr_nodc{ii,jj} = [];
        end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Step 2: organize errors and measures by serial position,
            % quadrant and image type across runs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% we first want to store some of the key variables as calculated above
% across runs
stats1.resp_type(ii,:) = p_resp; % response type
stats1.p_err(ii,:) = p_err; % participants x trial
stats1.pt_err(ii,:) = pt_err;
stats1.a_err(ii,:) = a_err;
stats1.a_rt(ii,:) = a_rt;
stats1.p_rt(ii,:) = p_rt;
stats1.pt_rt(ii,:) = pt_rt;
stats1.ord(ii,:) = p_ord;
stats1.quad(ii,:) = p_quad;
stats1.kept_trial(ii,:) = p_rm;
stats1.img(ii,:) = p_img;
stats1.cue(ii,:) = p_task;
stats1.ecc(ii,:) = p_ecc;
stats1.p_xend_pix(ii,:) = p_xend_pix;
stats1.p_yend_pix(ii,:) = p_yend_pix;
% stats1.p_xend_deg(ii,:) = p_xend;
% stats1.p_yend_deg(ii,:) = p_yend;
stats1.xtar_pix(ii,:) = p_xtar_pix;
stats1.ytar_pix(ii,:) = p_ytar_pix;
% stats1.xtar_deg(ii,:) = p_xtar;
% stats1.ytar_deg(ii,:) = p_ytar;
stats1.xntar_pix(ii,:,:) = p_xntar_pix; % participants x trial x number of non-targets 
stats1.yntar_pix(ii,:,:) = p_yntar_pix;

%%%%%%%%%%%% Organized by Serial Position %%%%%%%%%%%%%%%%%%%%%%%%%

% preallocate
for kk = 1:4 % for each ordinal position
    for hh = 1:2 % for each task (cue type)
        % find trials where serial position is kk, cue type is hh-1, and do not need to be revmoed
    ind = find(p_ord==kk & p_task==hh-1 & p_rm==1); 
    % primary saccade for any stimulus, error to target, N participant x 8 conditions
    % (serial order 1-4 for quadarnt cue, 1-4 for order cue)
p_err_new(:,kk+(hh-1)*4) = [p_err(ind);nan(16-length(ind),1)]; % at most 16 trials per serial position

% most accurate saccade error to target
a_err_new(:,kk+(hh-1)*4) = [a_err(ind);nan(16-length(ind),1)];
%n_sac_new(:,kk+(hh-1)*4) = [n_sac(ind);nan(16-length(ind),1)];

% closest saccade latency
a_rt_new(:,kk+(hh-1)*4) = [a_rt(ind);nan(16-length(ind),1)];

% primary saccade for any stimulus, latency
p_rt_new(:,kk+(hh-1)*4) = [p_rt(ind);nan(16-length(ind),1)];
%n_tr(1,kk+(hh-1)*4) = sum(~isnan(p_err(ind)));

% primary saccade for target, error to target and latency
pt_err_new(:,kk+(hh-1)*4) = [pt_err(ind);nan(16-length(ind),1)];
pt_rt_new(:,kk+(hh-1)*4) = [pt_rt(ind);nan(16-length(ind),1)];

% primary saccade endpoint x and y location in dva and pixel
p_xend_new(:,kk+(hh-1)*4) = [p_xend(ind);nan(16-length(ind),1)];
p_yend_new(:,kk+(hh-1)*4) = [p_yend(ind);nan(16-length(ind),1)];
p_xend_pix_new(:,kk+(hh-1)*4) = [p_xend_pix(ind);nan(16-length(ind),1)];
p_yend_pix_new(:,kk+(hh-1)*4) = [p_yend_pix(ind);nan(16-length(ind),1)];

% target x and y location in dva and pixel
p_xtar_new(:,kk+(hh-1)*4) = [p_xtar(ind);nan(16-length(ind),1)];
p_ytar_new(:,kk+(hh-1)*4) = [p_ytar(ind);nan(16-length(ind),1)];
p_xtar_pix_new(:,kk+(hh-1)*4) = [p_xtar_pix(ind);nan(16-length(ind),1)];
p_ytar_pix_new(:,kk+(hh-1)*4) = [p_ytar_pix(ind);nan(16-length(ind),1)];

% primary saccade (for any stimulus) angular, quadrant, and serial position error
p_theta_err_new(:,kk+(hh-1)*4) = [p_theta_err(ind);nan(16-length(ind),1)];
p_quad_err_new(:,kk+(hh-1)*4) = [p_quad_err(ind);nan(16-length(ind),1)];
p_ord_err_new(:,kk+(hh-1)*4) = [p_ord_err(ind);nan(16-length(ind),1)];

% most accurate saccade angular, quadrant, and serial position error
a_theta_err_new(:,kk+(hh-1)*4) = [a_theta_err(ind);nan(16-length(ind),1)];
a_quad_err_new(:,kk+(hh-1)*4) = [a_quad_err(ind);nan(16-length(ind),1)];
a_ord_err_new(:,kk+(hh-1)*4) = [a_ord_err(ind);nan(16-length(ind),1)];

% error by eccentricity
% trials with small eccentricity and large eccentricity (ind1 and ind2
% respectively)
ind1 = find(p_ord==kk & p_task==hh-1 & p_ecc==3 & p_rm == 1);
ind2 = find(p_ord==kk & p_task==hh-1 & p_ecc==5.5 & p_rm == 1);

% primary saccade (for any stimulus), error and latency, small and large eccentricity
p_err_secc(:,kk+(hh-1)*4) = [p_err(ind1);nan(8-length(ind1),1)]; % at most 8 trials per serial position per eccentricity per cue type
p_err_lecc(:,kk+(hh-1)*4) = [p_err(ind2);nan(8-length(ind2),1)];
p_rt_secc(:,kk+(hh-1)*4) = [p_rt(ind1);nan(8-length(ind1),1)];
p_rt_lecc(:,kk+(hh-1)*4) = [p_rt(ind2);nan(8-length(ind2),1)];

% primary saccade (for target), error and latency, small and large eccentricity
pt_err_secc(:,kk+(hh-1)*4) = [pt_err(ind1);nan(8-length(ind1),1)];
pt_err_lecc(:,kk+(hh-1)*4) = [pt_err(ind2);nan(8-length(ind2),1)];
pt_rt_secc(:,kk+(hh-1)*4) = [pt_rt(ind1);nan(8-length(ind1),1)];
pt_rt_lecc(:,kk+(hh-1)*4) = [pt_rt(ind2);nan(8-length(ind2),1)];

% most accurate saccade, error, small and large eccentricity
a_err_secc(:,kk+(hh-1)*4) = [a_err(ind1);nan(8-length(ind1),1)];
a_err_lecc(:,kk+(hh-1)*4) = [a_err(ind2);nan(8-length(ind2),1)];

a_rt_secc(:,kk+(hh-1)*4) = [a_rt(ind1);nan(8-length(ind1),1)];
a_rt_lecc(:,kk+(hh-1)*4) = [a_rt(ind2);nan(8-length(ind2),1)];

% primary saccade (for any stimulus), endpoint x and y location in pixel, small and
% large eccentricity
p_xend_pix_new_secc(:,kk+(hh-1)*4) = [p_xend_pix(ind1);nan(8-length(ind1),1)];
p_yend_pix_new_secc(:,kk+(hh-1)*4) = [p_yend_pix(ind1);nan(8-length(ind1),1)];
p_xend_pix_new_lecc(:,kk+(hh-1)*4) = [p_xend_pix(ind2);nan(8-length(ind2),1)];
p_yend_pix_new_lecc(:,kk+(hh-1)*4) = [p_yend_pix(ind2);nan(8-length(ind2),1)];

% target x and y location in pixel, small and large eccentricity
p_xtar_pix_new_secc(:,kk+(hh-1)*4) = [p_xtar_pix(ind1);nan(8-length(ind1),1)];
p_ytar_pix_new_secc(:,kk+(hh-1)*4) = [p_ytar_pix(ind1);nan(8-length(ind1),1)];
p_xtar_pix_new_lecc(:,kk+(hh-1)*4) = [p_xtar_pix(ind2);nan(8-length(ind2),1)];
p_ytar_pix_new_lecc(:,kk+(hh-1)*4) = [p_ytar_pix(ind2);nan(8-length(ind2),1)];
    end
end

% Store these variables in 'stats'
% Number of participant x Number of trials x Number of experimental
% conditions
% Here we have 8 conditions,
% serial order 1-4 for quadrant cue, and 1-4 for order cue
stats1.p_err_ord(ii,:,:) = p_err_new; 
stats1.a_err_ord(ii,:,:) = a_err_new;
stats1.a_rt_ord(ii,:,:) = a_rt_new;
%stats1.n_sac_ord(ii,:,:) = n_sac_new;
stats1.p_rt_ord(ii,:,:) = p_rt_new;
stats1.pt_err_ord(ii,:,:) = pt_err_new;
stats1.pt_rt_ord(ii,:,:) = pt_rt_new;
%stats1.n_tr_ord(ii,:,:) = n_tr;

stats1.p_xend_ord(ii,:,:) = p_xend_new;
stats1.p_yend_ord(ii,:,:) = p_yend_new;
stats1.p_xtar_ord(ii,:,:) = p_xtar_new;
stats1.p_ytar_ord(ii,:,:) = p_ytar_new;

stats1.p_xend_pix_ord(ii,:,:) = p_xend_pix_new;
stats1.p_yend_pix_ord(ii,:,:) = p_yend_pix_new;
stats1.p_xtar_pix_ord(ii,:,:) = p_xtar_pix_new;
stats1.p_ytar_pix_ord(ii,:,:) = p_ytar_pix_new;

stats1.p_xend_pix_ord_secc(ii,:,:) = p_xend_pix_new_secc;
stats1.p_xend_pix_ord_lecc(ii,:,:) = p_xend_pix_new_lecc;
stats1.p_yend_pix_ord_secc(ii,:,:) = p_yend_pix_new_secc;
stats1.p_yend_pix_ord_lecc(ii,:,:) = p_yend_pix_new_lecc;
stats1.p_xtar_pix_ord_secc(ii,:,:) = p_xtar_pix_new_secc;
stats1.p_xtar_pix_ord_lecc(ii,:,:) = p_xtar_pix_new_lecc;
stats1.p_ytar_pix_ord_secc(ii,:,:) = p_ytar_pix_new_secc;
stats1.p_ytar_pix_ord_lecc(ii,:,:) = p_ytar_pix_new_lecc;

stats1.p_err_ord_secc(ii,:,:) = p_err_secc;
stats1.p_err_ord_lecc(ii,:,:) = p_err_lecc;
stats1.p_rt_ord_secc(ii,:,:) = p_rt_secc;
stats1.p_rt_ord_lecc(ii,:,:) = p_rt_lecc;

stats1.pt_err_ord_secc(ii,:,:) = pt_err_secc;
stats1.pt_err_ord_lecc(ii,:,:) = pt_err_lecc;
stats1.pt_rt_ord_secc(ii,:,:) = pt_rt_secc;
stats1.pt_rt_ord_lecc(ii,:,:) = pt_rt_lecc;

stats1.a_err_ord_secc(ii,:,:) = a_err_secc;
stats1.a_err_ord_lecc(ii,:,:) = a_err_lecc;
stats1.a_rt_ord_secc(ii,:,:) = a_rt_secc;
stats1.a_rt_ord_lecc(ii,:,:) = a_rt_lecc;
% calculate systematic and unsystematic error for primary saccade (for any
% stimulus)
[stats1.p_sys_ord(ii,:),stats1.p_unsys_ord(ii,:)] = cal_sys_unsys(p_xend_pix_new,p_yend_pix_new,p_xtar_pix_new,p_ytar_pix_new,edf);
[stats1.p_sys_ord_secc(ii,:),stats1.p_unsys_ord_secc(ii,:)] = cal_sys_unsys(...
    p_xend_pix_new_secc,p_yend_pix_new_secc,p_xtar_pix_new_secc,p_ytar_pix_new_secc,edf);
[stats1.p_sys_ord_lecc(ii,:),stats1.p_unsys_ord_lecc(ii,:)] = cal_sys_unsys(...
    p_xend_pix_new_lecc,p_yend_pix_new_lecc,p_xtar_pix_new_lecc,p_ytar_pix_new_lecc,edf);
 
stats1.p_theta_err_ord(ii,:,:) = p_theta_err_new;
stats1.p_quad_err_ord(ii,:,:) = p_quad_err_new;
stats1.p_ord_err_ord(ii,:,:) = p_ord_err_new;

stats1.a_theta_err_ord(ii,:,:) = p_theta_err_new;
stats1.a_quad_err_ord(ii,:,:) = p_quad_err_new;
stats1.a_ord_err_ord(ii,:,:) = p_ord_err_new;

%%%%%%%%%%%% Organized by Quadrants %%%%%%%%%%%%%%%%%%%%%%%%%
for kk = 1:4 % for each quadrant
    for hh = 1:2 % for each task (cue type)
        % find trials where quadrant is kk, cue type is hh-1, and do not need to be revmoed
    ind = find(p_quad==kk & p_task==hh-1 & p_rm==1); 
    % primary saccade for any stimulus, error to target, N participant x 8 conditions
    % (serial order 1-4 for quadarnt cue, 1-4 for order cue)
p_err_new(:,kk+(hh-1)*4) = [p_err(ind);nan(16-length(ind),1)]; % at most 16 trials per serial position

% most accurate saccade error to target
a_err_new(:,kk+(hh-1)*4) = [a_err(ind);nan(16-length(ind),1)];
%n_sac_new(:,kk+(hh-1)*4) = [n_sac(ind);nan(16-length(ind),1)];

% closest saccade latency
a_rt_new(:,kk+(hh-1)*4) = [a_rt(ind);nan(16-length(ind),1)];

% primary saccade for any stimulus, latency
p_rt_new(:,kk+(hh-1)*4) = [p_rt(ind);nan(16-length(ind),1)];
%n_tr(1,kk+(hh-1)*4) = sum(~isnan(p_err(ind)));

% primary saccade for target, error to target and latency
pt_err_new(:,kk+(hh-1)*4) = [pt_err(ind);nan(16-length(ind),1)];
pt_rt_new(:,kk+(hh-1)*4) = [pt_rt(ind);nan(16-length(ind),1)];

% primary saccade endpoint x and y location in dva and pixel
p_xend_new(:,kk+(hh-1)*4) = [p_xend(ind);nan(16-length(ind),1)];
p_yend_new(:,kk+(hh-1)*4) = [p_yend(ind);nan(16-length(ind),1)];
p_xend_pix_new(:,kk+(hh-1)*4) = [p_xend_pix(ind);nan(16-length(ind),1)];
p_yend_pix_new(:,kk+(hh-1)*4) = [p_yend_pix(ind);nan(16-length(ind),1)];

% target x and y location in dva and pixel
p_xtar_new(:,kk+(hh-1)*4) = [p_xtar(ind);nan(16-length(ind),1)];
p_ytar_new(:,kk+(hh-1)*4) = [p_ytar(ind);nan(16-length(ind),1)];
p_xtar_pix_new(:,kk+(hh-1)*4) = [p_xtar_pix(ind);nan(16-length(ind),1)];
p_ytar_pix_new(:,kk+(hh-1)*4) = [p_ytar_pix(ind);nan(16-length(ind),1)];

% primary saccade (for any stimulus) angular, quadrant, and serial position error
p_theta_err_new(:,kk+(hh-1)*4) = [p_theta_err(ind);nan(16-length(ind),1)];
p_quad_err_new(:,kk+(hh-1)*4) = [p_quad_err(ind);nan(16-length(ind),1)];
p_ord_err_new(:,kk+(hh-1)*4) = [p_ord_err(ind);nan(16-length(ind),1)];

% most accurate saccade angular, quadrant, and serial position error
a_theta_err_new(:,kk+(hh-1)*4) = [a_theta_err(ind);nan(16-length(ind),1)];
a_quad_err_new(:,kk+(hh-1)*4) = [a_quad_err(ind);nan(16-length(ind),1)];
a_ord_err_new(:,kk+(hh-1)*4) = [a_ord_err(ind);nan(16-length(ind),1)];

% error by eccentricity
% trials with small eccentricity and large eccentricity (ind1 and ind2
% respectively)
ind1 = find(p_quad==kk & p_task==hh-1 & p_ecc==3 & p_rm == 1);
ind2 = find(p_quad==kk & p_task==hh-1 & p_ecc==5.5 & p_rm == 1);

% primary saccade (for any stimulus), error and latency, small and large eccentricity
p_err_secc(:,kk+(hh-1)*4) = [p_err(ind1);nan(8-length(ind1),1)]; % at most 8 trials per serial position per eccentricity per cue type
p_err_lecc(:,kk+(hh-1)*4) = [p_err(ind2);nan(8-length(ind2),1)];
p_rt_secc(:,kk+(hh-1)*4) = [p_rt(ind1);nan(8-length(ind1),1)];
p_rt_lecc(:,kk+(hh-1)*4) = [p_rt(ind2);nan(8-length(ind2),1)];

% primary saccade (for target), error and latency, small and large eccentricity
pt_err_secc(:,kk+(hh-1)*4) = [pt_err(ind1);nan(8-length(ind1),1)];
pt_err_lecc(:,kk+(hh-1)*4) = [pt_err(ind2);nan(8-length(ind2),1)];
pt_rt_secc(:,kk+(hh-1)*4) = [pt_rt(ind1);nan(8-length(ind1),1)];
pt_rt_lecc(:,kk+(hh-1)*4) = [pt_rt(ind2);nan(8-length(ind2),1)];

% most accurate saccade, error, small and large eccentricity
a_err_secc(:,kk+(hh-1)*4) = [a_err(ind1);nan(8-length(ind1),1)];
a_err_lecc(:,kk+(hh-1)*4) = [a_err(ind2);nan(8-length(ind2),1)];
a_rt_secc(:,kk+(hh-1)*4) = [a_rt(ind1);nan(8-length(ind1),1)];
a_rt_lecc(:,kk+(hh-1)*4) = [a_rt(ind2);nan(8-length(ind2),1)];

% primary saccade (for any stimulus), endpoint x and y location in pixel, small and
% large eccentricity
p_xend_pix_new_secc(:,kk+(hh-1)*4) = [p_xend_pix(ind1);nan(8-length(ind1),1)];
p_yend_pix_new_secc(:,kk+(hh-1)*4) = [p_yend_pix(ind1);nan(8-length(ind1),1)];
p_xend_pix_new_lecc(:,kk+(hh-1)*4) = [p_xend_pix(ind2);nan(8-length(ind2),1)];
p_yend_pix_new_lecc(:,kk+(hh-1)*4) = [p_yend_pix(ind2);nan(8-length(ind2),1)];

% target x and y location in pixel, small and large eccentricity
p_xtar_pix_new_secc(:,kk+(hh-1)*4) = [p_xtar_pix(ind1);nan(8-length(ind1),1)];
p_ytar_pix_new_secc(:,kk+(hh-1)*4) = [p_ytar_pix(ind1);nan(8-length(ind1),1)];
p_xtar_pix_new_lecc(:,kk+(hh-1)*4) = [p_xtar_pix(ind2);nan(8-length(ind2),1)];
p_ytar_pix_new_lecc(:,kk+(hh-1)*4) = [p_ytar_pix(ind2);nan(8-length(ind2),1)];
    end
end

% Store these variables in 'stats'
% Number of participant x Number of trials x Number of experimental
% conditions
% Here we have 8 conditions,
% serial order 1-4 for quadrant cue, and 1-4 for order cue
stats1.p_err_quad(ii,:,:) = p_err_new; 
stats1.a_err_quad(ii,:,:) = a_err_new;
stats1.a_rt_quad(ii,:,:) = a_rt_new;
%stats1.n_sac_quad(ii,:,:) = n_sac_new;
stats1.p_rt_quad(ii,:,:) = p_rt_new;
stats1.pt_err_quad(ii,:,:) = pt_err_new;
stats1.pt_rt_quad(ii,:,:) = pt_rt_new;
%stats1.n_tr_quad(ii,:,:) = n_tr;

stats1.p_xend_quad(ii,:,:) = p_xend_new;
stats1.p_yend_quad(ii,:,:) = p_yend_new;
stats1.p_xtar_quad(ii,:,:) = p_xtar_new;
stats1.p_ytar_quad(ii,:,:) = p_ytar_new;

stats1.p_xend_pix_quad(ii,:,:) = p_xend_pix_new;
stats1.p_yend_pix_quad(ii,:,:) = p_yend_pix_new;
stats1.p_xtar_pix_quad(ii,:,:) = p_xtar_pix_new;
stats1.p_ytar_pix_quad(ii,:,:) = p_ytar_pix_new;

stats1.p_xend_pix_quad_secc(ii,:,:) = p_xend_pix_new_secc;
stats1.p_xend_pix_quad_lecc(ii,:,:) = p_xend_pix_new_lecc;
stats1.p_yend_pix_quad_secc(ii,:,:) = p_yend_pix_new_secc;
stats1.p_yend_pix_quad_lecc(ii,:,:) = p_yend_pix_new_lecc;
stats1.p_xtar_pix_quad_secc(ii,:,:) = p_xtar_pix_new_secc;
stats1.p_xtar_pix_quad_lecc(ii,:,:) = p_xtar_pix_new_lecc;
stats1.p_ytar_pix_quad_secc(ii,:,:) = p_ytar_pix_new_secc;
stats1.p_ytar_pix_quad_lecc(ii,:,:) = p_ytar_pix_new_lecc;

stats1.p_err_quad_secc(ii,:,:) = p_err_secc;
stats1.p_err_quad_lecc(ii,:,:) = p_err_lecc;
stats1.p_rt_quad_secc(ii,:,:) = p_rt_secc;
stats1.p_rt_quad_lecc(ii,:,:) = p_rt_lecc;

stats1.pt_err_quad_secc(ii,:,:) = pt_err_secc;
stats1.pt_err_quad_lecc(ii,:,:) = pt_err_lecc;
stats1.pt_rt_quad_secc(ii,:,:) = pt_rt_secc;
stats1.pt_rt_quad_lecc(ii,:,:) = pt_rt_lecc;

stats1.a_err_quad_secc(ii,:,:) = a_err_secc;
stats1.a_err_quad_lecc(ii,:,:) = a_err_lecc;
stats1.a_rt_quad_secc(ii,:,:) = a_rt_secc;
stats1.a_rt_quad_lecc(ii,:,:) = a_rt_lecc;

% calculate systematic and unsystematic error for primary saccade (for any
% stimulus)
[stats1.p_sys_quad(ii,:),stats1.p_unsys_quad(ii,:)] = cal_sys_unsys(p_xend_pix_new,p_yend_pix_new,p_xtar_pix_new,p_ytar_pix_new,edf);
[stats1.p_sys_quad_secc(ii,:),stats1.p_unsys_quad_secc(ii,:)] = cal_sys_unsys(...
    p_xend_pix_new_secc,p_yend_pix_new_secc,p_xtar_pix_new_secc,p_ytar_pix_new_secc,edf);
[stats1.p_sys_quad_lecc(ii,:),stats1.p_unsys_quad_lecc(ii,:)] = cal_sys_unsys(...
    p_xend_pix_new_lecc,p_yend_pix_new_lecc,p_xtar_pix_new_lecc,p_ytar_pix_new_lecc,edf);
 
stats1.p_theta_err_quad(ii,:,:) = p_theta_err_new;
stats1.p_quad_err_quad(ii,:,:) = p_quad_err_new;
stats1.p_ord_err_quad(ii,:,:) = p_ord_err_new;

stats1.a_theta_err_quad(ii,:,:) = p_theta_err_new;
stats1.a_quad_err_quad(ii,:,:) = p_quad_err_new;
stats1.a_ord_err_quad(ii,:,:) = p_ord_err_new;

%%%%%%%%%%%% Organized by image type %%%%%%%%%%%%%%%%%%%%%%%%%
for kk = 1:4 % for each image type
    for hh = 1:2 % for each task (cue type)
        % find trials where quadrant is kk, cue type is hh-1, and do not need to be revmoed
    ind = find(p_img==kk & p_task==hh-1 & p_rm==1); 
    % primary saccade for any stimulus, error to target, N participant x 8 conditions
    % (serial order 1-4 for quadarnt cue, 1-4 for order cue)
p_err_new(:,kk+(hh-1)*4) = [p_err(ind);nan(16-length(ind),1)]; % at most 16 trials per serial position

% most accurate saccade error to target
a_err_new(:,kk+(hh-1)*4) = [a_err(ind);nan(16-length(ind),1)];
%n_sac_new(:,kk+(hh-1)*4) = [n_sac(ind);nan(16-length(ind),1)];

% closest saccade latency
a_rt_new(:,kk+(hh-1)*4) = [a_rt(ind);nan(16-length(ind),1)];

% primary saccade for any stimulus, latency
p_rt_new(:,kk+(hh-1)*4) = [p_rt(ind);nan(16-length(ind),1)];
%n_tr(1,kk+(hh-1)*4) = sum(~isnan(p_err(ind)));

% primary saccade for target, error to target and latency
pt_err_new(:,kk+(hh-1)*4) = [pt_err(ind);nan(16-length(ind),1)];
pt_rt_new(:,kk+(hh-1)*4) = [pt_rt(ind);nan(16-length(ind),1)];

% primary saccade endpoint x and y location in dva and pixel
p_xend_new(:,kk+(hh-1)*4) = [p_xend(ind);nan(16-length(ind),1)];
p_yend_new(:,kk+(hh-1)*4) = [p_yend(ind);nan(16-length(ind),1)];
p_xend_pix_new(:,kk+(hh-1)*4) = [p_xend_pix(ind);nan(16-length(ind),1)];
p_yend_pix_new(:,kk+(hh-1)*4) = [p_yend_pix(ind);nan(16-length(ind),1)];

% target x and y location in dva and pixel
p_xtar_new(:,kk+(hh-1)*4) = [p_xtar(ind);nan(16-length(ind),1)];
p_ytar_new(:,kk+(hh-1)*4) = [p_ytar(ind);nan(16-length(ind),1)];
p_xtar_pix_new(:,kk+(hh-1)*4) = [p_xtar_pix(ind);nan(16-length(ind),1)];
p_ytar_pix_new(:,kk+(hh-1)*4) = [p_ytar_pix(ind);nan(16-length(ind),1)];

% primary saccade (for any stimulus) angular, quadrant, and serial position error
p_theta_err_new(:,kk+(hh-1)*4) = [p_theta_err(ind);nan(16-length(ind),1)];
p_quad_err_new(:,kk+(hh-1)*4) = [p_quad_err(ind);nan(16-length(ind),1)];
p_ord_err_new(:,kk+(hh-1)*4) = [p_ord_err(ind);nan(16-length(ind),1)];

% most accurate saccade angular, quadrant, and serial position error
a_theta_err_new(:,kk+(hh-1)*4) = [a_theta_err(ind);nan(16-length(ind),1)];
a_quad_err_new(:,kk+(hh-1)*4) = [a_quad_err(ind);nan(16-length(ind),1)];
a_ord_err_new(:,kk+(hh-1)*4) = [a_ord_err(ind);nan(16-length(ind),1)];

% error by eccentricity
% trials with small eccentricity and large eccentricity (ind1 and ind2
% respectively)
ind1 = find(p_img==kk & p_task==hh-1 & p_ecc==3 & p_rm == 1);
ind2 = find(p_img==kk & p_task==hh-1 & p_ecc==5.5 & p_rm == 1);

% primary saccade (for any stimulus), error and latency, small and large eccentricity
p_err_secc(:,kk+(hh-1)*4) = [p_err(ind1);nan(8-length(ind1),1)]; % at most 8 trials per serial position per eccentricity per cue type
p_err_lecc(:,kk+(hh-1)*4) = [p_err(ind2);nan(8-length(ind2),1)];
p_rt_secc(:,kk+(hh-1)*4) = [p_rt(ind1);nan(8-length(ind1),1)];
p_rt_lecc(:,kk+(hh-1)*4) = [p_rt(ind2);nan(8-length(ind2),1)];

% primary saccade (for target), error and latency, small and large eccentricity
pt_err_secc(:,kk+(hh-1)*4) = [pt_err(ind1);nan(8-length(ind1),1)];
pt_err_lecc(:,kk+(hh-1)*4) = [pt_err(ind2);nan(8-length(ind2),1)];
pt_rt_secc(:,kk+(hh-1)*4) = [pt_rt(ind1);nan(8-length(ind1),1)];
pt_rt_lecc(:,kk+(hh-1)*4) = [pt_rt(ind2);nan(8-length(ind2),1)];

% most accurate saccade, error, small and large eccentricity
a_err_secc(:,kk+(hh-1)*4) = [a_err(ind1);nan(8-length(ind1),1)];
a_err_lecc(:,kk+(hh-1)*4) = [a_err(ind2);nan(8-length(ind2),1)];
a_rt_secc(:,kk+(hh-1)*4) = [a_rt(ind1);nan(8-length(ind1),1)];
a_rt_lecc(:,kk+(hh-1)*4) = [a_rt(ind2);nan(8-length(ind2),1)];

% primary saccade (for any stimulus), endpoint x and y location in pixel, small and
% large eccentricity
p_xend_pix_new_secc(:,kk+(hh-1)*4) = [p_xend_pix(ind1);nan(8-length(ind1),1)];
p_yend_pix_new_secc(:,kk+(hh-1)*4) = [p_yend_pix(ind1);nan(8-length(ind1),1)];
p_xend_pix_new_lecc(:,kk+(hh-1)*4) = [p_xend_pix(ind2);nan(8-length(ind2),1)];
p_yend_pix_new_lecc(:,kk+(hh-1)*4) = [p_yend_pix(ind2);nan(8-length(ind2),1)];

% target x and y location in pixel, small and large eccentricity
p_xtar_pix_new_secc(:,kk+(hh-1)*4) = [p_xtar_pix(ind1);nan(8-length(ind1),1)];
p_ytar_pix_new_secc(:,kk+(hh-1)*4) = [p_ytar_pix(ind1);nan(8-length(ind1),1)];
p_xtar_pix_new_lecc(:,kk+(hh-1)*4) = [p_xtar_pix(ind2);nan(8-length(ind2),1)];
p_ytar_pix_new_lecc(:,kk+(hh-1)*4) = [p_ytar_pix(ind2);nan(8-length(ind2),1)];
    end
end

% Store these variables in 'stats'
% Number of participant x Number of trials x Number of experimental
% conditions
% Here we have 8 conditions,
% serial order 1-4 for quadrant cue, and 1-4 for order cue
stats1.p_err_img(ii,:,:) = p_err_new; 
stats1.a_err_img(ii,:,:) = a_err_new;
stats1.a_rt_img(ii,:,:) = a_rt_new;

%stats1.n_sac_img(ii,:,:) = n_sac_new;
stats1.p_rt_img(ii,:,:) = p_rt_new;
stats1.pt_err_img(ii,:,:) = pt_err_new;
stats1.pt_rt_img(ii,:,:) = pt_rt_new;
%stats1.n_tr_img(ii,:,:) = n_tr;

stats1.p_xend_img(ii,:,:) = p_xend_new;
stats1.p_yend_img(ii,:,:) = p_yend_new;
stats1.p_xtar_img(ii,:,:) = p_xtar_new;
stats1.p_ytar_img(ii,:,:) = p_ytar_new;

stats1.p_xend_pix_img(ii,:,:) = p_xend_pix_new;
stats1.p_yend_pix_img(ii,:,:) = p_yend_pix_new;
stats1.p_xtar_pix_img(ii,:,:) = p_xtar_pix_new;
stats1.p_ytar_pix_img(ii,:,:) = p_ytar_pix_new;

stats1.p_xend_pix_img_secc(ii,:,:) = p_xend_pix_new_secc;
stats1.p_xend_pix_img_lecc(ii,:,:) = p_xend_pix_new_lecc;
stats1.p_yend_pix_img_secc(ii,:,:) = p_yend_pix_new_secc;
stats1.p_yend_pix_img_lecc(ii,:,:) = p_yend_pix_new_lecc;
stats1.p_xtar_pix_img_secc(ii,:,:) = p_xtar_pix_new_secc;
stats1.p_xtar_pix_img_lecc(ii,:,:) = p_xtar_pix_new_lecc;
stats1.p_ytar_pix_img_secc(ii,:,:) = p_ytar_pix_new_secc;
stats1.p_ytar_pix_img_lecc(ii,:,:) = p_ytar_pix_new_lecc;

stats1.p_err_img_secc(ii,:,:) = p_err_secc;
stats1.p_err_img_lecc(ii,:,:) = p_err_lecc;
stats1.p_rt_img_secc(ii,:,:) = p_rt_secc;
stats1.p_rt_img_lecc(ii,:,:) = p_rt_lecc;

stats1.pt_err_img_secc(ii,:,:) = pt_err_secc;
stats1.pt_err_img_lecc(ii,:,:) = pt_err_lecc;
stats1.pt_rt_img_secc(ii,:,:) = pt_rt_secc;
stats1.pt_rt_img_lecc(ii,:,:) = pt_rt_lecc;

stats1.a_err_img_secc(ii,:,:) = a_err_secc;
stats1.a_err_img_lecc(ii,:,:) = a_err_lecc;
stats1.a_rt_img_secc(ii,:,:) = a_rt_secc;
stats1.a_rt_img_lecc(ii,:,:) = a_rt_lecc;

% calculate systematic and unsystematic error for primary saccade (for any
% stimulus)
[stats1.p_sys_img(ii,:),stats1.p_unsys_img(ii,:)] = cal_sys_unsys(p_xend_pix_new,p_yend_pix_new,p_xtar_pix_new,p_ytar_pix_new,edf);
[stats1.p_sys_img_secc(ii,:),stats1.p_unsys_img_secc(ii,:)] = cal_sys_unsys(...
    p_xend_pix_new_secc,p_yend_pix_new_secc,p_xtar_pix_new_secc,p_ytar_pix_new_secc,edf);
[stats1.p_sys_img_lecc(ii,:),stats1.p_unsys_img_lecc(ii,:)] = cal_sys_unsys(...
    p_xend_pix_new_lecc,p_yend_pix_new_lecc,p_xtar_pix_new_lecc,p_ytar_pix_new_lecc,edf);
 
stats1.p_theta_err_img(ii,:,:) = p_theta_err_new;
stats1.p_quad_err_img(ii,:,:) = p_quad_err_new;
stats1.p_ord_err_img(ii,:,:) = p_ord_err_new;

stats1.a_theta_err_img(ii,:,:) = p_theta_err_new;
stats1.a_quad_err_img(ii,:,:) = p_quad_err_new;
stats1.a_ord_err_img(ii,:,:) = p_ord_err_new;
end

% store the subject id and run id
stats1.subs_id = subs_id1;
stats1.runs_id = {'serial position 1, quadrant cue',...
    'serial position 2, quadrant cue',...
    'serial position 3, quadrant cue',...
    'serial position 4, quadrant cue',...
    'serial position 1, order cue',...
    'serial position 2, order cue',...
    'serial position 3, order cue',...
    'serial position 4, order cue'};

% %% After further processing and quality check, we want to store the quality check results
% 
% qa.ntr_breakfix = ntr_breakfix;
% qa.ntr_err = ntr_err;
% qa.ntr_rm = ntr_rm;
% qa.ntr_nodc = ntr_nodc;
% %qa.ntr_err_nodc = ntr_err_nodc;
% qa.perc = perc;
% qa.nblink = nblink;
% qa.tr_breakfix = tr_breakfix;
% qa.tr_err = tr_err;
% qa.tr_rm = tr_rm;
% qa.tr_nodc = tr_nodc;
% 
% % RECALCULATE NUMBER OF TRIALS REMOVED
% tr_err_nodc = cellfun(@(x,y) unique([x;y']),qa.tr_err,qa.tr_nodc,'UniformOutput',false);
% ntr_err_nodc = cellfun(@(x) length(x),tr_err_nodc);
% qa.ntr_err_nodc = ntr_err_nodc;
% qa.tr_err_nodc = tr_err_nodc;
% 
% % Which participant to exclude?
% tr_exc = sum(qa.ntr_err_nodc,2,'omitnan')/128*100;
% mean(tr_exc,'omitnan')
% median(tr_exc,'omitnan')
% sum(tr_exc > 50)
% std(tr_exc,0,'omitnan')
% %final_sind2 = find(tr_exc <= 50);
% 
% tr_exc_quad = sum(qa.ntr_err_nodc(:,1:4),2,'omitnan')/64*100;
% tr_exc_ord = sum(qa.ntr_err_nodc(:,5:8),2,'omitnan')/64*100;
% sum(tr_exc_quad > 50 | tr_exc_ord > 50)
% final_sind1 = find(tr_exc_quad <= 50 & tr_exc_ord <= 50);
% 
% % plot number of trials excluded
% figure(1);clf
% histogram(tr_exc)
% 
% 
% % save the qa
% save(fullfile(output,'qa_exp1_longdelay.mat'),'qa','final_sind1')

%% Save the stats
% stats1.final_sind = final_sind1;
save(fullfile(output,'stats_exp1_longdelay.mat'),'stats1')

% 
% %%
% cd(output)
% save('data1_apr.mat', '-v7.3','data1','data_set1','subs_id1','subs_id2','runs_id1','runs_id2')
% 
% 
% %% Let's calculate data quality for v3_short
% clearvars -except data1 data2 data_dir1 data_dir2 data_set1 data_set2 file_dirs1 file_dirs2 subs_id1 subs_id2 runs_id1 runs_id2 subs1 subs2 runs1 runs2 output
% param = nan(length(subs_id1),length(runs_id1),16,3);
% for ii = 1:length(subs_id1)
%     for jj = 1:length(runs_id1)
%         if ~isempty(data1{ii,jj})
%         end
%     end
% end
% qa.perc_mean = mean(perc,2);
% qa.nblink_mean = mean(nblink,2)/16;
% qa.nnan_all = sum(ntr_rm,2);
% qa.perc_run = perc;
% qa.nblink_run = nblink;
% %qa.nnan_run = nnan;
% % qa.tr_breakfix = tr_breakfix;
% % qa.tr_err = tr_err;
% % qa.tr_rm = tr_rm;
% qa.ntr_breakfix = ntr_breakfix;
% qa.ntr_err = ntr_err;
% qa.ntr_rm = ntr_rm;
% qa.ntr_nodc = ntr_nodc;
% 
% 
% qa.nnan_order = nnan_order;
% qa.nnan_quad = nnan_quad;
% qa.nnan_img = nnan_img;
% 
% qa1 = qa;
% %save('qa_long_delay.mat',"qa1")
% tbl = table(subs_id1',qa.perc_mean,qa.nblink_mean,qa.nnan_all,'VariableNames',{'ID','% Sample Loss','Avg. No. Blinks Per Trial','No. Trials Excluded'});
% writetable(tbl,'data_quality_longdelay.csv')
% 
% %% Plot number of trials breaking fixation - v3short
% clearvars -except data1 data2 data_dir1 data_dir2 data_set1 data_set2 file_dirs1 file_dirs2 subs_id1 subs_id2 runs_id1 runs_id2 subs1 subs2 runs1 runs2 output qa1 qa2
% data_breakfix = qa1.ntr_breakfix; % perc of trials breaking fixation
% data_rmtrial = qa1.ntr_rm; % perc of trials removed
% data_nanerr = qa1.ntr_err; % perc of trials with error
% data_nodc = qa1.ntr_nodc; % perc of trials without drift correction
% 
% % figure(1);clf
% % set(gcf, 'Position',  [0, 0, 800, 600])
% % boxWithDots3(data_breakfix/16*100,{'Q1','Q2','Q3','Q4','O1','O2','O3','O4'},0,...
% %     [repmat([0 0.4470 0.7410],4,1);repmat([0.8500 0.3250 0.0980],4,1)],[],[]) % data_breakfix
% % xlabel('Run')
% % ylabel('Percentage of trials breaking fixation (%)')
% % ylim([0 80])
% % box off
% % saveas(gcf,'breakfix_long_delay.jpg')
% %
% % figure(2);clf
% % set(gcf, 'Position',  [0, 0, 800, 600])
% % boxWithDots3(data_nanerr/16*100,{'Q1','Q2','Q3','Q4','O1','O2','O3','O4'},0,...
% %     [repmat([0 0.4470 0.7410],4,1);repmat([0.8500 0.3250 0.0980],4,1)],[],[]) % data_breakfix
% % xlabel('Run')
% % ylabel('Percentage of trials with large errors (%)')
% % ylim([0 80])
% % box off
% % saveas(gcf,'nanerr_long_delay.jpg')
% %
% % figure(3);clf
% % set(gcf, 'Position',  [0, 0, 800, 600])
% % boxWithDots3(data_rmtrial/16*100,{'Q1','Q2','Q3','Q4','O1','O2','O3','O4'},0,...
% %     [repmat([0 0.4470 0.7410],4,1);repmat([0.8500 0.3250 0.0980],4,1)],[],[]) % data_breakfix
% % xlabel('Run')
% % ylabel('Percentage of trials removed from analysis (%)')
% % ylim([0 80])
% % box off
% % saveas(gcf,'rmtrial_long_delay.jpg')
% 
% % percentage of trials removed, across tasks
% quad_trial = sum(data_rmtrial(:,1:4),2)/64*100;
% quad_trial_breakfix = sum(data_breakfix(:,1:4),2)/64*100;
% quad_trial_nanerr = sum(data_nanerr(:,1:4),2)/64*100;
% quad_trial_nodc = sum(data_nodc(:,1:4),2)/64*100;
% 
% order_trial = sum(data_rmtrial(:,5:8),2)/64*100;
% order_trial_breakfix = sum(data_breakfix(:,5:8),2)/64*100;
% order_trial_nanerr = sum(data_nanerr(:,5:8),2)/64*100;
% order_trial_nodc = sum(data_nodc(:,5:8),2)/64*100;
% 
% [H,P,CI] = ttest(order_trial-quad_trial);
% figure(1);clf
% set(gcf, 'Position',  [0, 0, 800, 600])
% boxWithDots3([quad_trial order_trial],{'Quadrant Cue','Order Cue'},0, ...
%     [[0 0.4470 0.7410];[0.8500 0.3250 0.0980]],[P],{[1,2]})
% ylabel('Percentage of trials removed (%)')
% box off
% saveas(gcf,'rmtrial_task_long_delay.jpg')
% 
% [H,P,CI] = ttest(order_trial_breakfix-quad_trial_breakfix);
% figure(2);clf
% set(gcf, 'Position',  [0, 0, 800, 600])
% boxWithDots3([quad_trial_breakfix order_trial_breakfix],{'Quadrant Cue','Order Cue'},0, ...
%     [[0 0.4470 0.7410];[0.8500 0.3250 0.0980]],[P],{[1,2]})
% ylabel('Percentage of trials breaking fixation (%)')
% box off
% saveas(gcf,'breaktrial_task_long_delay.jpg')
% 
% [H,P,CI] = ttest(order_trial_nanerr-quad_trial_nanerr);
% figure(3);clf
% set(gcf, 'Position',  [0, 0, 800, 600])
% boxWithDots3([quad_trial_nanerr order_trial_nanerr],{'Quadrant Cue','Order Cue'},0, ...
%     [[0 0.4470 0.7410];[0.8500 0.3250 0.0980]],[P],{[1,2]})
% ylabel('Percentage of trials with large errors (%)')
% box off
% saveas(gcf,'errtrial_task_long_delay.jpg')
% 
% [H,P,CI] = ttest(order_trial_nodc-quad_trial_nodc);
% figure(4);clf
% set(gcf, 'Position',  [0, 0, 800, 600])
% boxWithDots3([quad_trial_nodc order_trial_nodc],{'Quadrant Cue','Order Cue'},0, ...
%     [[0 0.4470 0.7410];[0.8500 0.3250 0.0980]],[P],{[1,2]})
% ylabel('Percentage of trials without drift correction (%)')
% box off
% saveas(gcf,'nodctrial_task_long_delay.jpg')
% 
% %% Exclude participants
% % Exclude subjects based on number of trials excluded
% sub_rm = find(quad_trial > 50 | order_trial > 50 | isnan(sum(data_rmtrial,2)));
% % low quality data
% sum(quad_trial > 50 | order_trial > 50)
% % didn't complete
% sum(isnan(sum(data_rmtrial,2)))
% 
% data_rmtrial(sub_rm,:) = nan;
% data_breakfix(sub_rm,:) = nan;
% data_nanerr(sub_rm,:) = nan;
% 
% quad_trial = sum(data_rmtrial(:,1:4),2)/64*100;
% order_trial = sum(data_rmtrial(:,5:8),2)/64*100;
% [H,P,CI] = ttest(order_trial-quad_trial);
% 
% figure(5);clf
% set(gcf, 'Position',  [0, 0, 800, 600])
% boxWithDots3([quad_trial order_trial],{'Quadrant Cue','Order Cue'},0, ...
%     [[0 0.4470 0.7410];[0.8500 0.3250 0.0980]],[P],{[1,2]})
% ylabel('Percentage of trials removed (%)')
% box off
% saveas(gcf,'rmtrial_task_exc_long_delay.jpg')
% 
% % finally, which subjects are included in the analysis
% final_sind1 = find(~isnan(data_rmtrial(:,1)));
% 
% save('qa_long_delay.mat',"qa1","final_sind1")
% 
% %% Double check the bad participant
% bad_subj = sub_rm(4);
% data_nanerr(bad_subj,:)
% data_breakfix(bad_subj,:)
% data_nodc(bad_subj,:)
% subs_id1(bad_subj)
% 
% %miniEye_2024
