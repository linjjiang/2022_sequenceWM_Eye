function edf = select_saccades_exp1_exp2(edf,set)
% select specific saccade to analyze
% primary saccade:
% The first saccade that enters the AOI of any stimulus AND its latency is
% longer than 50 ms (not express saccade)

% most accurate saccade:
% The closest saccade to the target with a latency longer than 50 ms

% Changes I made on 4/16/2024:
% Clearly define bad trials vs. miss trials
% Bad trials: There is an artifact within the response window (1200 ms
% after cue onset)
% Miss trials: There is no response or a response that is too late,
% a response that is too small or too large (far away from any stimulus), 
% or a response to the non-target instead of the target
% Good trials: there is a response to the target within the designated
% response window and there is no artifact within the response window
% So primary saccade for all, e.g., primary_sac_* will contain all the
% trials including miss trials and good trials
% We will also have a separate variable to store the type of the trial:
% resp_type: it is a 16 x 1 cell, each cell contains the type of the
% responses indicated by the following:
% 0 - bad trial (with an artifact within the first 500 ms of the response window)
% 1 - no baseline during fixation
% 2 - miss trial - no response 
% 3 - miss trial - undefined
% 4 - early response (within 50 ms)
% 5 - large baseline during response (saccade start point >= 1.5 dva)
% 6 - target response
% 7 - non-target response
% 8 - outside of any region of interest
% 9 - others, undefined

% initiate variables
resp_type = 10*ones(edf.samples.ntrial,1);
primary_sac_tar_ind = nan(edf.samples.ntrial,1); % primary saccade to target, index
primary_sac_tar_err = nan(edf.samples.ntrial,1); % primary saccade to target, error to target
primary_sac_tar_err_all = nan(edf.samples.ntrial,4); % primary saccade to target, error to all stimuli
primary_sac_ind = nan(edf.samples.ntrial,1); % primary saccade to any (include misses and saccade to nontargets)
primary_sac_err = nan(edf.samples.ntrial,1); % primary saccade to any, error to target
primary_sac_err_all = nan(edf.samples.ntrial,4); % primary saccade to any, error to all stimuli
acc_sac_tar_ind = nan(edf.samples.ntrial,1); % most accurate saccade to target, index
acc_sac_tar_err = nan(edf.samples.ntrial,1); % most accurate saccade to target, error to target
acc_sac_tar_err_all = nan(edf.samples.ntrial,4); % most accurate saccade to target, error to all stimuli
acc_sac_tar_rt = nan(edf.samples.ntrial,1); % most accurate saccade latency
primary_sac_tm = nan(edf.samples.ntrial,1);
primary_sac_rt = nan(edf.samples.ntrial,1);
primary_sac_tar_tm = nan(edf.samples.ntrial,1);
primary_sac_tar_rt = nan(edf.samples.ntrial,1);

for ii = 1:edf.samples.ntrial
    % all the probe saccade
    ind = find(edf.events.sac_dc.trial == ii & (edf.events.sac_dc.msg_srt == set.sac.msg(1) | edf.events.sac_dc.msg_srt == set.sac.msg(2)));
    %          & edf.events.sac.amp >= 1 & edf.events.sac.dur >= 50 & ...
    
    % 1. Is there an artifact during the first 500 ms of the response window (0-500 ms)?
    % get artifact index
    idx_resp_window = edf.events.msg.ind_srt(ii,11) : ...
        edf.events.msg.ind_srt(ii,11) + double(500/(1000/edf.record.sample_rate)) - 1; % index of the response window
    % if no message in that trial -> others, undefined
    if all(isnan(idx_resp_window))
        % others, undefined
        continue;
    end
    % is there an artifact?
    artifact = edf.samples.is_artifact(idx_resp_window);
    if any(artifact) % yes
        resp_type(ii) = 0; % Skip the trial, assign the code 0 (artifact)
        continue;
    end
    
    % 2. Does this trial have a baseline?
    if ismember(ii,edf.trackloss.nodc_trial) % yes - if it is a member of no drift correction (baseline) trial
        resp_type(ii) = 1; % no baseline, assign the code 1 (no baseline)
        continue;
    end
    
    % 3. Is there a saccade during the response window (0-1500 ms)?
    if isempty(ind) % No saccade
        % 4. Are there fixations and all fixations are within 1.5 dva from the center cross?
        indfix = find(edf.samples.is_fixation == 1 & edf.samples.trial == ii & (edf.samples.msg == set.sac.msg(1) | edf.samples.msg == set.sac.msg(2)));
        if ~isempty(indfix) % there are fixations
            pos_fix = sqrt(edf.samples.x_deg_clean_dc(indfix,set.eye).^2 + edf.samples.y_deg_clean_dc(indfix,set.eye).^2);
            if all(pos_fix <= 1.5) % all fixation samples are within 1.5 dva from the center
                % This is a miss trial - no response
                resp_type(ii) = 2; % assign code 2 (miss trial - no response)
                continue;
            else % else, undefined
                
                resp_type(ii) = 3; % code 3 (miss trial - undefined)
                continue;
            end
        else % no fixations at all
            
            resp_type(ii) = 3; % code 3 (miss trial - undefined)
            continue;
        end
        
    else % there is a saccade detected
        % calculate saccade info
        % start position of saccades
        x_srt = edf.events.sac_dc.x_srt(ind);
        y_srt = edf.events.sac_dc.y_srt(ind);
        
        % end position of saccades
        x_end = edf.events.sac_dc.x_end(ind);
        y_end = edf.events.sac_dc.y_end(ind);
        
        % endpoint fixation position of saccades
        x_end_fix = edf.events.sac_dc.endfix_avg_x(ind)';
        y_end_fix = edf.events.sac_dc.endfix_avg_y(ind)';
        
        % target locations
        x_tar = edf.param.tarx_deg;
        y_tar = edf.param.tary_deg;
        
        % all stimuli locations
        x_stim = edf.param.stimx_deg;
        y_stim = edf.param.stimy_deg;
        
        % saccade amplitude
        amp = edf.events.sac_dc.amp(ind);
        
        % saccade onset
        sac_onset = edf.samples.time(edf.events.sac_dc.ind_srt(ind));
        
        % cue onset
        cue_onset = edf.events.msg.time(ii,11);
        
        % saccade latency
        sac_rt = sac_onset - cue_onset;
        
        % saccade endpoint euclidean distance to the stimuli (all stimuli)
        % rows - saccades
        % columns - always 4 columns, corresponding distance to 4 stimuli
        clear dist
        for tt = 1:4 % 4 items
            dist(:,tt) = sqrt((x_stim(ii,tt) - x_end_fix).^2 + (y_stim(ii,tt)- y_end_fix).^2);
        end
        
        % distance between the start position of the saccade and the fixation
        % cross
        clear dist_srt
        dist_srt = sqrt(x_srt.^2 + y_srt.^2);
        dist_srt = repmat(dist_srt,1,4);
        
        % which item (column) is the target
        col_tar = edf.param.probe_order(ii);
        
        % which item (column) is the nontarget
        col_ntar = 1:4; col_ntar(col_tar) = []; % column of nontarget
        
        % extract saccade endpoint distance (error) to the target
        err_to_tar = dist(:,col_tar);
        
        % extract saccade endpoint error to nontarget
        err_to_ntar = dist(:,col_ntar);
        
        % extract saccade endpoint distance to the fixation
        err_to_fix = sqrt((x_end_fix).^2 + (y_end_fix).^2);
        
        % now we've calculated all the saccade-related metrics. Let's start
        % figuring out which saccade to select!
        %     ind_large_amp = find(amp > 1);
        %
        %     if ~isempty(ind_large_amp)
        %     ind_first = ind_large_amp(1); % select the first saccade during the response window with an amplitude of at least 1 dva
        %     else % if all saccades are very small during the response window
        %         % assign code 10, small saccades
        %         resp_type(ii) = 10;
        %         % we select the whatever first saccade is there
        %         ind_first = 1;
        %     end
        ind_first = 1; % select the very first saccade
        % retrieve the first saccade endpoint fixation distance to all stimuli
        first_sac_end_stim = dist(ind_first,:);
        
        % Does the current first saccade have an endpoint fixation?
        if isnan(first_sac_end_stim(1)) % there is no endpoint fixation
            resp_type(ii) = 9; % no endpoint fixation
            continue;
        end
        
        % Is the first saccade endpoint outside of ROI of any stimulus?
        if all(first_sac_end_stim > 3.5)
            resp_type(ii) = 8; % outside
        end
        
        % Is the first saccade endpoint within 3.5 dva from the non-target?
        first_sac_end_ntar = err_to_ntar(ind_first,:);
        if any(first_sac_end_ntar <= 3.5) % yes
            % it is a non-target trial
            resp_type(ii) = 7;
        end
        
        % Is the first saccade endpoint within 3.5 dva from the target?
        first_sac_end_tar = err_to_tar(ind_first);
        if first_sac_end_tar <= 3.5 % yes
            % it is a good trial
            resp_type(ii) = 6; % target response
            % store the target saccade response
            % store the index, error to target, and error to target
            % stimulus
            primary_sac_tar_ind(ii,1) = ind(1); % this stores saccade to target
            primary_sac_tar_err(ii,1) = err_to_tar(1);
            primary_sac_tar_err_all(ii,:) = dist(1,:);
            primary_sac_tar_tm(ii,1) = sac_onset(1);
            primary_sac_tar_rt(ii,1) = sac_rt(1);
        end
        
        % Does the first saccade start point more than 2 dva from the center?
        first_sac_srt = dist_srt(ind_first,1);
        if first_sac_srt > 2 % yes
            % assign code 5 (large baseline)
            resp_type(ii) = 5;
            % we will record the saccade info later
        end
        
        % Does the first saccade starts within 50 ms of the cue onset (0-50 ms)?
        first_sac_rt = sac_rt(ind_first); % first saccade latency
        if first_sac_rt <= 50 % yes
            resp_type(ii) = 4; % assign code 4 - early response
            % we will record the saccade info later
        end
        
        % store the primary saccade information
        primary_sac_ind(ii,1) = ind(ind_first); % this stores saccade to any stimulus including target
        % so if we found a saccade to the target AOI, primary_sac_ind is
        % equivalent to primary_sac_tar_ind
        primary_sac_err(ii,1) = err_to_tar(ind_first); % error to target
        primary_sac_err_all(ii,:) = dist(ind_first,:); % distance to all 4 stimuli
        primary_sac_tm(ii,1) = sac_onset(ind_first);
        primary_sac_rt(ii,1) = sac_rt(ind_first);
        
        
        % we also want to retrieve the most accurate saccade (closest to the target)
        ind_a_potential = find(err_to_fix >= set.sac.aend_err & err_to_tar <= set.sac.end_err);
        % exclude this saccade if its endpoint is 1.5 dva within the fixation cross
        % and its endpoint is more than 3.5 dva from the target
        err_to_tar_left = err_to_tar(ind_a_potential); % error to the target of saccades with its endpoint position more than 1.5 dva from the center cross
        if ~isempty(err_to_tar_left) % if we could find such saccade
            [err,idx_a] = min(err_to_tar_left); % find one with the minimal error to the target
            acc_sac_tar_ind(ii,1) = ind(ind_a_potential(idx_a)'); % retrieve this saccade index
            acc_sac_tar_err(ii,1) = err; % error to target
            acc_sac_tar_err_all(ii,:) = dist(ind_a_potential(idx_a),:); % error to all stimuli
            acc_sac_tar_rt(ii,1) = sac_rt(ind_a_potential(idx_a));
        end
    end
end

% Store those values
edf.cal.resp_type = resp_type;
edf.cal.primary_sac_ind = primary_sac_ind;
edf.cal.primary_sac_err = primary_sac_err;
edf.cal.primary_sac_err_all = primary_sac_err_all;

edf.cal.primary_sac_tar_ind = primary_sac_tar_ind;
edf.cal.primary_sac_tar_err = primary_sac_tar_err;
edf.cal.primary_sac_tar_err_all = primary_sac_tar_err_all;

%edf.cal.num_sac = num_sac; % number of saccades per order
edf.cal.acc_sac_tar_ind = acc_sac_tar_ind;
edf.cal.acc_sac_tar_err = acc_sac_tar_err;
edf.cal.acc_sac_tar_err_all = acc_sac_tar_err_all;
edf.cal.acc_sac_tar_rt = acc_sac_tar_rt;

% number of trials to be excluded (primary_sac_err is nan)
edf.cal.num_nan = sum(isnan(primary_sac_err));

for ii = 1:length(primary_sac_ind)
    if ~isnan(primary_sac_ind(ii))
        edf.cal.primary_sac_xsrt(ii,1) = edf.events.sac_dc.x_srt(primary_sac_ind(ii));
        edf.cal.primary_sac_ysrt(ii,1) = edf.events.sac_dc.y_srt(primary_sac_ind(ii));
        edf.cal.primary_sac_xend(ii,1) = edf.events.sac_dc.endfix_avg_x(primary_sac_ind(ii));
        edf.cal.primary_sac_yend(ii,1) = edf.events.sac_dc.endfix_avg_y(primary_sac_ind(ii));
        
        edf.cal.primary_sac_xsrt_pix(ii,1) = edf.events.sac_dc.x_srt_pix(primary_sac_ind(ii));
        edf.cal.primary_sac_ysrt_pix(ii,1) = edf.events.sac_dc.y_srt_pix(primary_sac_ind(ii));
        edf.cal.primary_sac_xend_pix(ii,1) = edf.events.sac_dc.endfix_avg_x_pix(primary_sac_ind(ii));
        edf.cal.primary_sac_yend_pix(ii,1) = edf.events.sac_dc.endfix_avg_y_pix(primary_sac_ind(ii));
        
    else
        edf.cal.primary_sac_xsrt(ii,1) = nan;
        edf.cal.primary_sac_ysrt(ii,1) = nan;
        edf.cal.primary_sac_xend(ii,1) = nan;
        edf.cal.primary_sac_yend(ii,1) = nan;
        
        edf.cal.primary_sac_xsrt_pix(ii,1) = nan;
        edf.cal.primary_sac_ysrt_pix(ii,1) = nan;
        edf.cal.primary_sac_xend_pix(ii,1) = nan;
        edf.cal.primary_sac_yend_pix(ii,1) = nan;
    end
    
    if ~isnan(primary_sac_tar_ind(ii))
        edf.cal.primary_sac_tar_xsrt(ii,1) = edf.events.sac_dc.x_srt(primary_sac_tar_ind(ii));
        edf.cal.primary_sac_tar_ysrt(ii,1) = edf.events.sac_dc.y_srt(primary_sac_tar_ind(ii));
        edf.cal.primary_sac_tar_xend(ii,1) = edf.events.sac_dc.endfix_avg_x(primary_sac_tar_ind(ii));
        edf.cal.primary_sac_tar_yend(ii,1) = edf.events.sac_dc.endfix_avg_y(primary_sac_tar_ind(ii));
        
        edf.cal.primary_sac_tar_xsrt_pix(ii,1) = edf.events.sac_dc.x_srt_pix(primary_sac_tar_ind(ii));
        edf.cal.primary_sac_tar_ysrt_pix(ii,1) = edf.events.sac_dc.y_srt_pix(primary_sac_tar_ind(ii));
        edf.cal.primary_sac_tar_xend_pix(ii,1) = edf.events.sac_dc.endfix_avg_x_pix(primary_sac_tar_ind(ii));
        edf.cal.primary_sac_tar_yend_pix(ii,1) = edf.events.sac_dc.endfix_avg_y_pix(primary_sac_tar_ind(ii));
        
    else
        edf.cal.primary_sac_tar_xsrt(ii,1) = nan;
        edf.cal.primary_sac_tar_ysrt(ii,1) = nan;
        edf.cal.primary_sac_tar_xend(ii,1) = nan;
        edf.cal.primary_sac_tar_yend(ii,1) = nan;
        
        edf.cal.primary_sac_tar_xsrt_pix(ii,1) = nan;
        edf.cal.primary_sac_tar_ysrt_pix(ii,1) = nan;
        edf.cal.primary_sac_tar_xend_pix(ii,1) = nan;
        edf.cal.primary_sac_tar_yend_pix(ii,1) = nan;
    end
    
    if ~isnan(acc_sac_tar_ind(ii))
        edf.cal.acc_sac_tar_xsrt(ii,1) = edf.events.sac_dc.x_srt(acc_sac_tar_ind(ii));
        edf.cal.acc_sac_tar_ysrt(ii,1) = edf.events.sac_dc.y_srt(acc_sac_tar_ind(ii));
        edf.cal.acc_sac_tar_xend(ii,1) = edf.events.sac_dc.endfix_avg_x(acc_sac_tar_ind(ii));
        edf.cal.acc_sac_tar_yend(ii,1) = edf.events.sac_dc.endfix_avg_y(acc_sac_tar_ind(ii));
        
        edf.cal.acc_sac_tar_xsrt_pix(ii,1) = edf.events.sac_dc.x_srt_pix(acc_sac_tar_ind(ii));
        edf.cal.acc_sac_tar_ysrt_pix(ii,1) = edf.events.sac_dc.y_srt_pix(acc_sac_tar_ind(ii));
        edf.cal.acc_sac_tar_xend_pix(ii,1) = edf.events.sac_dc.endfix_avg_x_pix(acc_sac_tar_ind(ii));
        edf.cal.acc_sac_tar_yend_pix(ii,1) = edf.events.sac_dc.endfix_avg_y_pix(acc_sac_tar_ind(ii));
        
    else
        edf.cal.acc_sac_tar_xsrt(ii,1) = nan;
        edf.cal.acc_sac_tar_ysrt(ii,1) = nan;
        edf.cal.acc_sac_tar_xend(ii,1) = nan;
        edf.cal.acc_sac_tar_yend(ii,1) = nan;
        
        edf.cal.acc_sac_tar_xsrt_pix(ii,1) = nan;
        edf.cal.acc_sac_tar_ysrt_pix(ii,1) = nan;
        edf.cal.acc_sac_tar_xend_pix(ii,1) = nan;
        edf.cal.acc_sac_tar_yend_pix(ii,1) = nan;
    end
end

edf.cal.primary_sac_tm = primary_sac_tm;
edf.cal.primary_sac_rt = primary_sac_rt;
edf.cal.primary_sac_tar_tm = primary_sac_tar_tm;
edf.cal.primary_sac_tar_rt = primary_sac_tar_rt;

end